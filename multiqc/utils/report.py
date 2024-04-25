"""MultiQC report module. Holds the output from each
module. Is available to subsequent modules. Contains
helper functions to generate markup for report."""

import fnmatch
import inspect
import io
import json
import logging
import mimetypes
import os
import re
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union, List, Optional, Any

import math
import rich
import rich.progress
import yaml
from tqdm import tqdm

from multiqc.utils import lzstring
from multiqc.utils import log

from . import config
from .util_functions import replace_defaultdicts, is_running_in_notebook, no_unicode

logger = logging.getLogger(__name__)

initialized = False

tmp_dir: Optional[str] = None

# Uninitialised global variables for static typing
multiqc_command: str
modules_output: List  # List of BaseMultiqcModule objects
general_stats_html: str
lint_errors: List[str]
num_flat_plots: int
saved_raw_data: Dict[str, Dict[str, Any]]  # Indexed by unique key, then sample name
last_found_file: Optional[str]
runtimes: Dict[str, Union[float, Dict]]
file_search_stats: Dict[str, int]
searchfiles: List
files: Dict

# Fields below is kept between interactive runs
data_sources: Dict[str, Dict[str, Dict]]
html_ids: List[str]
plot_data: Dict = dict()
general_stats_data: List[Dict]
general_stats_headers: List[Dict]
# Map of Software tools to a set of unique version strings
software_versions: Dict[str, Dict[str, List]]


def init():
    # Set up global variables shared across modules. Inside a function so that the global
    # vars are reset if MultiQC is run more than once within a single session / environment.
    global initialized
    global tmp_dir
    global multiqc_command
    global modules_output
    global general_stats_html
    global lint_errors
    global num_flat_plots
    global saved_raw_data
    global last_found_file
    global runtimes
    global data_sources
    global html_ids
    global plot_data
    global general_stats_data
    global general_stats_headers
    global software_versions

    # Create new temporary directory for module data exports
    initialized = True
    tmp_dir = tempfile.mkdtemp()
    logger.debug(f"Using temporary directory: {tmp_dir}")
    multiqc_command = ""
    modules_output = []
    general_stats_html = ""
    lint_errors = []
    num_flat_plots = 0
    saved_raw_data = dict()
    last_found_file = None
    runtimes = {
        "total": 0.0,
        "total_sp": 0.0,
        "total_mods": 0.0,
        "total_compression": 0.0,
        "sp": defaultdict(),
        "mods": defaultdict(),
    }
    data_sources = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    html_ids = []
    plot_data = dict()
    general_stats_data = []
    general_stats_headers = []
    software_versions = defaultdict(lambda: defaultdict(list))

    reset_file_search()


def reset():
    """
    Reset interactive session.
    """
    init()


def reset_file_search():
    """
    Reset the file search session. Useful in interactive session to call in the end
    of parse_logs(), so the next run is not affected by the previous one.
    """
    global searchfiles
    global files
    global file_search_stats
    searchfiles = list()
    files = dict()  # Discovered files for each search key
    file_search_stats = {
        "skipped_symlinks": 0,
        "skipped_not_a_file": 0,
        "skipped_ignore_pattern": 0,
        "skipped_filesize_limit": 0,
        "skipped_module_specific_max_filesize": 0,
        "skipped_no_match": 0,
        "skipped_directory_fn_ignore_dirs": 0,
        "skipped_file_contents_search_errors": 0,
    }


def is_searching_in_source_dir(path: Path) -> bool:
    """
    Checks whether MultiQC is searching for files in the source code folder
    """
    multiqc_installation_dir_files = [
        "LICENSE",
        "CHANGELOG.md",
        "Dockerfile",
        "MANIFEST.in",
        ".gitmodules",
        "README.md",
        "setup.py",
        ".gitignore",
    ]

    filenames = [f.name for f in path.iterdir() if f.is_file()]

    if len(filenames) > 0 and all([fn in filenames for fn in multiqc_installation_dir_files]):
        logger.error(f"Error: MultiQC is running in source code directory! {path}")
        logger.warning("Please see the docs for how to use MultiQC: https://multiqc.info/docs/#running-multiqc")
        return True
    else:
        return False


def handle_analysis_path(item: Path):
    """
    Branching logic to handle analysis paths (directories and files)
    Walks directory trees recursively calling pathlib's `path.iterdir()`.
    Guaranteed to work correctly with symlinks even on non-POSIX compliant filesystems.
    """
    if item.is_symlink() and config.ignore_symlinks:
        file_search_stats["skipped_symlinks"] += 1
        return
    elif item.is_file():
        searchfiles.append([item.name, os.fspath(item.parent)])
    elif item.is_dir():
        # Skip directory if it matches ignore patterns
        d_matches = any(d for d in config.fn_ignore_dirs if item.match(d.rstrip(os.sep)))
        p_matches = any(p for p in config.fn_ignore_paths if item.match(p.rstrip(os.sep)))
        if d_matches or p_matches:
            file_search_stats["skipped_directory_fn_ignore_dirs"] += 1
            return

        # Check not running in install directory
        if is_searching_in_source_dir(item):
            return

        for item in item.iterdir():
            handle_analysis_path(item)


def get_filelist(run_module_names):
    """
    Go through all supplied search directories and assembly a master
    list of files to search. Then fire search functions for each file.
    """
    # Prep search patterns
    spatterns = [{}, {}, {}, {}, {}, {}, {}]
    runtimes["sp"] = defaultdict()
    ignored_patterns = []
    skipped_patterns = []
    for key, sps in config.sp.items():
        mod_name = key.split("/", 1)[0]
        if mod_name.lower() not in [m.lower() for m in run_module_names]:
            ignored_patterns.append(key)
            continue
        files[key] = list()
        if not isinstance(sps, list):
            sps = [sps]

        # Warn if we have any unrecognised search pattern keys
        expected_sp_keys = [
            "fn",
            "fn_re",
            "contents",
            "contents_re",
            "num_lines",
            "shared",
            "skip",
            "max_filesize",
            "exclude_fn",
            "exclude_fn_re",
            "exclude_contents",
            "exclude_contents_re",
        ]
        unrecognised_keys = [y for x in sps for y in x.keys() if y not in expected_sp_keys]
        if len(unrecognised_keys) > 0:
            logger.warning(f"Unrecognised search pattern keys for '{key}': {', '.join(unrecognised_keys)}")

        # Check if we are skipping this search key
        if any([x.get("skip") for x in sps]):
            skipped_patterns.append(key)

        # Split search patterns according to speed of execution.
        if any([x for x in sps if "contents_re" in x]):
            if any([x for x in sps if "num_lines" in x]):
                spatterns[4][key] = sps
            elif any([x for x in sps if "max_filesize" in x]):
                spatterns[5][key] = sps
            else:
                spatterns[6][key] = sps
        elif any([x for x in sps if "contents" in x]):
            if any([x for x in sps if "num_lines" in x]):
                spatterns[1][key] = sps
            elif any([x for x in sps if "max_filesize" in x]):
                spatterns[2][key] = sps
            else:
                spatterns[3][key] = sps
        else:
            spatterns[0][key] = sps

    if len(ignored_patterns) > 0:
        logger.debug(f"Ignored {len(ignored_patterns)} search patterns as didn't match running modules.")

    if len(skipped_patterns) > 0:
        logger.info(f"Skipping {len(skipped_patterns)} file search patterns")
        logger.debug(f"Skipping search patterns: {', '.join(skipped_patterns)}")

    def add_file(fn, root):
        """
        Function applied to each file found when walking the analysis
        directories. Runs through all search patterns and returns True
        if a match is found.
        """
        f = {"fn": fn, "root": root}

        # Check that this is a file and not a pipe or anything weird
        if not os.path.isfile(os.path.join(root, fn)):
            file_search_stats["skipped_not_a_file"] += 1
            return False

        # Check that we don't want to ignore this file
        i_matches = [n for n in config.fn_ignore_files if fnmatch.fnmatch(fn, n)]
        if len(i_matches) > 0:
            file_search_stats["skipped_ignore_pattern"] += 1
            return False

        # Limit search to small files, to avoid 30GB FastQ files etc.
        try:
            f["filesize"] = os.path.getsize(os.path.join(root, fn))
        except (IOError, OSError, ValueError, UnicodeDecodeError):
            logger.debug(f"Couldn't read file when checking filesize: {fn}")
        else:
            if f["filesize"] > config.log_filesize_limit:
                file_search_stats["skipped_filesize_limit"] += 1
                return False

        # Use mimetypes to exclude binary files where possible
        if not re.match(r".+_mqc\.(png|jpg|jpeg)", f["fn"]) and config.ignore_images:
            (ftype, encoding) = mimetypes.guess_type(os.path.join(f["root"], f["fn"]))
            if encoding is not None:
                return False
            if ftype is not None and ftype.startswith("image"):
                return False

        # Test file for each search pattern
        file_matched = False
        for patterns in spatterns:
            for key, sps in patterns.items():
                start = time.time()
                for sp in sps:
                    if search_file(sp, f, key):
                        # Check that we shouldn't exclude this file
                        if not exclude_file(sp, f):
                            # Looks good! Remember this file
                            files[key].append(f)
                            file_search_stats[key] = file_search_stats.get(key, 0) + 1
                            file_matched = True
                        # Don't keep searching this file for other modules
                        if not sp.get("shared", False) and key not in config.filesearch_file_shared:
                            runtimes["sp"][key] = runtimes["sp"].get(key, 0) + (time.time() - start)
                            return True
                        # Don't look at other patterns for this module
                        break
                runtimes["sp"][key] = runtimes["sp"].get(key, 0) + (time.time() - start)

        return file_matched

    # Go through the analysis directories and get file list
    total_sp_starttime = time.time()
    for path in config.analysis_dir:
        handle_analysis_path(Path(path))

    if is_running_in_notebook() or no_unicode():
        PRINT_FNAME = False
        # ANSI escape code for dim text
        if not config.no_ansi:
            DIM = "\033[2m"
            BLUE = "\033[34m"
            RESET = "\033[0m"
        else:
            DIM = ""
            BLUE = ""
            RESET = ""

        bar_format = f"{BLUE}| {'searching':>17} {RESET}| " + "{bar:40} {percentage:3.0f}% {n_fmt}/{total_fmt}"
        if PRINT_FNAME:
            bar_format += " {desc}"

        # Set up the tqdm progress bar
        with tqdm(
            total=len(searchfiles),
            desc="Searching",
            unit="file",
            file=sys.stdout,
            disable=config.no_ansi or config.quiet,
            bar_format=bar_format,
        ) as pbar:
            for i, sf in enumerate(searchfiles):
                # Update the progress bar description with the file being searched
                if PRINT_FNAME and i % 100 == 0:
                    pbar.set_description_str(f"{DIM}{os.path.join(sf[1], sf[0])[-50:]}{RESET}")
                pbar.update(1)

                # Your file processing logic
                if not add_file(sf[0], sf[1]):
                    file_search_stats["skipped_no_match"] += 1

            # Clear the description after the loop is complete
            pbar.set_description("")
            pbar.refresh()
    else:
        progress_obj = rich.progress.Progress(
            "[blue]|[/]      ",
            rich.progress.SpinnerColumn(),
            "[blue]{task.description}[/] |",
            rich.progress.BarColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            "[green]{task.completed}/{task.total}",
            # "[dim]{task.fields[s_fn]}[/]",
            console=log.rich_console,
            disable=config.no_ansi or config.quiet,
        )
        with progress_obj as progress:
            mqc_task = progress.add_task("searching", total=len(searchfiles), s_fn="")
            for sf in searchfiles:
                progress.update(mqc_task, advance=1, s_fn=os.path.join(sf[1], sf[0])[-50:])
                if not add_file(sf[0], sf[1]):
                    file_search_stats["skipped_no_match"] += 1
            progress.update(mqc_task, s_fn="")

    runtimes["total_sp"] = time.time() - total_sp_starttime
    if config.profile_runtime:
        logger.info(f"Profile-runtime: Searching files took {runtimes['total_sp']:.2f}s")

    # Debug log summary about what we skipped
    summaries = []
    for key in sorted(file_search_stats, key=file_search_stats.get, reverse=True):
        if "skipped_" in key and file_search_stats[key] > 0:
            summaries.append(f"{key}: {file_search_stats[key]}")
    logger.debug(f"Summary of files that were skipped by the search: |{'|, |'.join(summaries)}|")


def search_file(pattern, f, module_key):
    """
    Function to searach a single file for a single search pattern.
    """

    global file_search_stats
    fn_matched = False
    contents_matched = False

    # Search pattern specific filesize limit
    if pattern.get("max_filesize") is not None and "filesize" in f:
        if f["filesize"] > pattern.get("max_filesize"):
            file_search_stats["skipped_module_specific_max_filesize"] += 1
            return False

    # Search by file name (glob)
    if pattern.get("fn") is not None:
        if fnmatch.fnmatch(f["fn"], pattern["fn"]):
            fn_matched = True
            if pattern.get("contents") is None and pattern.get("contents_re") is None:
                return True
        else:
            return False

    # Search by file name (regex)
    if pattern.get("fn_re") is not None:
        if re.match(pattern["fn_re"], f["fn"]):
            fn_matched = True
            if pattern.get("contents") is None and pattern.get("contents_re") is None:
                return True
        else:
            return False

    # Search by file contents
    repattern = None
    if pattern.get("contents") is not None or pattern.get("contents_re") is not None:
        if pattern.get("contents_re") is not None:
            repattern = re.compile(pattern["contents_re"])
        if "contents_lines" not in f or ("num_lines" in pattern and len(f["contents_lines"]) < pattern["num_lines"]):
            f["contents_lines"] = []
            file_path = os.path.join(f["root"], f["fn"])

            try:
                fh = io.open(file_path, "r", encoding="utf-8")
            except Exception as e:
                if config.report_readerrors:
                    logger.debug(f"Couldn't read file when looking for output: {file_path}, {e}")
                file_search_stats["skipped_file_contents_search_errors"] += 1
                return False
            else:
                try:
                    for i, line in enumerate(fh):
                        f["contents_lines"].append(line)
                        if i >= config.filesearch_lines_limit and i >= pattern.get("num_lines", 0):
                            break
                except UnicodeDecodeError as e:
                    if config.report_readerrors:
                        logger.debug(
                            f"Couldn't read file as a utf-8 text when looking for output: {file_path}, {e}. "
                            f"Usually because it's a binary file. But sometimes there are single non-unicode "
                            f"characters, so attempting reading while skipping such characters."
                        )
                    try:
                        with io.open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
                            for i, line in enumerate(fh):
                                f["contents_lines"].append(line)
                                if i >= config.filesearch_lines_limit and i >= pattern.get("num_lines", 0):
                                    break
                    except Exception as e:
                        if config.report_readerrors:
                            logger.debug(f"Still couldn't read the file, skipping: {file_path}, {e}")
                        file_search_stats["skipped_file_contents_search_errors"] += 1
                        return False
                    else:
                        if not f["contents_lines"]:
                            if config.report_readerrors:
                                logger.debug(f"No utf-8 lines were read from the file, skipping {file_path}")
                            file_search_stats["skipped_file_contents_search_errors"] += 1
                            return False
                except Exception as e:
                    if config.report_readerrors:
                        logger.debug(f"Couldn't read file when looking for output: {file_path}, {e}")
                    file_search_stats["skipped_file_contents_search_errors"] += 1
                    return False
            finally:
                fh.close()

        # Go through the parsed file contents
        for i, line in enumerate(f["contents_lines"]):
            # Break if we've searched enough lines for this pattern
            if pattern.get("num_lines") and i >= pattern.get("num_lines"):
                break

            # Search by file contents (string)
            if pattern.get("contents") is not None:
                if pattern["contents"] in line:
                    contents_matched = True
                    if pattern.get("fn") is None and pattern.get("fn_re") is None:
                        return True
                    break
            # Search by file contents (regex)
            elif pattern.get("contents_re") is not None:
                if re.search(repattern, line):
                    contents_matched = True
                    if pattern.get("fn") is None and pattern.get("fn_re") is None:
                        return True
                    break

    return fn_matched and contents_matched


def exclude_file(sp, f):
    """
    Exclude discovered files if they match the special exclude_
    search pattern keys
    """
    # Make everything a list if it isn't already
    for k in sp:
        if k in ["exclude_fn", "exclude_fn_re" "exclude_contents", "exclude_contents_re"]:
            if not isinstance(sp[k], list):
                sp[k] = [sp[k]]

    # Search by file name (glob)
    if "exclude_fn" in sp:
        for pat in sp["exclude_fn"]:
            if fnmatch.fnmatch(f["fn"], pat):
                return True

    # Search by file name (regex)
    if "exclude_fn_re" in sp:
        for pat in sp["exclude_fn_re"]:
            if re.match(pat, f["fn"]):
                return True

    # Search the contents of the file
    if "exclude_contents" in sp or "exclude_contents_re" in sp:
        # Compile regex patterns if we have any
        if "exclude_contents_re" in sp:
            sp["exclude_contents_re"] = [re.compile(pat) for pat in sp["exclude_contents_re"]]
        with io.open(os.path.join(f["root"], f["fn"]), "r", encoding="utf-8") as fh:
            for line in fh:
                if "exclude_contents" in sp:
                    for pat in sp["exclude_contents"]:
                        if pat in line:
                            return True
                if "exclude_contents_re" in sp:
                    for pat in sp["exclude_contents_re"]:
                        if re.search(pat, line):
                            return True
    return False


def data_sources_tofile():
    fn = f"multiqc_sources.{config.data_format_extensions[config.data_format]}"
    with io.open(os.path.join(config.data_dir, fn), "w", encoding="utf-8") as f:
        if config.data_format == "json":
            jsonstr = json.dumps(data_sources, indent=4, ensure_ascii=False)
            print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
        elif config.data_format == "yaml":
            yaml.dump(replace_defaultdicts(data_sources), f, default_flow_style=False)
        else:
            lines = [["Module", "Section", "Sample Name", "Source"]]
            for mod in data_sources:
                for sec in data_sources[mod]:
                    for s_name, source in data_sources[mod][sec].items():
                        lines.append([mod, sec, s_name, source])
            body = "\n".join(["\t".join(line) for line in lines])
            print(body.encode("utf-8", "ignore").decode("utf-8"), file=f)


def dois_tofile(modules_output):
    """Find all DOIs listed in report sections and write to a file"""
    # Collect DOIs
    dois = {"MultiQC": ["10.1093/bioinformatics/btw354"]}
    for mod in modules_output:
        if mod.doi is not None and mod.doi != "" and mod.doi != []:
            dois[mod.anchor] = mod.doi
    # Write to a file
    fn = f"multiqc_citations.{config.data_format_extensions[config.data_format]}"
    with io.open(os.path.join(config.data_dir, fn), "w", encoding="utf-8") as f:
        if config.data_format == "json":
            jsonstr = json.dumps(dois, indent=4, ensure_ascii=False)
            print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
        elif config.data_format == "yaml":
            yaml.dump(replace_defaultdicts(dois), f, default_flow_style=False)
        else:
            body = ""
            for mod, dois in dois.items():
                for doi in dois:
                    body += f"{doi}{' ' * (50 - len(doi))} # {mod}\n"
            print(body.encode("utf-8", "ignore").decode("utf-8"), file=f)


def save_htmlid(html_id, skiplint=False):
    """Take a HTML ID, sanitise for HTML, check for duplicates and save.
    Returns sanitised, unique ID"""
    global html_ids
    global lint_errors

    # Trailing whitespace
    html_id_clean = html_id.strip()

    # Trailing underscores
    html_id_clean = html_id_clean.strip("_")

    # Must begin with a letter
    if re.match(r"^[a-zA-Z]", html_id_clean) is None:
        html_id_clean = f"mqc_{html_id_clean}"

    # Replace illegal characters
    html_id_clean = re.sub("[^a-zA-Z0-9_-]+", "_", html_id_clean)

    # Validate if linting
    modname = ""
    codeline = ""
    if config.strict and not skiplint:
        callstack = inspect.stack()
        for n in callstack:
            if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
                callpath = n[1].split("multiqc/modules/", 1)[-1]
                modname = f">{callpath}< "
                codeline = n[4][0].strip()
                break
        if html_id != html_id_clean:
            errmsg = f"LINT: {modname}HTML ID was not clean ('{html_id}' -> '{html_id_clean}') ## {codeline}"
            logger.error(errmsg)
            lint_errors.append(errmsg)

    # Check for duplicates
    i = 1
    html_id_base = html_id_clean
    while html_id_clean in html_ids:
        if html_id_clean == "general_stats_table":
            raise ValueError("HTML ID 'general_stats_table' is reserved and cannot be used")
        html_id_clean = f"{html_id_base}-{i}"
        i += 1
        if config.strict and not skiplint:
            errmsg = f"LINT: {modname}HTML ID was a duplicate ({html_id_clean}) ## {codeline}"
            logger.error(errmsg)
            lint_errors.append(errmsg)

    # Remember and return
    html_ids.append(html_id_clean)
    return html_id_clean


def compress_json(data):
    """Take a Python data object. Convert to JSON and compress using lzstring"""
    json_string = json.dumps(data).encode("utf-8", "ignore").decode("utf-8")
    json_string = sanitise_json(json_string)
    x = lzstring.LZString()
    return x.compressToBase64(json_string)


def sanitise_json(json_string):
    """
    The Python json module uses a bunch of values which are valid JavaScript
    but invalid JSON. These crash the browser when parsing the JSON.
    Nothing in the MultiQC front-end uses these values, so instead we just
    do a find-and-replace for them and switch them with `null`, which works fine.

    Side effect: Any string values that include the word "Infinity"
    (case-sensitive) will have it switched for "null". Hopefully that doesn't happen
    a lot, otherwise we'll have to do this in a more complicated manner.
    """
    json_string = re.sub(r"\bNaN\b", "null", json_string)
    json_string = re.sub(r"\b-?Infinity\b", "null", json_string)
    return json_string


def data_tmp_dir() -> str:
    """
    Temporary directory to collect data files from running modules before copying to the final
    destination in multiqc.core.write_results
    """
    path = os.path.join(tmp_dir, "multiqc_data")
    os.makedirs(path, exist_ok=True)
    return path


def plots_tmp_dir() -> str:
    """
    Temporary directory to collect plot exports from running modules before copying to the final
    destination in multiqc.core.write_results
    """
    path = os.path.join(tmp_dir, "multiqc_plots")
    os.makedirs(path, exist_ok=True)
    return path


def write_data_file(
    data: Union[Dict[str, Union[Dict, List]], List[Dict]],
    fn: str,
    sort_cols=False,
    data_format=None,
):
    """
    Write a data file to the report directory. Will do nothing
    if config.data_dir is not set (i.e. when config.make_data_dir == False or config.filename == "stdout")
    :param: data - either: a 2D dict, first key sample name (row header),
        second key field (column header); a list of dicts; or a list of lists
    :param: fn - Desired filename. Directory will be prepended automatically.
    :param: sort_cols - Sort columns alphabetically
    :param: data_format - Output format. Defaults to config.data_format (usually tsv)
    :return: None
    """

    if config.data_dir is None:
        return

    # Get data format from config
    if data_format is None:
        data_format = config.data_format

    body = None
    # Some metrics can't be coerced to tab-separated output, test and handle exceptions
    if data_format in ["tsv", "csv"]:
        sep = "\t" if data_format == "tsv" else ","
        # Attempt to reshape data to tsv
        # noinspection PyBroadException
        try:
            # Get all headers from the data, except if data is a dictionary (i.e. has >1 dimensions)
            headers = []
            rows = []

            for d in data.values() if isinstance(data, dict) else data:
                if not d or (isinstance(d, list) and isinstance(d[0], dict)):
                    continue
                if isinstance(d, dict):
                    for h in d.keys():
                        if h not in headers:
                            headers.append(h)
            if headers:
                if sort_cols:
                    headers = sorted(headers)
                headers_str = [str(item) for item in headers]
                if isinstance(data, dict):
                    # Add Sample header as a first element
                    headers_str.insert(0, "Sample")
                rows.append(sep.join(headers_str))

            # The rest of the rows
            for key, d in sorted(data.items()) if isinstance(data, dict) else enumerate(data):
                # Make a list starting with the sample name, then each field in order of the header cols
                if headers:
                    line = [str(d.get(h, "")) for h in headers]
                else:
                    line = [
                        str(item)
                        for item in (d.values() if isinstance(d, dict) else (d if isinstance(d, list) else [d]))
                    ]
                if isinstance(data, dict):
                    # Add Sample header as a first element
                    line.insert(0, str(key))
                rows.append(sep.join(line))
            body = "\n".join(rows)

        except Exception as e:
            if config.development:
                raise
            data_format = "yaml"
            logger.debug(f"{fn} could not be saved as tsv/csv, falling back to YAML. {e}")

    # Add relevant file extension to filename, save file.
    fn = f"{fn}.{config.data_format_extensions[data_format]}"
    fpath = os.path.join(data_tmp_dir(), fn)
    assert data_tmp_dir() and os.path.exists(data_tmp_dir())
    with io.open(fpath, "w", encoding="utf-8") as f:
        if data_format == "json":
            jsonstr = dump_json(data, indent=4, ensure_ascii=False)
            print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
        elif data_format == "yaml":
            yaml.dump(replace_defaultdicts(data), f, default_flow_style=False)
        elif body:
            # Default - tab separated output
            print(body.encode("utf-8", "ignore").decode("utf-8"), file=f)
    logger.debug(f"Wrote data file {fn}")


def dump_json(data, **kwargs):
    """
    Recursively replace non-JSON-conforming NaNs and lambdas with None.
    Note that a custom JSONEncoder would have worked for lambdas, but not for NaNs: https://stackoverflow.com/a/28640141
    """

    # Recursively replace NaNs with None
    def replace_nan(obj):
        if isinstance(obj, dict):
            return {k: replace_nan(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [replace_nan(v) for v in obj]
        elif isinstance(obj, set):
            return {replace_nan(v) for v in obj}
        elif callable(obj):
            return None
        elif isinstance(obj, float) and math.isnan(obj):
            return None
        return obj

    return json.dumps(replace_nan(data), **kwargs)


def multiqc_dump_json():
    """
    Export the parsed data in memory to a JSON file.
    Used for MegaQC and other data export.
    WARNING: May be depreciated and removed in future versions.
    """
    exported_data = dict()
    export_vars = {
        "report": [
            "data_sources",
            "general_stats_data",
            "general_stats_headers",
            "multiqc_command",
            "plot_data",
            "saved_raw_data",
        ],
        "config": [
            "analysis_dir",
            "creation_date",
            "git_hash",
            "intro_text",
            "report_comment",
            "report_header_info",
            "script_path",
            "short_version",
            "subtitle",
            "title",
            "version",
            "output_dir",
        ],
    }
    for s in export_vars:
        for k in export_vars[s]:
            try:
                d = None
                if s == "config":
                    d = {f"{s}_{k}": getattr(config, k)}
                elif s == "report":
                    d = {f"{s}_{k}": getattr(sys.modules[__name__], k)}
                if d:
                    dump_json(d, ensure_ascii=False)  # Test that exporting to JSON works
                    exported_data.update(d)
            except (TypeError, KeyError, AttributeError) as e:
                logger.warning(f"Couldn't export data key '{s}.{k}': {e}")
        # Get the absolute paths of analysis directories
        exported_data["config_analysis_dir_abs"] = list()
        for d in exported_data.get("config_analysis_dir", []):
            try:
                exported_data["config_analysis_dir_abs"].append(os.path.abspath(d))
            except Exception:
                pass
    return exported_data
