"""MultiQC report module. Holds the output from each
module. Is available to subsequent modules. Contains
helper functions to generate markup for report."""

import base64
import fnmatch
import gzip
import inspect
import io
import json
import mimetypes
import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union, List, Optional

import rich
import rich.progress
import yaml


from . import config
from .util_functions import replace_defaultdicts, dump_json

logger = config.logger


# Uninitialised global variables for static typing
multiqc_command: str
modules_output: List  # List of BaseMultiqcModule objects
general_stats_data: List[Dict]
general_stats_headers: List[Dict]
general_stats_html: str
data_sources: Dict[str, Dict[str, Dict]]
plot_data: Dict
html_ids: List[str]
lint_errors: List[str]
num_hc_plots: int
num_mpl_plots: int
saved_raw_data: Dict
last_found_file: Optional[str]
runtimes: Dict[str, Union[float, Dict]]
file_search_stats: Dict[str, int]
searchfiles: List
files: Dict
software_versions: Dict[str, Dict[str, List]]


def init():
    # Set up global variables shared across modules. Inside a function so that the global
    # vars are reset if MultiQC is run more than once within a single session / environment.
    global multiqc_command
    global modules_output
    global general_stats_data
    global general_stats_headers
    global general_stats_html
    global data_sources
    global plot_data
    global html_ids
    global lint_errors
    global num_hc_plots
    global num_mpl_plots
    global saved_raw_data
    global last_found_file
    global runtimes
    global file_search_stats
    global searchfiles
    global files
    global software_versions

    multiqc_command = ""
    modules_output = list()
    general_stats_data = list()
    general_stats_headers = list()
    general_stats_html = ""
    data_sources = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    plot_data = dict()
    html_ids = list()
    lint_errors = list()
    num_hc_plots = 0
    num_mpl_plots = 0
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
    searchfiles = list()
    # Discovered files for each search key
    files = dict()
    # Map of Software tools to a set of unique version strings
    software_versions = defaultdict(lambda: defaultdict(list))


class SearchFile:
    """
    Wrap file handler and provide a lazy line iterator on it with caching.

    The class is a context manager with context being an open file handler.

    The `self.line_iterator` method is a distinct context manager over the file lines,
    and thanks for caching can be called multiple times for the same open file handler.
    """

    filename: str
    root: str
    path: str
    _filesize: Optional[int]

    def __init__(self, filename: str, root: str):
        self.filename = filename
        self.root = root
        self.path = os.path.join(root, filename)
        self._filesize = None
        self.content_lines = []
        self._filehandle = None

    @property
    def filesize(self) -> Optional[int]:
        if self._filesize is None:
            try:
                self._filesize = os.path.getsize(self.path)
            except (IOError, OSError, ValueError, UnicodeDecodeError):
                logger.debug(f"Couldn't read file when checking filesize: {self.filename}")
                self._filesize = None
        return self._filesize

    def line_iterator(self):
        """
        Iterate over `self.content_lines`, try to read more lines from file handler when needed.

        Essentially it's a lazy `f.readlines()` with caching.

        Can be called multiple times for the same open file handler.
        """
        for line in self.content_lines:
            yield line
        if self._filehandle is None:
            try:
                self._filehandle = io.open(self.path, "rt", encoding="utf-8")
            except Exception as e:
                if config.report_readerrors:
                    logger.debug(f"Couldn't read file when looking for output: {self.path}, {e}")
        try:
            for line in self._filehandle:
                self.content_lines.append(line)
                yield line

        except UnicodeDecodeError as e:
            if config.report_readerrors:
                logger.debug(
                    f"Couldn't read file as a utf-8 text when looking for output: {self.path}, {e}. "
                    f"Usually because it's a binary file. But sometimes there are single non-unicode "
                    f"characters, so attempting reading while skipping such characters."
                )
        except Exception as e:
            if config.report_readerrors:
                logger.debug(f"Couldn't read file when looking for output: {self.path}, {e}")
            raise
        else:
            # When no lines are parsed, self.content_lines should be empty
            if not self.content_lines and config.report_readerrors:
                logger.debug(f"No utf-8 lines were read from the file, skipping {self.path}")
            return  # No errors.
        self._filehandle.close()
        self._filehandle = io.open(self.path, "rt", encoding="utf-8", errors="ignore")
        try:
            if self.content_lines:
                # Skip all the lines that were already read
                self._filehandle.readlines(len(self.content_lines))
            for line in self._filehandle:
                self.content_lines.append(line)
                yield line
        except Exception as e:
            if config.report_readerrors:
                logger.debug(f"Still couldn't read the file, skipping: {self.path}, {e}")

        if not self.content_lines and config.report_readerrors:
            logger.debug(f"No utf-8 lines were read from the file, skipping {self.path}")

    def close(self):
        if self._filehandle:
            self._filehandle.close()
        self._filehandle = None

    def closed(self):
        return bool(self._filehandle)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def to_dict(self):
        return {"fn": self.filename, "root": self.root}


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

    def _sort_by_key(sps: Dict[str, List[Dict[str, str]]], key: str) -> Dict:
        """Sort search patterns dict by key like num_lines or max_filesize."""

        def _sort_key(kv) -> int:
            _, sps = kv
            # Since modules can have multiple search patterns, we sort by the slowest one:
            return max([sp.get(key, 0) for sp in sps])

        return dict(sorted(sps.items(), key=_sort_key))

    # Sort patterns for faster access. File searches with fewer lines or
    # smaller file sizes go first.
    sorted_spatterns = [{}, {}, {}, {}, {}, {}, {}]
    sorted_spatterns[0] = spatterns[0]  # Only filename matching
    sorted_spatterns[1] = _sort_by_key(spatterns[1], "num_lines")
    sorted_spatterns[2] = _sort_by_key(spatterns[2], "max_filesize")
    sorted_spatterns[3] = spatterns[3]
    sorted_spatterns[4] = _sort_by_key(spatterns[4], "num_lines")
    sorted_spatterns[5] = _sort_by_key(spatterns[5], "max_filesize")
    sorted_spatterns[6] = spatterns[6]
    spatterns = sorted_spatterns

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
        f = SearchFile(fn, root)

        # Check that this is a file and not a pipe or anything weird
        if not os.path.isfile(os.path.join(root, fn)):
            file_search_stats["skipped_not_a_file"] += 1
            return False

        # Check that we don't want to ignore this file
        i_matches = [n for n in config.fn_ignore_files if fnmatch.fnmatch(fn, n)]
        if len(i_matches) > 0:
            file_search_stats["skipped_ignore_pattern"] += 1
            return False

        if f.filesize is not None and f.filesize > config.log_filesize_limit:
            file_search_stats["skipped_filesize_limit"] += 1
            return False

        # Use mimetypes to exclude binary files where possible
        if not re.match(r".+_mqc\.(png|jpg|jpeg)", f.filename) and config.ignore_images:
            (ftype, encoding) = mimetypes.guess_type(f.path)
            if encoding is not None:
                return False
            if ftype is not None and ftype.startswith("image"):
                return False

        # Test file for each search pattern
        file_matched = False
        with f:  # Ensure any open filehandles are closed.
            for patterns in spatterns:
                for key, sps in patterns.items():
                    start = time.time()
                    for sp in sps:
                        if search_file(sp, f, key):
                            # Check that we shouldn't exclude this file
                            if not exclude_file(sp, f):
                                # Looks good! Remember this file
                                files[key].append(f.to_dict())
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

    # Search through collected files
    console = rich.console.Console(
        stderr=True,
        highlight=False,
        force_interactive=False if config.no_ansi else None,
        color_system=None if config.no_ansi else "auto",
    )
    progress_obj = rich.progress.Progress(
        "[blue]|[/]      ",
        rich.progress.SpinnerColumn(),
        "[blue]{task.description}[/] |",
        rich.progress.BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        "[green]{task.completed}/{task.total}",
        "[dim]{task.fields[s_fn]}",
        console=console,
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
    logger.debug(f"Summary of files that were skipped by the search: [{'] // ['.join(summaries)}]")


def search_file(pattern, f: SearchFile, module_key):
    """
    Function to searach a single file for a single search pattern.
    """

    global file_search_stats
    fn_matched = False
    contents_matched = False

    # Search pattern specific filesize limit
    max_filesize = pattern.get("max_filesize")
    if max_filesize is not None and f.filesize:
        if f.filesize > max_filesize:
            file_search_stats["skipped_module_specific_max_filesize"] += 1
            return False

    # Search by file name (glob)
    filename_pattern = pattern.get("fn")
    if filename_pattern is not None:
        if fnmatch.fnmatch(f.filename, filename_pattern):
            fn_matched = True
        else:
            return False

    # Search by file name (regex)
    filename_regex_pattern = pattern.get("fn_re")
    if filename_regex_pattern is not None:
        if re.match(filename_regex_pattern, f.filename):
            fn_matched = True
        else:
            return False

    contents_pattern = pattern.get("contents")
    contents_regex_pattern = pattern.get("contents_re")
    if fn_matched and contents_pattern is None and contents_regex_pattern is None:
        return True

    # Search by file contents
    if contents_pattern is not None or contents_regex_pattern is not None:
        if pattern.get("contents_re") is not None:
            repattern = re.compile(pattern["contents_re"])
        else:
            repattern = None
        num_lines = pattern.get("num_lines", config.filesearch_lines_limit)
        expected_contents = pattern.get("contents")
        try:
            # Zipping with range(num_lines) halts at the desired number.
            for i, line in zip(range(num_lines), f.line_iterator()):
                if expected_contents and expected_contents in line:
                    contents_matched = True
                    break
                if repattern and repattern.match(line):
                    contents_matched = True
                    break
        except Exception:
            file_search_stats["skipped_file_contents_search_errors"] += 1
            return False
        if filename_pattern is None and filename_regex_pattern is None and contents_matched:
            return True
        return fn_matched and contents_matched


def exclude_file(sp, f: SearchFile):
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
            if fnmatch.fnmatch(f.filename, pat):
                return True

    # Search by file name (regex)
    if "exclude_fn_re" in sp:
        for pat in sp["exclude_fn_re"]:
            if re.match(pat, f.filename):
                return True

    # Search the contents of the file
    if "exclude_contents" in sp or "exclude_contents_re" in sp:
        # Compile regex patterns if we have any
        if "exclude_contents_re" in sp:
            sp["exclude_contents_re"] = [re.compile(pat) for pat in sp["exclude_contents_re"]]
        for line in f.line_iterator():
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
            json.dump(data_sources, f, indent=4, ensure_ascii=False)
        elif config.data_format == "yaml":
            # Unlike JSON, YAML represents defaultdicts as objects, so need to convert
            # them to normal dicts
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
    with open(os.path.join(config.data_dir, fn), "w") as f:
        if config.data_format == "json":
            json.dump(dois, f, indent=4, ensure_ascii=False)
        elif config.data_format == "yaml":
            # Unlike JSON, YAML represents defaultdicts as objects, so need to convert
            # them to normal dicts
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
    """
    Take a Python data object. Convert to JSON and compress using gzip.
    Represent in base64 format.
    """
    json_string = dump_json(data)
    json_bytes = json_string.encode("utf-8")
    json_gzip = gzip.compress(json_bytes)
    base64_bytes = base64.b64encode(json_gzip)
    return base64_bytes.decode("ascii")
