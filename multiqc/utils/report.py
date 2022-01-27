#!/usr/bin/env python

""" MultiQC report module. Holds the output from each
module. Is available to subsequent modules. Contains
helper functions to generate markup for report. """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import fnmatch
import inspect
import io
import json
import sys
import threading
import traceback
import lzstring
import mimetypes
import os
import re
import rich.progress
import time
import yaml
import multiprocessing

from multiqc import config

logger = config.logger

# Treat defaultdict and OrderedDict as normal dicts for YAML output
from yaml.representer import Representer, SafeRepresenter

yaml.add_representer(defaultdict, Representer.represent_dict)
yaml.add_representer(OrderedDict, Representer.represent_dict)
try:
    yaml.add_representer(unicode, SafeRepresenter.represent_unicode)
except NameError:
    pass  # Python 3


# Set up global variables shared across modules
# Inside a function so that the global vars are reset if MultiQC is run more than once within a single session / environment
def init():
    global general_stats_data
    general_stats_data = list()

    global general_stats_headers
    general_stats_headers = list()

    global general_stats_html
    general_stats_html = ""

    global data_sources
    data_sources = defaultdict(lambda: defaultdict(lambda: defaultdict()))

    global plot_data
    plot_data = dict()

    global html_ids
    html_ids = list()

    global lint_errors
    lint_errors = list()

    global num_hc_plots
    num_hc_plots = 0

    global num_mpl_plots
    num_mpl_plots = 0

    global saved_raw_data
    saved_raw_data = dict()

    global last_found_file
    last_found_file = None

    global runtimes
    runtimes = {
        "total": 0,
        "total_sp": 0,
        "total_mods": 0,
        "total_compression": 0,
        "sp": defaultdict(),
        "mods": defaultdict(),
    }

    global file_search_stats
    file_search_stats = {
        "skipped_symlinks": 0,
        "skipped_not_a_file": 0,
        "skipped_ignore_pattern": 0,
        "skipped_filesize_limit": 0,
        "skipped_no_match": 0,
    }

    global searchfiles
    searchfiles = list()

    # Make a dict of discovered files for each seach key
    global files
    files = dict()


def get_filelist(run_module_names):
    global file_search_stats, runtimes, files
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
            logger.warning("Unrecognised search pattern keys for '{}': {}".format(key, ", ".join(unrecognised_keys)))

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
        logger.debug("Ignored {} search patterns as didn't match running modules.".format(len(ignored_patterns)))

    if len(skipped_patterns) > 0:
        logger.info("Skipping {} file search patterns".format(len(skipped_patterns)))
        logger.debug("Skipping search patterns: {}".format(", ".join(skipped_patterns)))

    # Go through the analysis directories and get file list
    multiqc_installation_dir_files = [
        "LICENSE",
        "CHANGELOG.md",
        "Dockerfile",
        "MANIFEST.in",
        ".gitmodules",
        "README.md",
        "CSP.txt",
        "setup.py",
        ".gitignore",
    ]

    total_sp_starttime = time.time()
    for path in config.analysis_dir:
        if os.path.islink(path) and config.ignore_symlinks:
            file_search_stats["skipped_symlinks"] += 1
            continue
        elif os.path.isfile(path):
            searchfiles.append([os.path.basename(path), os.path.dirname(path)])
        elif os.path.isdir(path):
            for root, dirnames, filenames in os.walk(path, followlinks=(not config.ignore_symlinks), topdown=True):
                bname = os.path.basename(root)

                # Skip any sub-directories matching ignore params
                orig_dirnames = dirnames[:]
                for n in config.fn_ignore_dirs:
                    dirnames[:] = [d for d in dirnames if not fnmatch.fnmatch(d, n.rstrip(os.sep))]
                    if len(orig_dirnames) != len(dirnames):
                        removed_dirs = [
                            os.path.join(root, d) for d in set(orig_dirnames).symmetric_difference(set(dirnames))
                        ]
                        logger.debug("Ignoring directory as matched fn_ignore_dirs: {}".format(", ".join(removed_dirs)))
                        orig_dirnames = dirnames[:]
                for n in config.fn_ignore_paths:
                    dirnames[:] = [d for d in dirnames if not fnmatch.fnmatch(os.path.join(root, d), n.rstrip(os.sep))]
                    if len(orig_dirnames) != len(dirnames):
                        removed_dirs = [
                            os.path.join(root, d) for d in set(orig_dirnames).symmetric_difference(set(dirnames))
                        ]
                        logger.debug(
                            "Ignoring directory as matched fn_ignore_paths: {}".format(", ".join(removed_dirs))
                        )

                # Skip *this* directory if matches ignore params
                d_matches = [n for n in config.fn_ignore_dirs if fnmatch.fnmatch(bname, n.rstrip(os.sep))]
                if len(d_matches) > 0:
                    logger.debug("Ignoring directory as matched fn_ignore_dirs: {}".format(bname))
                    continue
                p_matches = [n for n in config.fn_ignore_paths if fnmatch.fnmatch(root, n.rstrip(os.sep))]
                if len(p_matches) > 0:
                    logger.debug("Ignoring directory as matched fn_ignore_paths: {}".format(root))
                    continue

                # Sanity check - make sure that we're not just running in the installation directory
                if len(filenames) > 0 and all([fn in filenames for fn in multiqc_installation_dir_files]):
                    logger.error("Error: MultiQC is running in source code directory! {}".format(root))
                    logger.warning(
                        "Please see the docs for how to use MultiQC: https://multiqc.info/docs/#running-multiqc"
                    )
                    dirnames[:] = []
                    filenames[:] = []
                    continue

                # Search filenames in this directory
                for fn in filenames:
                    searchfiles.append([fn, root])

    file_name_search_end = time.time()
    runtimes["file_name_search"] = file_name_search_end - total_sp_starttime
    file_content_search_start = time.time()

    # Search through collected files
    progress_obj = rich.progress.Progress(
        "[blue]|[/]      ",
        rich.progress.SpinnerColumn(),
        "[blue]{task.description}[/] |",
        rich.progress.BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        "[green]{task.completed}/{task.total}",
        "[dim]{task.fields[s_fn]}",
    )

    class SerialFileAdder:
        """
        Container class for running the file addition serially
        """

        def __init__(self, files, file_search_stats, spatterns):
            self.files = files
            self.sp_runtimes = dict()
            self.file_search_stats = file_search_stats
            self.spatterns = spatterns

        # Bit of a hack...
        from .report import add_file

    class FileAdderProcess(multiprocessing.Process):
        def __init__(self, id, spatterns, results_pipe_s, path_pipe_r, _files, queue, prog_queue):
            """
            A worker process for running the file addition in parallel
            """
            super().__init__()
            self.id = id
            self.spatterns = spatterns
            self.results_pipe_s = results_pipe_s
            self.path_pipe_r = path_pipe_r
            self.need_work_queue = queue
            self.prog_queue = prog_queue
            self.files = _files
            self.sp_runtimes = dict()
            self.file_search_stats = {
                "skipped_no_match": 0,
                "skipped_not_a_file": 0,
                "skipped_ignore_pattern": 0,
                "skipped_filesize_limit": 0,
            }

        # Bit of a hack...
        from .report import add_file

        def run(self):
            try:
                chunk = []
                while True:
                    # Ask the main process for more work
                    self.need_work_queue.put(self.id)
                    signal = self.path_pipe_r.recv()
                    last_chunk_len = len(chunk)
                    if signal == -1:
                        # We received an exit signal from the main process
                        # and should notify the progress tracker
                        self.prog_queue.put((2, self.id, last_chunk_len))
                        break
                    chunk = self.path_pipe_r.recv()
                    self.prog_queue.put((1, self.id, last_chunk_len, len(chunk)))

                    # Iterate over the paths in the current chunk and update the results accordingly
                    for fn, root in chunk:
                        if not self.add_file(fn, root):
                            self.file_search_stats["skipped_no_match"] += 1
                        self.prog_queue.put((0, self.id, os.path.join(fn, root)))

                # Everything went well, send the data
                self.results_pipe_s.send(0)
                self.results_pipe_s.send((self.files, self.sp_runtimes, self.file_search_stats))

            except Exception as e:
                # Report the exception to the main process
                self.results_pipe_s.send(-1)
                self.results_pipe_s.send((e, traceback.format_exc()))

    class ProgressTracker(threading.Thread):
        """
        Thread for tracking the progress of the worker processes
        """

        def __init__(self, progress, prog_queue, n_procs, len_files):
            super().__init__()

            self.queue = prog_queue
            self.active = n_procs

            # Initialise the progress bars
            self.progress = progress
            self.main_task = self.progress.add_task("Total    ", total=len_files, s_fn="")
            self.proc_tasks = [self.progress.add_task(f"Process {i}", total=0, s_fn="") for i in range(1, n_procs + 1)]

        def run(self):
            while self.active > 0:
                v = self.queue.get()
                status = v[0]
                id = v[1]
                if status == 0:
                    # Update the worker process task
                    file_name = v[2]
                    self.progress.update(self.proc_tasks[id], advance=1, s_fn=file_name)
                elif status == 1:
                    # A worker process is finished with a chunk
                    last_chunk_len = v[2]
                    curr_chunk_len = v[3]
                    self.progress.update(self.main_task, advance=last_chunk_len)
                    self.progress.update(self.proc_tasks[id], total=curr_chunk_len, completed=0)
                elif status == 2:
                    # The process is finished, just update the main process task
                    last_chunk_len = v[2]
                    self.progress.update(self.main_task, advance=last_chunk_len)
                    self.active -= 1
                elif status == -1:
                    # We got an exception in one of the processes and should exit gracefully
                    return

    with progress_obj as progress:
        if config.cpus == 1:
            mqc_task = progress.add_task("searching", total=len(searchfiles), s_fn="")
            searcher = SerialFileAdder(files, file_search_stats, spatterns)
            for sf in searchfiles:
                progress.update(mqc_task, advance=1, s_fn=os.path.join(sf[1], sf[0])[-50:])
                if not searcher.add_file(sf[0], sf[1]):
                    file_search_stats["skipped_no_match"] += 1
            runtimes["sp"] = searcher.sp_runtimes
            files = searcher.files
            file_search_stats = searcher.file_search_stats
            progress.update(mqc_task, s_fn="")
        else:
            # Set up variables for parallel file search
            results_pipes = []
            path_pipes = []
            worker_procs = []

            # Set up queues for communication between the processes
            need_work_queue = multiprocessing.Queue()
            prog_queue = multiprocessing.Queue()

            # Initialise the worker processes
            for i in range(config.cpus):
                results_pipe_r, results_pipe_s = multiprocessing.Pipe(duplex=False)
                path_pipe_r, path_pipe_s = multiprocessing.Pipe(duplex=False)
                _files = {mod: [] for mod in config.sp.keys()}
                worker_proc = FileAdderProcess(
                    i, spatterns, results_pipe_s, path_pipe_r, _files, need_work_queue, prog_queue
                )
                worker_proc.daemon = True
                worker_proc.start()

                results_pipes.append(results_pipe_r)
                path_pipes.append(path_pipe_s)
                worker_procs.append(worker_proc)

            class Chunker:
                """
                Iterator returning chunks from the file path list based on file size
                """

                def __init__(self, paths, size_threshold):
                    self.paths = paths
                    self.threshold = size_threshold
                    self.i = 0

                def __iter__(self):
                    self.i = 0
                    return self

                def __next__(self):
                    if self.i >= len(self.paths):
                        raise StopIteration

                    curr_size = 0
                    chunk = []
                    while curr_size < self.threshold and self.i < len(self.paths):
                        path = self.paths[self.i]
                        chunk.append(path)
                        curr_size += os.path.getsize(os.path.join(path[1], path[0]))
                        self.i += 1
                    return chunk

            # Set up progress tracker thread
            progress_tracker = ProgressTracker(progress, prog_queue, config.cpus, len(searchfiles))
            progress_tracker.daemon = True
            progress_tracker.start()

            # We iterate over the file paths chunks
            # and send them to the workers
            for chunk in Chunker(searchfiles, 2 ** 21):
                id = need_work_queue.get()
                path_pipes[id].send(0)
                path_pipes[id].send(chunk)

            # Tell the worker processes to exit
            for _ in worker_procs:
                id = need_work_queue.get()
                path_pipes[id].send(-1)

            # Gather the data from the processes
            for pipe in results_pipes:
                exit_status = pipe.recv()
                if exit_status == 0:
                    # Gather the data from the process
                    _files, _sp_runtimes, _file_search_stats = pipe.recv()

                    for key, file_list in _files.items():
                        if key not in files:
                            files[key] = []
                        files[key].extend(file_list)

                    for cat, num in _file_search_stats.items():
                        if cat not in file_search_stats:
                            file_search_stats[cat] = num
                        else:
                            file_search_stats[cat] += num

                    sp_runtimes = runtimes["sp"]
                    for key, sp_runtime in _sp_runtimes.items():
                        sp_runtimes[key] = sp_runtime + sp_runtimes.get(key, 0)
                elif exit_status == -1:
                    # An exception has occurred in the worker

                    # Tell the progress tracker thread to exit
                    prog_queue.send((-1, -1))

                    # Fetch the exception information
                    e, tb_str = pipe.recv()

                    # We should use the worker's actual traceback object
                    # here, but traceback objects are not picklable.
                    logger.error("%s", tb_str)
                    raise e

            # Join the processes
            for worker_proc in worker_procs:
                worker_proc.join()

            # Finally, join the progress tracker thread
            progress_tracker.join()

    file_content_search_end = time.time()
    runtimes["file_content_search"] = file_content_search_end - file_content_search_start
    runtimes["total_sp"] = file_content_search_end - total_sp_starttime


def add_file(self, fn, root):
    """
    Function applied to each file found when walking the analysis
    directories. Runs through all search patterns and returns True
    if a match is found.
    """
    f = {"fn": fn, "root": root}

    # Check that this is a file and not a pipe or anything weird
    if not os.path.isfile(os.path.join(root, fn)):
        self.file_search_stats["skipped_not_a_file"] += 1
        return False

    # Check that we don't want to ignore this file
    i_matches = [n for n in config.fn_ignore_files if fnmatch.fnmatch(fn, n)]
    if len(i_matches) > 0:
        logger.debug("Ignoring file as matched an ignore pattern: {}".format(fn))
        self.file_search_stats["skipped_ignore_pattern"] += 1
        return False

    # Limit search to small files, to avoid 30GB FastQ files etc.
    try:
        f["filesize"] = os.path.getsize(os.path.join(root, fn))
    except (IOError, OSError, ValueError, UnicodeDecodeError):
        logger.debug("Couldn't read file when checking filesize: {}".format(fn))
    else:
        if f["filesize"] > config.log_filesize_limit:
            self.file_search_stats["skipped_filesize_limit"] += 1
            return False

    # Test file for each search pattern
    file_matched = False
    for patterns in self.spatterns:
        for key, sps in patterns.items():
            start = time.time()
            for sp in sps:
                if search_file(sp, f, key):
                    # Check that we shouldn't exclude this file
                    if not exclude_file(sp, f):
                        # Looks good! Remember this file
                        self.files[key].append(f)
                        self.file_search_stats[key] = file_search_stats.get(key, 0) + 1
                        file_matched = True
                    # Don't keep searching this file for other modules
                    if not sp.get("shared", False):
                        self.sp_runtimes[key] = self.sp_runtimes.get(key, 0) + (time.time() - start)
                        return True
                    # Don't look at other patterns for this module
                    else:
                        break
            self.sp_runtimes[key] = self.sp_runtimes.get(key, 0) + (time.time() - start)

    return file_matched


def search_file(pattern, f, module_key):
    """
    Function to search a single file for a single search pattern.
    """

    fn_matched = False
    contents_matched = False

    # Use mimetypes to exclude binary files where possible
    if not re.match(r".+_mqc\.(png|jpg|jpeg)", f["fn"]) and config.ignore_images:
        (ftype, encoding) = mimetypes.guess_type(os.path.join(f["root"], f["fn"]))
        if encoding is not None:
            return False
        if ftype is not None and ftype.startswith("image"):
            return False

    # Search pattern specific filesize limit
    if pattern.get("max_filesize") is not None and "filesize" in f:
        if f["filesize"] > pattern.get("max_filesize"):
            logger.debug(
                "File ignored by {} because it exceeded search pattern filesize limit: {}".format(module_key, f["fn"])
            )
            return False

    # Search by file name (glob)
    if pattern.get("fn") is not None:
        if fnmatch.fnmatch(f["fn"], pattern["fn"]):
            fn_matched = True
            if pattern.get("contents") is None and pattern.get("contents_re") is None:
                return True

    # Search by file name (regex)
    if pattern.get("fn_re") is not None:
        if re.match(pattern["fn_re"], f["fn"]):
            fn_matched = True
            if pattern.get("contents") is None and pattern.get("contents_re") is None:
                return True

    # Search by file contents
    if pattern.get("contents") is not None or pattern.get("contents_re") is not None:
        if pattern.get("contents_re") is not None:
            repattern = re.compile(pattern["contents_re"])

        try:
            with io.open(os.path.join(f["root"], f["fn"]), "r", encoding="utf-8") as f:
                l = 1
                for line in f:
                    # Search by file contents (string)
                    if pattern.get("contents") is not None:
                        if pattern["contents"] in line:
                            contents_matched = True
                            if pattern.get("fn") is None and pattern.get("fn_re") is None:
                                return True
                            break
                    # Search by file contents (regex)
                    elif "contents_re" in pattern:
                        if re.search(repattern, line):
                            contents_matched = True
                            if pattern.get("fn") is None and pattern.get("fn_re") is None:
                                return True
                            break
                    # Break if we've searched enough lines for this pattern
                    if pattern.get("num_lines") and l >= pattern.get("num_lines"):
                        break
                    l += 1
        except (IOError, OSError, ValueError, UnicodeDecodeError):
            if config.report_readerrors:
                logger.debug("Couldn't read file when looking for output: {}".format(f["fn"]))
                return False

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
    fn = "multiqc_sources.{}".format(config.data_format_extensions[config.data_format])
    with io.open(os.path.join(config.data_dir, fn), "w", encoding="utf-8") as f:
        if config.data_format == "json":
            jsonstr = json.dumps(data_sources, indent=4, ensure_ascii=False)
            print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
        elif config.data_format == "yaml":
            yaml.dump(data_sources, f, default_flow_style=False)
        else:
            lines = [["Module", "Section", "Sample Name", "Source"]]
            for mod in data_sources:
                for sec in data_sources[mod]:
                    for s_name, source in data_sources[mod][sec].items():
                        lines.append([mod, sec, s_name, source])
            body = "\n".join(["\t".join(l) for l in lines])
            print(body.encode("utf-8", "ignore").decode("utf-8"), file=f)


def dois_tofile():
    """Find all DOIs listed in report sections and write to a file"""
    # Collect DOIs
    dois = {"MultiQC": ["10.1093/bioinformatics/btw354"]}
    for mod in modules_output:
        if mod.doi is not None and mod.doi != []:
            dois[mod.anchor] = mod.doi
    # Write to a file
    fn = "multiqc_citations.{}".format(config.data_format_extensions[config.data_format])
    with io.open(os.path.join(config.data_dir, fn), "w", encoding="utf-8") as f:
        if config.data_format == "json":
            jsonstr = json.dumps(dois, indent=4, ensure_ascii=False)
            print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
        elif config.data_format == "yaml":
            yaml.dump(dois, f, default_flow_style=False)
        else:
            body = ""
            for mod, dois in dois.items():
                for doi in dois:
                    body += "{}{} # {}\n".format(doi, " " * (50 - len(doi)), mod)
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
        html_id_clean = "mqc_{}".format(html_id_clean)

    # Replace illegal characters
    html_id_clean = re.sub("[^a-zA-Z0-9_-]+", "_", html_id_clean)

    # Validate if linting
    if config.lint and not skiplint:
        modname = ""
        codeline = ""
        callstack = inspect.stack()
        for n in callstack:
            if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
                callpath = n[1].split("multiqc/modules/", 1)[-1]
                modname = ">{}< ".format(callpath)
                codeline = n[4][0].strip()
                break
    if config.lint and not skiplint and html_id != html_id_clean:
        errmsg = "LINT: {}HTML ID was not clean ('{}' -> '{}') ## {}".format(modname, html_id, html_id_clean, codeline)
        logger.error(errmsg)
        lint_errors.append(errmsg)

    # Check for duplicates
    i = 1
    html_id_base = html_id_clean
    while html_id_clean in html_ids:
        html_id_clean = "{}-{}".format(html_id_base, i)
        i += 1
        if config.lint and not skiplint:
            errmsg = "LINT: {}HTML ID was a duplicate ({}) ## {}".format(modname, html_id_clean, codeline)
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
