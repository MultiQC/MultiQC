#!/usr/bin/env python

""" MultiQC report module. Holds the output from each
module. Is available to subsequent modules. Contains
helper functions to generate markup for report. """


import fnmatch
import inspect
import io
import json
import mimetypes
import os
import re
import time
from collections import OrderedDict, defaultdict
from typing import Dict

import lzstring
import rich
import rich.progress
import yaml

from . import config

logger = config.logger

# Treat defaultdict and OrderedDict as normal dicts for YAML output
from yaml.representer import Representer

yaml.add_representer(defaultdict, Representer.represent_dict)
yaml.add_representer(OrderedDict, Representer.represent_dict)


class Report:
    def __init__(
        self,
        multiqc_command: str,
    ):
        self.multiqc_command = multiqc_command

        self.general_stats_data = list()
        self.general_stats_headers = list()
        self.general_stats_html = ""
        self.data_sources = defaultdict(lambda: defaultdict(lambda: defaultdict()))
        self.plot_data = dict()
        self.html_ids = list()
        self.lint_errors = list()
        self.num_hc_plots = 0
        self.num_mpl_plots = 0
        self.saved_raw_data = dict()
        self.last_found_file = None
        self.runtimes = {
            "total": 0,
            "total_sp": 0,
            "total_mods": 0,
            "total_compression": 0,
            "sp": defaultdict(),
            "mods": defaultdict(),
        }
        self.file_search_stats = {
            "skipped_symlinks": 0,
            "skipped_not_a_file": 0,
            "skipped_ignore_pattern": 0,
            "skipped_filesize_limit": 0,
            "skipped_module_specific_max_filesize": 0,
            "skipped_no_match": 0,
            "skipped_directory_fn_ignore_dirs": 0,
            "skipped_file_contents_search_errors": 0,
        }
        self.searchfiles = list()

        # Make a dict of discovered files for each search key
        self.files = dict()
        self._modules = []

    def add_module(self, mod):
        """
        Build the sections implemented in the module.
        """
        mod.report = self
        self._modules.append(mod)

        # Sanitise anchor ID and check for duplicates
        mod.anchor = self.save_htmlid(mod.anchor)

        # See if we have a user comment in the config
        if mod.anchor in config.section_comments:
            mod.comment = config.section_comments[mod.anchor]

        mod.build()

    def get_filelist(self, run_module_names):
        """
        Go through all supplied search directories and assembly a master
        list of files to search. Then fire search functions for each file.
        """
        # Prep search patterns
        spatterns = [{}, {}, {}, {}, {}, {}, {}]
        self.runtimes["sp"] = defaultdict()
        ignored_patterns = []
        skipped_patterns = []
        for key, sps in config.sp.items():
            mod_name = key.split("/", 1)[0]
            if mod_name.lower() not in [m.lower() for m in run_module_names]:
                ignored_patterns.append(key)
                continue
            self.files[key] = list()
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
                logger.warning(
                    "Unrecognised search pattern keys for '{}': {}".format(key, ", ".join(unrecognised_keys))
                )

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

        def add_file(fn, root):
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
                        if self.search_file(sp, f, key):
                            # Check that we shouldn't exclude this file
                            if not self.exclude_file(sp, f):
                                # Looks good! Remember this file
                                self.files[key].append(f)
                                self.file_search_stats[key] = self.file_search_stats.get(key, 0) + 1
                                file_matched = True
                            # Don't keep searching this file for other modules
                            if not sp.get("shared", False):
                                self.runtimes["sp"][key] = self.runtimes["sp"].get(key, 0) + (time.time() - start)
                                return True
                            # Don't look at other patterns for this module
                            else:
                                break
                    self.runtimes["sp"][key] = self.runtimes["sp"].get(key, 0) + (time.time() - start)

            return file_matched

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
                self.file_search_stats["skipped_symlinks"] += 1
                continue
            elif os.path.isfile(path):
                self.searchfiles.append([os.path.basename(path), os.path.dirname(path)])
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
                            self.file_search_stats["skipped_directory_fn_ignore_dirs"] += len(removed_dirs)
                            orig_dirnames = dirnames[:]
                    for n in config.fn_ignore_paths:
                        dirnames[:] = [
                            d for d in dirnames if not fnmatch.fnmatch(os.path.join(root, d), n.rstrip(os.sep))
                        ]
                        if len(orig_dirnames) != len(dirnames):
                            removed_dirs = [
                                os.path.join(root, d) for d in set(orig_dirnames).symmetric_difference(set(dirnames))
                            ]
                            self.file_search_stats["skipped_directory_fn_ignore_dirs"] += len(removed_dirs)

                    # Skip *this* directory if matches ignore params
                    d_matches = [n for n in config.fn_ignore_dirs if fnmatch.fnmatch(bname, n.rstrip(os.sep))]
                    if len(d_matches) > 0:
                        self.file_search_stats["skipped_directory_fn_ignore_dirs"] += 1
                        continue
                    p_matches = [n for n in config.fn_ignore_paths if fnmatch.fnmatch(root, n.rstrip(os.sep))]
                    if len(p_matches) > 0:
                        self.file_search_stats["skipped_directory_fn_ignore_dirs"] += 1
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
                        self.searchfiles.append([fn, root])

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
            mqc_task = progress.add_task("searching", total=len(self.searchfiles), s_fn="")
            for sf in self.searchfiles:
                progress.update(mqc_task, advance=1, s_fn=os.path.join(sf[1], sf[0])[-50:])
                if not add_file(sf[0], sf[1]):
                    self.file_search_stats["skipped_no_match"] += 1
            progress.update(mqc_task, s_fn="")

        self.runtimes["total_sp"] = time.time() - total_sp_starttime
        if config.profile_runtime:
            logger.info(f"Profile-runtime: Searching files took {self.runtimes['total_sp']:.2f}s")

        # Debug log summary about what we skipped
        summaries = []
        for key in sorted(self.file_search_stats, key=self.file_search_stats.get, reverse=True):
            if "skipped_" in key and self.file_search_stats[key] > 0:
                summaries.append(f"{key}: {self.file_search_stats[key]}")
        logger.debug(f"Summary of files that were skipped by the search: [{'] // ['.join(summaries)}]")

    def search_file(self, pattern, f, module_key):
        """
        Function to search a single file for a single search pattern.
        """

        fn_matched = False
        contents_matched = False

        # Search pattern specific filesize limit
        if pattern.get("max_filesize") is not None and "filesize" in f:
            if f["filesize"] > pattern.get("max_filesize"):
                self.file_search_stats["skipped_module_specific_max_filesize"] += 1
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
        if pattern.get("contents") is not None or pattern.get("contents_re") is not None:
            if pattern.get("contents_re") is not None:
                repattern = re.compile(pattern["contents_re"])
            if "contents_lines" not in f or (
                "num_lines" in pattern and len(f["contents_lines"]) < pattern["num_lines"]
            ):
                f["contents_lines"] = []
                file_path = os.path.join(f["root"], f["fn"])
                try:
                    with io.open(file_path, "r", encoding="utf-8") as fh:
                        for i, line in enumerate(fh):
                            f["contents_lines"].append(line)
                            if i >= config.filesearch_lines_limit and i >= pattern.get("num_lines", 0):
                                break
                # Can't open file - usually because it's a binary file, and we're reading as utf-8
                except (IOError, OSError, ValueError, UnicodeDecodeError) as e:
                    if config.report_readerrors:
                        logger.debug(f"Couldn't read file when looking for output: {file_path}, {e}")
                    self.file_search_stats["skipped_file_contents_search_errors"] += 1
                    return False

            # Go through the parsed file contents
            for i, line in enumerate(f["contents_lines"]):
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
                # Break if we've searched enough lines for this pattern
                if pattern.get("num_lines") and i >= pattern.get("num_lines"):
                    break

        return fn_matched and contents_matched

    def exclude_file(self, sp, f):
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

    def data_sources_tofile(self):
        fn = "multiqc_sources.{}".format(config.data_format_extensions[config.data_format])
        with io.open(os.path.join(config.data_dir, fn), "w", encoding="utf-8") as f:
            if config.data_format == "json":
                jsonstr = json.dumps(self.data_sources, indent=4, ensure_ascii=False)
                print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
            elif config.data_format == "yaml":
                yaml.dump(self.data_sources, f, default_flow_style=False)
            else:
                lines = [["Module", "Section", "Sample Name", "Source"]]
                for mod in self.data_sources:
                    for sec in self.data_sources[mod]:
                        for s_name, source in self.data_sources[mod][sec].items():
                            lines.append([mod, sec, s_name, source])
                body = "\n".join(["\t".join(l) for l in lines])
                print(body.encode("utf-8", "ignore").decode("utf-8"), file=f)

    def dois_tofile(self):
        """Find all DOIs listed in report sections and write to a file"""
        # Collect DOIs
        dois = {"MultiQC": ["10.1093/bioinformatics/btw354"]}
        for mod in self.get_modules():
            if mod.doi is not None and mod.doi != "" and mod.doi != []:
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

    def save_htmlid(self, html_id, skiplint=False):
        """Take a HTML ID, sanitise for HTML, check for duplicates and save.
        Returns sanitised, unique ID"""
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
            if html_id != html_id_clean:
                errmsg = "LINT: {}HTML ID was not clean ('{}' -> '{}') ## {}".format(
                    modname, html_id, html_id_clean, codeline
                )
                logger.error(errmsg)
                self.lint_errors.append(errmsg)

        # Check for duplicates
        i = 1
        html_id_base = html_id_clean
        while html_id_clean in self.html_ids:
            html_id_clean = "{}-{}".format(html_id_base, i)
            i += 1
            if config.lint and not skiplint:
                errmsg = "LINT: {}HTML ID was a duplicate ({}) ## {}".format(modname, html_id_clean, codeline)
                logger.error(errmsg)
                self.lint_errors.append(errmsg)

        # Remember and return
        self.html_ids.append(html_id_clean)
        return html_id_clean

    def sort_modules(self, report_section_order: Dict[str, Dict]):
        """
        Sort the report module output according to the config.
        """
        section_id_order = {}
        idx = 10
        for mod in reversed(self._modules):
            section_id_order[mod.anchor] = idx
            idx += 10
        for anchor, ss in report_section_order.items():
            if anchor not in section_id_order.keys():
                logger.debug("Reordering sections: anchor '{}' not found.".format(anchor))
                continue
            if ss.get("order") is not None:
                section_id_order[anchor] = ss["order"]
            if ss.get("after") in section_id_order.keys():
                section_id_order[anchor] = section_id_order[ss["after"]] + 1
            if ss.get("before") in section_id_order.keys():
                section_id_order[anchor] = section_id_order[ss["before"]] - 1
        sorted_ids = sorted(section_id_order, key=section_id_order.get)
        self._modules = [mod for i in reversed(sorted_ids) for mod in self._modules if mod.anchor == i]

        # Now, sort the report sections within a module
        for midx, mod in enumerate(self._modules):
            section_id_order = {}
            # Get a list of the section anchors
            idx = 10
            for s in mod.sections:
                section_id_order[s["anchor"]] = idx
                idx += 10
            # Go through each section to be reordered
            for anchor, ss in report_section_order.items():
                # Section to be moved is not in this module
                if anchor not in section_id_order.keys():
                    logger.debug("Reordering sections: anchor '{}' not found for module '{}'.".format(anchor, mod.name))
                    continue
                if ss == "remove":
                    section_id_order[anchor] = False
                    continue
                if ss.get("order") is not None:
                    section_id_order[anchor] = ss["order"]
                if ss.get("after") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["after"]] + 1
                if ss.get("before") in section_id_order.keys():
                    section_id_order[anchor] = section_id_order[ss["before"]] - 1
            # Remove module sections
            section_id_order = {s: o for s, o in section_id_order.items() if o is not False}
            # Sort the module sections
            sorted_ids = sorted(section_id_order, key=section_id_order.get)
            self._modules[midx].sections = [s for i in sorted_ids for s in mod.sections if s["anchor"] == i]

    def get_modules(self):
        return self._modules


def compress_json(data):
    """Take a Python data object. Convert to JSON and compress using lzstring"""
    json_string = json.dumps(data).encode("utf-8", "ignore").decode("utf-8")
    json_string = _sanitise_json(json_string)
    x = lzstring.LZString()
    return x.compressToBase64(json_string)


def _sanitise_json(json_string):
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
