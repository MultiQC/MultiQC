#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

from __future__ import print_function
from collections import OrderedDict
import io
import fnmatch
import logging
import markdown
import mimetypes
import os
import re
import textwrap

from multiqc.utils import report, config, util_functions

logger = logging.getLogger(__name__)


class BaseMultiqcModule(object):
    def __init__(
        self,
        name="base",
        anchor="base",
        target=None,
        href=None,
        info=None,
        comment=None,
        extra=None,
        autoformat=True,
        autoformat_type="markdown",
    ):

        # Custom options from user config that can overwrite base module values
        mod_cust_config = getattr(self, "mod_cust_config", {})
        self.name = mod_cust_config.get("name", name)
        self.anchor = mod_cust_config.get("anchor", anchor)
        target = mod_cust_config.get("target", target)
        self.href = mod_cust_config.get("href", href)
        self.info = mod_cust_config.get("info", info)
        self.comment = mod_cust_config.get("comment", comment)
        self.extra = mod_cust_config.get("extra", extra)
        # Specific module level config to overwrite (e.g. config.bcftools, config.fastqc)
        config.update({anchor: mod_cust_config.get("custom_config", {})})

        # Sanitise anchor ID and check for duplicates
        self.anchor = report.save_htmlid(self.anchor)

        # See if we have a user comment in the config
        if self.anchor in config.section_comments:
            self.comment = config.section_comments[self.anchor]

        if self.info is None:
            self.info = ""
        if self.extra is None:
            self.extra = ""
        if target is None:
            target = self.name
        if self.href is not None:
            self.mname = '<a href="{}" target="_blank">{}</a>'.format(self.href, target)
        else:
            self.mname = target
        if self.href or self.info or self.extra:
            self.intro = "<p>{} {}</p>{}".format(self.mname, self.info, self.extra)

        # Format the markdown strings
        if autoformat:
            if self.comment is not None:
                self.comment = textwrap.dedent(self.comment)
                if autoformat_type == "markdown":
                    self.comment = markdown.markdown(self.comment)

        self.sections = list()

    def find_log_files(self, sp_key, filecontents=True, filehandles=False):
        """
        Return matches log files of interest.
        :param sp_key: Search pattern key specified in config
        :param filehandles: Set to true to return a file handle instead of slurped file contents
        :return: Yields a dict with filename (fn), root directory (root), cleaned sample name
                 generated from the filename (s_name) and either the file contents or file handle
                 for the current matched file (f).
                 As yield is used, the results can be iterated over without loading all files at once
        """

        # Pick up path filters if specified.
        # Allows modules to be called multiple times with different sets of files
        path_filters = getattr(self, "mod_cust_config", {}).get("path_filters")
        path_filters_exclude = getattr(self, "mod_cust_config", {}).get("path_filters_exclude")

        # Old, depreciated syntax support. Likely to be removed in a future version.
        if isinstance(sp_key, dict):
            report.files[self.name] = list()
            for sf in report.searchfiles:
                if report.search_file(sp_key, {"fn": sf[0], "root": sf[1]}, module_key=None):
                    report.files[self.name].append({"fn": sf[0], "root": sf[1]})
            sp_key = self.name
            logwarn = "Depreciation Warning: {} - Please use new style for find_log_files()".format(self.name)
            if len(report.files[self.name]) > 0:
                logger.warning(logwarn)
            else:
                logger.debug(logwarn)
        elif not isinstance(sp_key, str):
            logger.warning("Did not understand find_log_files() search key")
            return

        for f in report.files[sp_key]:
            # Make a note of the filename so that we can report it if something crashes
            report.last_found_file = os.path.join(f["root"], f["fn"])

            # Filter out files based on exclusion patterns
            if path_filters_exclude and len(path_filters_exclude) > 0:
                exlusion_hits = (fnmatch.fnmatch(report.last_found_file, pfe) for pfe in path_filters_exclude)
                if any(exlusion_hits):
                    logger.debug(
                        "{} - Skipping '{}' as it matched the path_filters_exclude for '{}'".format(
                            sp_key, f["fn"], self.name
                        )
                    )
                    continue

            # Filter out files based on inclusion patterns
            if path_filters and len(path_filters) > 0:
                inclusion_hits = (fnmatch.fnmatch(report.last_found_file, pf) for pf in path_filters)
                if not any(inclusion_hits):
                    logger.debug(
                        "{} - Skipping '{}' as it didn't match the path_filters for '{}'".format(
                            sp_key, f["fn"], self.name
                        )
                    )
                    continue
                else:
                    logger.debug(
                        "{} - Selecting '{}' as it matched the path_filters for '{}'".format(sp_key, f["fn"], self.name)
                    )

            # Make a sample name from the filename
            f["sp_key"] = sp_key
            f["s_name"] = self.clean_s_name(f["fn"], f)
            if filehandles or filecontents:
                try:
                    # Custom content module can now handle image files
                    (ftype, encoding) = mimetypes.guess_type(os.path.join(f["root"], f["fn"]))
                    if ftype is not None and ftype.startswith("image"):
                        with io.open(os.path.join(f["root"], f["fn"]), "rb") as fh:
                            # always return file handles
                            f["f"] = fh
                            yield f
                    else:
                        # Everything else - should be all text files
                        with io.open(os.path.join(f["root"], f["fn"]), "r", encoding="utf-8") as fh:
                            if filehandles:
                                f["f"] = fh
                                yield f
                            elif filecontents:
                                f["f"] = fh.read()
                                yield f
                except (IOError, OSError, ValueError, UnicodeDecodeError) as e:
                    if config.report_readerrors:
                        logger.debug("Couldn't open filehandle when returning file: {}\n{}".format(f["fn"], e))
                        f["f"] = None
            else:
                yield f

    def add_section(
        self,
        name=None,
        anchor=None,
        description="",
        comment="",
        helptext="",
        plot="",
        content="",
        autoformat=True,
        autoformat_type="markdown",
    ):
        """Add a section to the module report output"""

        # Default anchor
        if anchor is None:
            if name is not None:
                nid = name.lower().strip().replace(" ", "-")
                anchor = "{}-{}".format(self.anchor, nid)
            else:
                sl = len(self.sections) + 1
                anchor = "{}-section-{}".format(self.anchor, sl)

        # Append custom module anchor to the section if set
        mod_cust_config = getattr(self, "mod_cust_config", {})
        if "anchor" in mod_cust_config:
            anchor = "{}_{}".format(mod_cust_config["anchor"], anchor)

        # Sanitise anchor ID and check for duplicates
        anchor = report.save_htmlid(anchor)

        # Skip if user has a config to remove this module section
        if anchor in config.remove_sections:
            logger.debug("Skipping section '{}' because specified in user config".format(anchor))
            return

        # See if we have a user comment in the config
        if anchor in config.section_comments:
            comment = config.section_comments[anchor]

        # Format the content
        if autoformat:
            if len(description) > 0:
                description = textwrap.dedent(description)
                if autoformat_type == "markdown":
                    description = markdown.markdown(description)
            if len(comment) > 0:
                comment = textwrap.dedent(comment)
                if autoformat_type == "markdown":
                    comment = markdown.markdown(comment)
            if len(helptext) > 0:
                helptext = textwrap.dedent(helptext)
                if autoformat_type == "markdown":
                    helptext = markdown.markdown(helptext)

        # Strip excess whitespace
        description = description.strip()
        comment = comment.strip()
        helptext = helptext.strip()

        self.sections.append(
            {
                "name": name,
                "anchor": anchor,
                "description": description,
                "comment": comment,
                "helptext": helptext,
                "plot": plot,
                "content": content,
                "print_section": any(
                    [n is not None and len(n) > 0 for n in [description, comment, helptext, plot, content]]
                ),
            }
        )

    def clean_s_name(self, s_name, f=None, root=None, filename=None, seach_pattern_key=None):
        """Helper function to take a long file name and strip it
        back to a clean sample name. Somewhat arbitrary.
        :param s_name: The sample name to clean
        :param root: The directory path that this file is within
        :config.prepend_dirs: boolean, whether to prepend dir name to s_name
        :return: The cleaned sample name, ready to be used
        """
        s_name_original = s_name

        # Backwards compatability - if f is a string, it's probably the root (this used to be the second argument)
        if isinstance(f, str):
            root = f
            f = None

        # Set string variables from f if it was a dict from find_log_files()
        if isinstance(f, dict):
            if "root" in f and root is None:
                root = f["root"]
            if "fn" in f and filename is None:
                filename = f["fn"]
            if "sp_key" in f and seach_pattern_key is None:
                seach_pattern_key = f["sp_key"]

        # For modules setting s_name from file contents, set s_name back to the filename
        # (if wanted in the config)
        if filename is not None and (
            config.use_filename_as_sample_name is True
            or (
                isinstance(config.use_filename_as_sample_name, list)
                and seach_pattern_key is not None
                and seach_pattern_key in config.use_filename_as_sample_name
            )
        ):
            s_name = filename

        # Set root to empty string if not known
        if root is None:
            root = ""

        # if s_name comes from file contents, it may have a file path
        # For consistency with other modules, we keep just the basename
        s_name = os.path.basename(s_name)

        # Prepend sample name with directory
        if config.prepend_dirs:
            sep = config.prepend_dirs_sep
            root = root.lstrip(".{}".format(os.sep))
            dirs = [d.strip() for d in root.split(os.sep) if d.strip() != ""]
            if config.prepend_dirs_depth != 0:
                d_idx = config.prepend_dirs_depth * -1
                if config.prepend_dirs_depth > 0:
                    dirs = dirs[d_idx:]
                else:
                    dirs = dirs[:d_idx]
            if len(dirs) > 0:
                s_name = "{}{}{}".format(sep.join(dirs), sep, s_name)

        if config.fn_clean_sample_names:
            # Split then take first section to remove everything after these matches
            for ext in config.fn_clean_exts:
                # Check if this config is limited to a module
                if "module" in ext:
                    if type(ext["module"]) is str:
                        ext["module"] = [ext["module"]]
                    if not any([m == self.anchor for m in ext["module"]]):
                        continue

                # Go through different filter types
                if type(ext) is str:
                    ext = {"type": "truncate", "pattern": ext}
                if ext.get("type") == "truncate":
                    s_name = s_name.split(ext["pattern"], 1)[0]
                elif ext.get("type") in ("remove", "replace"):
                    if ext["type"] == "replace":
                        logger.warning(
                            "use 'config.fn_clean_sample_names.remove' instead "
                            "of 'config.fn_clean_sample_names.replace' [deprecated]"
                        )
                    s_name = s_name.replace(ext["pattern"], "")
                elif ext.get("type") == "regex":
                    s_name = re.sub(ext["pattern"], "", s_name)
                elif ext.get("type") == "regex_keep":
                    match = re.search(ext["pattern"], s_name)
                    s_name = match.group() if match else s_name
                elif ext.get("type") is None:
                    logger.error('config.fn_clean_exts config was missing "type" key: {}'.format(ext))
                else:
                    logger.error("Unrecognised config.fn_clean_exts type: {}".format(ext.get("type")))
            # Trim off characters at the end of names
            for chrs in config.fn_clean_trim:
                if s_name.endswith(chrs):
                    s_name = s_name[: -len(chrs)]
                if s_name.startswith(chrs):
                    s_name = s_name[len(chrs) :]

        # Remove trailing whitespace
        s_name = s_name.strip()

        # If we cleaned back to an empty string, just use the original value
        if s_name == "":
            s_name = s_name_original

        # Do any hard replacements that are set with --replace-names
        if config.sample_names_replace:
            for s_name_search, s_name_replace in config.sample_names_replace.items():
                try:
                    # Skip if we're looking for exact matches only
                    if config.sample_names_replace_exact:
                        # Simple strings
                        if not config.sample_names_replace_regex and s_name != s_name_search:
                            continue
                        # regexes
                        if config.sample_names_replace_regex and not re.fullmatch(s_name_search, s_name):
                            continue
                    # Replace - regex
                    if config.sample_names_replace_regex:
                        s_name = re.sub(s_name_search, s_name_replace, s_name)
                    # Replace - simple string
                    else:
                        # Complete name swap
                        if config.sample_names_replace_complete:
                            if s_name_search in s_name:
                                s_name = s_name_replace
                        # Partial substring replace
                        else:
                            s_name = s_name.replace(s_name_search, s_name_replace)
                except re.error as e:
                    logger.error("Error with sample name replacement regex: {}".format(e))

        return s_name

    def ignore_samples(self, data):
        """Strip out samples which match `sample_names_ignore`"""
        try:
            if isinstance(data, OrderedDict):
                newdata = OrderedDict()
            elif isinstance(data, dict):
                newdata = dict()
            else:
                return data
            for s_name, v in data.items():
                if not self.is_ignore_sample(s_name):
                    newdata[s_name] = v
            return newdata
        except (TypeError, AttributeError):
            return data

    def is_ignore_sample(self, s_name):
        """Should a sample name be ignored?"""
        glob_match = any(fnmatch.fnmatch(s_name, sn) for sn in config.sample_names_ignore)
        re_match = any(re.match(sn, s_name) for sn in config.sample_names_ignore_re)
        return glob_match or re_match

    def general_stats_addcols(self, data, headers=None, namespace=None):
        """Helper function to add to the General Statistics variable.
        Adds to report.general_stats and does not return anything. Fills
        in required config variables if not supplied.
        :param data: A dict with the data. First key should be sample name,
                     then the data key, then the data.
        :param headers: Dict / OrderedDict with information for the headers,
                        such as colour scales, min and max values etc.
                        See docs/writing_python.md for more information.
        :return: None
        """
        if headers is None:
            headers = {}
        # Use the module namespace as the name if not supplied
        if namespace is None:
            namespace = self.name

        # Guess the column headers from the data if not supplied
        if headers is None or len(headers) == 0:
            hs = set()
            for d in data.values():
                hs.update(d.keys())
            hs = list(hs)
            hs.sort()
            headers = OrderedDict()
            for k in hs:
                headers[k] = dict()

        # Add the module name to the description if not already done
        keys = headers.keys()
        for k in keys:
            if "namespace" not in headers[k]:
                headers[k]["namespace"] = namespace
            if "description" not in headers[k]:
                headers[k]["description"] = headers[k].get("title", k)

        # Append to report.general_stats for later assembly into table
        report.general_stats_data.append(data)
        report.general_stats_headers.append(headers)

    def add_data_source(self, f=None, s_name=None, source=None, module=None, section=None):
        try:
            if module is None:
                module = self.name
            if section is None:
                section = "all_sections"
            if s_name is None:
                s_name = f["s_name"]
            if source is None:
                source = os.path.abspath(os.path.join(f["root"], f["fn"]))
            report.data_sources[module][section][s_name] = source
        except AttributeError:
            logger.warning("Tried to add data source for {}, but was missing fields data".format(self.name))

    def write_data_file(self, data, fn, sort_cols=False, data_format=None):
        """Saves raw data to a dictionary for downstream use, then redirects
        to report.write_data_file() to create the file in the report directory"""

        # Append custom module anchor if set
        mod_cust_config = getattr(self, "mod_cust_config", {})
        if "anchor" in mod_cust_config:
            fn = "{}_{}".format(fn, mod_cust_config["anchor"])

        # Generate a unique filename if the file already exists (running module multiple times)
        i = 1
        base_fn = fn
        while fn in report.saved_raw_data:
            fn = "{}_{}".format(base_fn, i)
            i += 1

        # Save the file
        report.saved_raw_data[fn] = data
        util_functions.write_data_file(data, fn, sort_cols, data_format)

    ##################################################
    #### DEPRECATED FORWARDERS
    def plot_bargraph(self, data, cats=None, pconfig=None):
        """Depreciated function. Forwards to new location."""
        from multiqc.plots import bargraph

        if pconfig is None:
            pconfig = {}
        return bargraph.plot(data, cats, pconfig)

    def plot_xy_data(self, data, pconfig=None):
        """Depreciated function. Forwards to new location."""
        from multiqc.plots import linegraph

        if pconfig is None:
            pconfig = {}
        return linegraph.plot(data, pconfig)
