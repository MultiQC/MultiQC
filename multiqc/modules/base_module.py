""" MultiQC modules base class, contains helper functions """

from typing import List, Union, Optional, Dict

import fnmatch
import io
import itertools
import logging
import mimetypes
import os
import re
import textwrap
from collections import defaultdict

import markdown

from multiqc.utils import config, report, software_versions, util_functions

logger = logging.getLogger(__name__)


class ModuleNoSamplesFound(Exception):
    """Module checked all input files but couldn't find any data to use"""


class BaseMultiqcModule:
    # Custom options from user config that can overwrite base module values
    mod_cust_config: Dict = {}

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
        doi: Optional[str] = None,
    ):
        # Custom options from user config that can overwrite base module values
        self.name = self.mod_cust_config.get("name", name)
        self.anchor = self.mod_cust_config.get("anchor", anchor)
        target = self.mod_cust_config.get("target", target)
        self.href = self.mod_cust_config.get("href", href)
        self.info = self.mod_cust_config.get("info", info)
        self.comment = self.mod_cust_config.get("comment", comment)
        self.extra = self.mod_cust_config.get("extra", extra)
        self.doi = self.mod_cust_config.get("doi", (doi or []))

        # List of software version(s) for module. Don't append directly, use add_software_version()
        self.versions = defaultdict(list)

        # Specific module level config to overwrite (e.g. config.bcftools, config.fastqc)
        config.update({anchor: self.mod_cust_config.get("custom_config", {})})

        # Sanitise anchor ID and check for duplicates
        self.anchor = report.save_htmlid(self.anchor)

        # See if we have a user comment in the config
        if self.anchor in config.section_comments:
            self.comment = config.section_comments[self.anchor]

        if self.info is None:
            self.info = ""
        # Always finish with a ".", as we may add a DOI after the intro.
        if len(self.info) > 0 and self.info[-1] != ".":
            self.info += "."
        if self.extra is None:
            self.extra = ""
        self.doi_link = ""
        if isinstance(self.doi, str):
            self.doi = [self.doi]
        self.doi = [i for i in self.doi if i != ""]
        if len(self.doi) > 0:
            doi_links = []
            for doi in self.doi:
                # Build the HTML link for the DOI
                doi_links.append(
                    f' <a class="module-doi" data-doi="{doi}" data-toggle="popover" href="https://doi.org/{doi}" target="_blank">{doi}</a>'
                )
            self.doi_link = '<em class="text-muted small" style="margin-left: 1rem;">DOI: {}.</em>'.format(
                "; ".join(doi_links)
            )

        if target is None:
            target = self.name
        if self.href is not None:
            self.mname = f'<a href="{self.href}" target="_blank">{target}</a>'
        else:
            self.mname = target
        if self.href or self.info or self.extra or self.doi_link:
            self.intro = f"<p>{self.mname} {self.info}{self.doi_link}</p>{self.extra}"

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
        path_filters: Union[str, List[str]] = self.mod_cust_config.get("path_filters", [])
        path_filters_exclude: Union[str, List[str]] = self.mod_cust_config.get("path_filters_exclude", [])
        if isinstance(path_filters, str):
            path_filters: List[str] = [path_filters]
        if isinstance(path_filters_exclude, str):
            path_filters_exclude: List[str] = [path_filters_exclude]

        # Old, depreciated syntax support. Likely to be removed in a future version.
        if isinstance(sp_key, dict):
            report.files[self.name] = list()
            for sf in report.searchfiles:
                if report.search_file(sp_key, {"fn": sf[0], "root": sf[1]}, module_key=None):
                    report.files[self.name].append({"fn": sf[0], "root": sf[1]})
            sp_key = self.name
            logwarn = f"Depreciation Warning: {self.name} - Please use new style for find_log_files()"
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
                # Try both the given path and also the path prefixed with the analysis dirs
                exlusion_hits = itertools.chain(
                    (fnmatch.fnmatch(report.last_found_file, pfe) for pfe in path_filters_exclude),
                    *(
                        (
                            fnmatch.fnmatch(report.last_found_file, os.path.join(analysis_dir, pfe))
                            for pfe in path_filters_exclude
                        )
                        for analysis_dir in config.analysis_dir
                    ),
                )
                if any(exlusion_hits):
                    logger.debug(
                        f"{sp_key} - Skipping '{report.last_found_file}' as it matched the path_filters_exclude for '{self.name}'"
                    )
                    continue

            # Filter out files based on inclusion patterns
            if path_filters and len(path_filters) > 0:
                # Try both the given path and also the path prefixed with the analyis dirs
                inclusion_hits = itertools.chain(
                    (fnmatch.fnmatch(report.last_found_file, pf) for pf in path_filters),
                    *(
                        (fnmatch.fnmatch(report.last_found_file, os.path.join(analysis_dir, pf)) for pf in path_filters)
                        for analysis_dir in config.analysis_dir
                    ),
                )
                if not any(inclusion_hits):
                    logger.debug(
                        f"{sp_key} - Skipping '{report.last_found_file}' as it didn't match the path_filters for '{self.name}'"
                    )
                    continue
                else:
                    logger.debug(
                        f"{sp_key} - Selecting '{report.last_found_file}' as it matched the path_filters for '{self.name}'"
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
                                try:
                                    f["f"] = fh.read()
                                except UnicodeDecodeError as e:
                                    logger.debug(
                                        f"Couldn't read file as utf-8: {f['fn']}, will attempt to skip non-unicode characters\n{e}"
                                    )
                                    try:
                                        with io.open(
                                            os.path.join(f["root"], f["fn"]), "r", encoding="utf-8", errors="ignore"
                                        ) as fh_ignoring:
                                            f["f"] = fh_ignoring.read()
                                    except Exception as e:
                                        logger.debug(f"Still couldn't read file: {f['fn']}\n{e}")
                                        f["f"] = None
                                    finally:
                                        fh.close()
                                yield f
                except (IOError, OSError, ValueError, UnicodeDecodeError) as e:
                    logger.debug(f"Couldn't open filehandle when returning file: {f['fn']}\n{e}")
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
                anchor = f"{self.anchor}-{nid}"
            else:
                sl = len(self.sections) + 1
                anchor = f"{self.anchor}-section-{sl}"

        # Append custom module anchor to the section if set
        if "anchor" in self.mod_cust_config:
            anchor = f"{self.mod_cust_config['anchor']}_{anchor}"

        # Sanitise anchor ID and check for duplicates
        anchor = report.save_htmlid(anchor)

        # Skip if user has a config to remove this module section
        if anchor in config.remove_sections:
            logger.debug(f"Skipping section '{anchor}' because specified in user config")
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

    @staticmethod
    def _clean_fastq_pair(r1: str, r2: str) -> Optional[str]:
        """
        Try trimming r1 and r2 as paired FASTQ file names.
        """
        # Try trimming the conventional illumina suffix with a tail 001 ending. Refs:
        # https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
        # https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/fastq-input#:~:text=10x%20pipelines%20need%20files%20named,individual%20who%20demultiplexed%20your%20flowcell.
        cleaned_r1 = re.sub(r"_R1_\d{3}$", "", r1)
        cleaned_r2 = re.sub(r"_R2_\d{3}$", "", r2)
        if cleaned_r1 == cleaned_r2:  # trimmed successfully
            return cleaned_r1

        # Try removing _R1 and _R2 from the middle.
        cleaned_r1 = re.sub(r"_R1_", "_", r1)
        cleaned_r2 = re.sub(r"_R2_", "_", r2)
        if cleaned_r1 == cleaned_r2:  # trimmed successfully
            return cleaned_r1

        # Try trimming other variations from the end (-R1, _r1, _1, .1, etc).
        cleaned_r1 = re.sub(r"([_.-][rR]?1)?$", "", r1)
        cleaned_r2 = re.sub(r"([_.-][rR]?2)?$", "", r2)
        if cleaned_r1 == cleaned_r2:  # trimmed successfully
            return cleaned_r1

        return None

    def clean_s_name(self, s_name: Union[str, List[str]], f=None, root=None, filename=None, seach_pattern_key=None):
        """
        Helper function to take a long file name(s) and strip back to one clean sample name. Somewhat arbitrary.
        :param s_name: The sample name(s) to clean.
        :param root: The directory path that this file is within
        :config.prepend_dirs: boolean, whether to prepend dir name to s_name
        :return: The cleaned sample name, ready to be used
        """
        if isinstance(s_name, list):
            if len(s_name) == 0:
                raise ValueError("Empty list of sample names passed to clean_s_name()")

            # Extract a sample name from a list of file names (for example, FASTQ pairs).
            # Each name is cleaned separately first:
            clean_names = [
                self.clean_s_name(sn, f=f, root=root, filename=filename, seach_pattern_key=seach_pattern_key)
                for sn in s_name
            ]
            if len(set(clean_names)) == 1:
                # All the same, returning the first one.
                return clean_names[0]

            if len(clean_names) == 2:
                # Checking if it's a FASTQ pair.
                fastq_s_name = self._clean_fastq_pair(*clean_names)
                if fastq_s_name is not None:
                    return fastq_s_name

            # Couldn't clean as FASTQ. Just concatenating the clean names.
            return "_".join(clean_names)

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
            root = root.lstrip(f".{os.sep}")
            dirs = [d.strip() for d in root.split(os.sep) if d.strip() != ""]
            if config.prepend_dirs_depth != 0:
                d_idx = config.prepend_dirs_depth * -1
                if config.prepend_dirs_depth > 0:
                    dirs = dirs[d_idx:]
                else:
                    dirs = dirs[:d_idx]
            if len(dirs) > 0:
                s_name = f"{sep.join(dirs)}{sep}{s_name}"

        if config.fn_clean_sample_names:
            # Split then take first section to remove everything after these matches
            for ext in config.fn_clean_exts:
                # Check if this config is limited to a module
                if "module" in ext:
                    if isinstance(ext["module"], str):
                        ext["module"] = [ext["module"]]
                    if not any([m == self.anchor for m in ext["module"]]):
                        continue

                # Go through different filter types
                if isinstance(ext, str):
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
                    logger.error(f'config.fn_clean_exts config was missing "type" key: {ext}')
                else:
                    logger.error(f"Unrecognised config.fn_clean_exts type: {ext.get('type')}")
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
                    logger.error(f"Error with sample name replacement regex: {e}")

        return s_name

    def ignore_samples(self, data):
        """Strip out samples which match `sample_names_ignore`"""
        try:
            if isinstance(data, dict):
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
        :param headers: Dict with information for the headers,
                        such as colour scales, min and max values etc.
                        See docs/writing_python.md for more information.
        :param namespace: Append to the module name in the table column description.
                          Can be e.g. a submodule name.
        :return: None
        """
        if headers is None:
            headers = {}
        # Deepish copy of headers so that we can modify it in place
        headers = {k: v.copy() for k, v in headers.items()}

        # Guess the column headers from the data if not supplied
        if headers is None or len(headers) == 0:
            hs = set()
            for d in data.values():
                hs.update(d.keys())
            hs = list(hs)
            hs.sort()
            headers = dict()
            for k in hs:
                headers[k] = dict()

        # Add the module name to the description if not already done
        keys = headers.keys()
        for k in keys:
            # Prepend the namespace displayed in the table with the module name
            namespace = headers[k].get("namespace", namespace)
            headers[k]["namespace"] = self.name
            if namespace:
                headers[k]["namespace"] = self.name + ": " + namespace
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
            logger.warning(f"Tried to add data source for {self.name}, but was missing fields data")

    def add_software_version(self, version: str = None, sample: str = None, software_name: str = None):
        """Save software versions for module."""
        # Don't add if version is None. This allows every module to call this function
        # even those without a version to add. This is useful to check that all modules
        # are calling this function.
        if version is None:
            return

        # Don't add if version detection is disabled
        if config.disable_version_detection:
            return

        # Don't add if sample is ignored
        if sample is not None and self.is_ignore_sample(sample):
            return

        # Use module name as software name if not specified
        if software_name is None:
            software_name = self.name

        # Check if version string is PEP 440 compliant to enable version normalization and proper ordering.
        # Otherwise, use raw string is used for version.
        # - https://peps.python.org/pep-0440/
        version = software_versions.parse_version(version)

        if version in self.versions[software_name]:
            return

        self.versions[software_name].append(version)

        # Sort version in order newest --> oldest
        self.versions[software_name] = software_versions.sort_versions(self.versions[software_name])

        # Update version list for report section.
        group_name = self.name
        report.software_versions[group_name][software_name] = self.versions[software_name]

    def write_data_file(self, data, fn, sort_cols=False, data_format=None):
        """Saves raw data to a dictionary for downstream use, then redirects
        to report.write_data_file() to create the file in the report directory"""

        # Append custom module anchor if set
        if "anchor" in self.mod_cust_config:
            fn = f"{fn}_{self.mod_cust_config['anchor']}"

        # Generate a unique filename if the file already exists (running module multiple times)
        i = 1
        base_fn = fn
        while fn in report.saved_raw_data:
            fn = f"{base_fn}_{i}"
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
