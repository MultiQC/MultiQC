"""
MultiQC modules base class, contains helper functions
"""

import dataclasses
from pathlib import Path
from typing import List, Union, Optional, Dict, Any, cast, Tuple

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
import packaging.version

from multiqc.plots.plotly.plot import Plot
from multiqc import config, report
from multiqc.core import software_versions

logger = logging.getLogger(__name__)


class ModuleNoSamplesFound(Exception):
    """Module checked all input files but couldn't find any data to use"""


@dataclasses.dataclass
class Section:
    name: str
    anchor: str
    description: str
    module: str
    comment: str = ""
    helptext: str = ""
    content_before_plot: str = ""
    content: str = ""
    plot: str = ""
    print_section: bool = True
    plot_id: Optional[str] = None


class BaseMultiqcModule:
    # Custom options from user config that can overwrite base module values
    mod_cust_config: Dict = {}
    mod_id = None

    def __init__(
        self,
        name="base",
        anchor="base",
        target=None,
        href: Union[str, List[str], None] = None,
        info=None,
        comment=None,
        extra=None,
        autoformat=True,
        autoformat_type="markdown",
        doi: Optional[Union[str, List[str]]] = None,
    ):
        # Custom options from user config that can overwrite base module values
        self.name = self.mod_cust_config.get("name", name)
        self.id = self.mod_id if self.mod_id else anchor  # cannot be overwritten for repeated modules with path_filters
        self.anchor = self.mod_cust_config.get("anchor", anchor)
        self.href = self.mod_cust_config.get("href", [href] if isinstance(href, str) else href or [])
        self.info = self.mod_cust_config.get("info", info)
        self.comment = self.mod_cust_config.get("comment", comment)
        self.extra = self.mod_cust_config.get("extra", extra)
        self.doi = self.mod_cust_config.get("doi", [doi] if isinstance(doi, str) else doi or [])

        # List of software version(s) for module. Don't append directly, use add_software_version()
        self.versions: Dict[str, List[Tuple[Optional[packaging.version.Version], str]]] = defaultdict(list)

        # Specific module level config to overwrite (e.g. config.bcftools, config.fastqc)
        config.update({self.id: self.mod_cust_config.get("custom_config", {})})

        # Sanitise anchor ID and check for duplicates
        self.anchor = report.save_htmlid(self.anchor)

        # See if we have a user comment in the config
        if self.anchor in config.section_comments:
            self.comment = config.section_comments[self.anchor]

        if self.info is None:
            self.info = ""
        self.info = self.info.strip().strip(".")
        # Legacy: if self.info starts with a lowercase letter, prepend the module name to it
        if self.info and self.info[0].islower():
            self.info = f"{self.name} {self.info}"

        if self.extra is None:
            self.extra = ""

        if isinstance(self.href, str):
            self.href = [self.href]
        self.href = [i for i in self.href if i != ""]

        if isinstance(self.doi, str):
            self.doi = [self.doi]
        self.doi = [i for i in self.doi if i != ""]

        self.intro = self._get_intro()

        # Format the markdown strings
        if autoformat:
            if self.comment is not None:
                self.comment = textwrap.dedent(self.comment)
                if autoformat_type == "markdown":
                    self.comment = markdown.markdown(self.comment)

        self.sections: List[Section] = []

        self.hidden = False

        self.__saved_raw_data: Dict[str, Dict[str, Any]] = dict()  # Saved raw data. Identical to report.saved_raw_data

        self.css: Dict[str, str] = dict()
        self.js: Dict[str, str] = dict()

        # Get list of all base attributes, so we clean up any added by child modules
        self._base_attributes = [k for k in dir(self)]

    def _get_intro(self):
        doi_html = ""
        if len(self.doi) > 0:
            doi_links = []
            for doi in self.doi:
                # Build the HTML link for the DOI
                doi_links.append(
                    f' <a class="module-doi" data-doi="{doi}" data-toggle="popover" href="https://doi.org/{doi}" target="_blank">{doi}</a>'
                )
            doi_html = '<em class="text-muted small" style="margin-left: 1rem;">DOI: {}</em>'.format(
                "; ".join(doi_links)
            )

        url_link = ""
        if len(self.href) > 0:
            url_links = []
            for url in self.href:
                url_links.append(f'<a href="{url}" target="_blank">{url.strip("/")}</a>')
            url_link = '<em class="text-muted small" style="margin-left: 1rem;">URL: {}</em>'.format(
                "; ".join(url_links)
            )

        info = (self.info + ".") if self.info else ""
        return f"<p>{info}{url_link}{doi_html}</p>{self.extra}"

    def clean_child_attributes(self):
        """
        Clean up non-base attribute to save memory. If the attribute is added in subclass,
        but absent in the base class BaseMultiqcModule, delete it.
        """
        for key in list(self.__dict__.keys()):
            if key not in self._base_attributes and not key.startswith("_"):
                logger.debug(f"{self.anchor}: deleting attribute self.{key}")
                delattr(self, key)

    @property
    def saved_raw_data(self):
        """
        Wrapper to give access to private __saved_raw_data. We could have just called __saved_raw_data without the
        underscore: saved_raw_data, and that would work just fine. But users might override saved_raw_data in
        their child modules, and we would lose that data.
        """
        return self.__saved_raw_data

    def find_log_files(self, sp_key, filecontents=True, filehandles=False):
        """
        Return matches log files of interest.
        :param sp_key: Search pattern key specified in config
        :param filecontents: f["f"] will contain raw file contents
        :param filehandles: f["f"] will be the file handle
        :return: Yields a dict with filename (fn), root directory (root), cleaned sample name
                 generated from the filename (s_name) and either the file contents or file handle
                 for the current matched file (f).
                 As yield is used, the results can be iterated over without loading all files at once
        """

        # Pick up path filters if specified.
        # Allows modules to be called multiple times with different sets of files
        def get_path_filters(key: str) -> List[str]:
            pfs: List[str] = []
            val = self.mod_cust_config.get(key, [])
            values = val if isinstance(val, list) else [val]
            pf: str
            for pf in values:
                if pf.startswith("./"):
                    pf = pf[2:]
                pfs.append(pf)
            return pfs

        path_filters: List[str] = get_path_filters("path_filters")
        path_filters_exclude: List[str] = get_path_filters("path_filters_exclude")

        if not isinstance(sp_key, str):
            logger.warning(f"The find_log_files() search key must be a string, got {type(sp_key)}: {sp_key}")
            return

        for f in report.files.get(sp_key, []):
            # Make a note of the filename so that we can report it if something crashes
            last_found_file: str = os.path.join(f["root"], f["fn"])
            report.last_found_file = last_found_file

            # Filter out files based on exclusion patterns
            if path_filters_exclude and len(path_filters_exclude) > 0:
                # Try both the given path and also the path prefixed with the analysis dirs
                exclusion_hits = itertools.chain(
                    (fnmatch.fnmatch(last_found_file, pfe) for pfe in path_filters_exclude),
                    *(
                        (
                            fnmatch.fnmatch(last_found_file, os.path.join(analysis_dir, pfe))
                            for pfe in path_filters_exclude
                        )
                        for analysis_dir in report.analysis_files
                    ),
                )
                if any(exclusion_hits):
                    logger.debug(
                        f"{sp_key} - Skipping '{report.last_found_file}' as it matched the path_filters_exclude for '{self.name}'"
                    )
                    continue

            # Filter out files based on inclusion patterns
            if path_filters and len(path_filters) > 0:
                # Try both the given path and also the path prefixed with the analyis dirs
                inclusion_hits = itertools.chain(
                    (fnmatch.fnmatch(last_found_file, pf) for pf in path_filters),
                    *(
                        (fnmatch.fnmatch(last_found_file, os.path.join(analysis_dir, pf)) for pf in path_filters)
                        for analysis_dir in report.analysis_files
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
                    fh: io.IOBase  # make mypy happy
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
        content_before_plot="",
        plot: Optional[Union[Plot, str]] = None,
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

        section = Section(
            name=name,
            anchor=anchor,
            description=description,
            module=self.name,
            comment=comment,
            helptext=helptext,
            content_before_plot=content_before_plot,
            content=content,
            print_section=any([description, comment, helptext, content_before_plot, plot, content]),
        )

        if plot is not None:
            if isinstance(plot, Plot):
                section.plot_id = plot.id
                # separately keeping track of Plot objects to be rendered further
                report.plot_by_id[plot.id] = plot
            elif isinstance(plot, str):
                section.plot = plot

        # self.sections is passed into Jinja template:
        self.sections.append(section)

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

    def clean_s_name(
        self,
        s_name: Union[str, List[str]],
        f: Optional[Union[Dict, str]] = None,
        root: Optional[str] = None,
        filename: Optional[str] = None,
        search_pattern_key: Optional[str] = None,
    ) -> str:
        """
        Helper function to take a long file name(s) and strip back to one clean sample name. Somewhat arbitrary.
        """
        if isinstance(s_name, list):
            if len(s_name) == 0:
                raise ValueError("Empty list of sample names passed to clean_s_name()")

            # Extract a sample name from a list of file names (for example, FASTQ pairs).
            # Each name is cleaned separately first:
            clean_names = [
                self.clean_s_name(sn, f=f, root=root, filename=filename, search_pattern_key=search_pattern_key)
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

        sname: str = cast(str, s_name)
        sname_original = sname

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
            if "sp_key" in f and search_pattern_key is None:
                search_pattern_key = f["sp_key"]

        # For modules setting s_name from file contents, set s_name back to the filename
        # (if wanted in the config)
        if filename is not None and (
            config.use_filename_as_sample_name is True
            or (
                isinstance(config.use_filename_as_sample_name, list)
                and search_pattern_key is not None
                and search_pattern_key in config.use_filename_as_sample_name
            )
        ):
            sname = filename

        # if s_name comes from file contents, it may have a file path
        # For consistency with other modules, we keep just the basename
        sname = os.path.basename(sname)

        # Prepend sample name with directory
        if config.prepend_dirs:
            sep = config.prepend_dirs_sep
            dirs = [d.strip() for d in (Path(root).parts if root else []) if d.strip() != ""]
            if config.prepend_dirs_depth != 0:
                d_idx = config.prepend_dirs_depth * -1
                if config.prepend_dirs_depth > 0:
                    dirs = dirs[d_idx:]
                else:
                    dirs = dirs[:d_idx]
            if len(dirs) > 0:
                sname = f"{sep.join(dirs)}{sep}{sname}"

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
                    sname = sname.split(ext["pattern"], 1)[0]
                elif ext.get("type") in ("remove", "replace"):
                    if ext["type"] == "replace":
                        logger.warning(
                            "use 'config.fn_clean_sample_names.remove' instead "
                            "of 'config.fn_clean_sample_names.replace' [deprecated]"
                        )
                    sname = sname.replace(ext["pattern"], "")
                elif ext.get("type") == "regex":
                    sname = re.sub(ext["pattern"], "", sname)
                elif ext.get("type") == "regex_keep":
                    match = re.search(ext["pattern"], sname)
                    sname = match.group() if match else sname
                elif ext.get("type") is None:
                    logger.error(f'config.fn_clean_exts config was missing "type" key: {ext}')
                else:
                    logger.error(f"Unrecognised config.fn_clean_exts type: {ext.get('type')}")
            # Trim off characters at the end of names
            for chrs in config.fn_clean_trim:
                if sname.endswith(chrs):
                    sname = sname[: -len(chrs)]
                if sname.startswith(chrs):
                    sname = sname[len(chrs) :]

        # Remove trailing whitespace
        sname = sname.strip()

        # If we cleaned back to an empty string, just use the original value
        if sname == "":
            sname = sname_original

        # Do any hard replacements that are set with --replace-names
        if config.sample_names_replace:
            for s_name_search, s_name_replace in config.sample_names_replace.items():
                try:
                    # Skip if we're looking for exact matches only
                    if config.sample_names_replace_exact:
                        # Simple strings
                        if not config.sample_names_replace_regex and sname != s_name_search:
                            continue
                        # regexes
                        if config.sample_names_replace_regex and not re.fullmatch(s_name_search, sname):
                            continue
                    # Replace - regex
                    if config.sample_names_replace_regex:
                        sname = re.sub(s_name_search, s_name_replace, sname)
                    # Replace - simple string
                    else:
                        # Complete name swap
                        if config.sample_names_replace_complete:
                            if s_name_search in sname:
                                sname = s_name_replace
                        # Partial substring replace
                        else:
                            sname = sname.replace(s_name_search, s_name_replace)
                except re.error as e:
                    logger.error(f"Error with sample name replacement regex: {e}")

        return sname

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
            headers = dict()
            for k in sorted(hs):
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
        if s_name is not None and self.is_ignore_sample(s_name):
            return
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

    def add_software_version(
        self, version: Optional[str] = None, sample: Optional[str] = None, software_name: Optional[str] = None
    ):
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
        parsed_version = software_versions.parse_version(version)
        if (parsed_version, version) in self.versions[software_name]:
            return

        self.versions[software_name].append((parsed_version, version))

        # Sort version in order newest --> oldest
        self.versions[software_name] = software_versions.sort_versions(self.versions[software_name])

        # Update version list for report section.
        group_name = self.name
        report.software_versions[group_name][software_name] = [v for _, v in self.versions[software_name]]

    def write_data_file(self, data, fn, sort_cols=False, data_format=None):
        """Saves raw data to a dictionary for downstream use, then redirects
        to report.write_data_file() to create the file in the report directory"""

        # Append custom module anchor if set
        if self.mod_cust_config.get("anchor"):
            fn = f"{fn}_{self.mod_cust_config['anchor']}"

        # Generate a unique filename if the file already exists (running module multiple times)
        i = 1
        base_fn = fn
        while fn in report.saved_raw_data:
            fn = f"{base_fn}_{i}"
            i += 1

        if config.preserve_module_raw_data:
            report.saved_raw_data[fn] = data
            # Keep also in the module instance, so it's possible to map back data to specific module
            self.__saved_raw_data[fn] = data

        # Save the file
        report.write_data_file(data, fn, sort_cols, data_format)

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
