"""
MultiQC modules base class, contains helper functions
"""

import dataclasses
import fnmatch
import io
import itertools
import logging
import mimetypes
import os
import re
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Literal,
    Optional,
    Set,
    Tuple,
    TypeVar,
    Union,
    cast,
    overload,
)

import markdown
import packaging.version

from multiqc import config, report, validation
from multiqc.config import CleanPatternT
from multiqc.core import software_versions
from multiqc.core.strict_helpers import lint_error
from multiqc.plots.plotly.plot import Plot
from multiqc.plots.table_object import (
    ColumnDict,
    ColumnKey,
    InputRow,
    SampleGroup,
    SampleName,
    ValueT,
)
from multiqc.types import Anchor, FileDict, LoadedFileDict, ModuleId, SectionId

logger = logging.getLogger(__name__)


class ModuleNoSamplesFound(Exception):
    """Module checked all input files but couldn't find any data to use"""


@dataclasses.dataclass
class SampleNameMeta:
    original_name: SampleName
    trimmed_name: Optional[SampleName] = None
    trimmed_suffixes: List[str] = dataclasses.field(default_factory=list)
    group: Optional[SampleGroup] = None
    labels: List[str] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class Section:
    name: str
    anchor: Anchor
    id: SectionId  # unlike anchor, doesn't have to be different from the module or plot ids
    description: str
    module: str
    comment: str = ""
    helptext: str = ""
    content_before_plot: str = ""
    content: str = ""
    plot: str = ""
    print_section: bool = True
    plot_anchor: Optional[Anchor] = None


ExtraFunctionType = Callable[[InputRow, List[Tuple[Optional[str], SampleName, SampleName]]], None]

DataT = TypeVar("DataT")
SampleNameT = TypeVar("SampleNameT", str, SampleName)


@dataclasses.dataclass
class SampleGroupingConfig:
    cols_to_weighted_average: Optional[List[Tuple[ColumnKey, ColumnKey]]] = None
    cols_to_average: Optional[List[ColumnKey]] = None
    cols_to_sum: Optional[List[ColumnKey]] = None
    extra_functions: Optional[List[ExtraFunctionType]] = dataclasses.field(default_factory=list)


class BaseMultiqcModule:
    # Custom options from user config that can overwrite base module values
    mod_cust_config: Dict[str, Any] = {}
    mod_id: Optional[ModuleId] = None

    def __init__(
        self,
        name: str = "base",
        anchor: Union[Anchor, str] = Anchor("base"),
        target: Optional[str] = None,
        href: Union[str, List[str], None] = None,
        info: Optional[str] = None,
        comment: Optional[str] = None,
        extra: Optional[str] = None,
        autoformat: bool = True,
        autoformat_type: str = "markdown",
        doi: Optional[Union[str, List[str]]] = None,
    ):
        validation.reset()

        # Custom options from user config that can overwrite base module values
        self.name: str = name
        _cust_name = self.mod_cust_config.get("name")
        if _cust_name is not None:
            self.name = str(_cust_name)

        # cannot be overwritten for repeated modules with path_filters:
        self.id: ModuleId = ModuleId(self.mod_id or anchor)

        self.anchor: Anchor = Anchor(anchor)
        _cust_anchor = self.mod_cust_config.get("anchor")
        if _cust_anchor is not None:
            self.anchor = Anchor(str(_cust_anchor))

        self.info: str = info or ""
        _cust_info = self.mod_cust_config.get("info")
        if _cust_info is not None:
            self.info = str(_cust_info)

        self.comment: str = comment or ""
        _cust_comment = self.mod_cust_config.get("comment")
        if _cust_comment is not None:
            self.comment = str(_cust_comment)

        self.extra: str = extra or ""
        _cust_extra = self.mod_cust_config.get("extra")
        if _cust_extra is not None:
            self.extra = str(_cust_extra)

        self.href: List[str] = [href] if isinstance(href, str) else href or []
        _cust_href = self.mod_cust_config.get("href")
        if _cust_href is not None:
            if isinstance(_cust_href, str):
                self.href = [_cust_href]
            elif isinstance(_cust_href, list):
                self.href = [str(h) for h in _cust_href]

        self.doi: List[str] = [doi] if isinstance(doi, str) else doi or []
        _cust_doi = self.mod_cust_config.get("doi")
        if _cust_doi is not None:
            if isinstance(_cust_doi, str):
                self.doi = [_cust_doi]
            elif isinstance(_cust_doi, list):
                self.doi = [str(d) for d in _cust_doi]

        self.skip_generalstats = True if self.mod_cust_config.get("generalstats") is False else False

        # List of software version(s) for module. Don't append directly, use add_software_version()
        self.versions: Dict[str, List[Tuple[Optional[packaging.version.Version], str]]] = defaultdict(list)

        # Specific module level config to overwrite (e.g. config.bcftools, config.fastqc)
        config.update({self.id: self.mod_cust_config.get("custom_config", {})})

        # Sanitise anchor ID and check for duplicates
        self.anchor = Anchor(report.save_htmlid(str(self.anchor)))

        # See if we have a user comment in the config
        if _config_section_comment := config.section_comments.get(str(self.anchor)):
            self.comment = _config_section_comment

        self.info = self.info.strip().strip(".")
        # Legacy: if self.info starts with a lowercase letter, prepend the module name to it
        if self.info and self.info[0].islower():
            self.info = f"{self.name} {self.info}"

        if isinstance(self.href, str):
            self.href = [self.href]
        self.href = [i for i in self.href if i != ""]

        if isinstance(self.doi, str):
            self.doi = [self.doi]
        self.doi = [i for i in self.doi if i != ""]

        self.intro = self._get_intro()

        # Format the markdown strings
        if autoformat and self.comment:
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

        self.sample_names: List[SampleNameMeta] = []

    def _get_intro(self):
        doi_html = ""
        if len(self.doi) > 0:
            doi_links: List[str] = []
            for doi in self.doi:
                # Build the HTML link for the DOI
                doi_links.append(
                    f' <a class="module-doi" data-doi="{doi}" data-toggle="popover" href="https://doi.org/{doi}" '
                    f'target="_blank">{doi}</a>'
                )
            doi_html = '<em class="text-muted small" style="margin-left: 1rem;">DOI: {}</em>'.format(
                "; ".join(doi_links)
            )

        url_link = ""
        if len(self.href) > 0:
            url_links: List[str] = []
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

    @overload
    def find_log_files(
        self, sp_key: str, filecontents: Literal[True] = True, filehandles: bool = False
    ) -> Iterable[LoadedFileDict[str]]: ...

    @overload
    def find_log_files(
        self, sp_key: str, filecontents: Literal[False] = False, filehandles: Literal[True] = True
    ) -> Iterable[LoadedFileDict[io.BufferedReader]]: ...

    @overload
    def find_log_files(
        self, sp_key: str, filecontents: Literal[False] = False, filehandles: Literal[False] = False
    ) -> Iterable[LoadedFileDict[None]]: ...

    def find_log_files(
        self, sp_key: str, filecontents: bool = True, filehandles: bool = False
    ) -> Union[
        Iterable[LoadedFileDict[str]],
        Iterable[LoadedFileDict[io.BufferedReader]],
        Iterable[LoadedFileDict[io.TextIOWrapper]],
        Iterable[LoadedFileDict[None]],
    ]:
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
            values: List[str] = val if isinstance(val, list) else [val]
            pf: str
            for pf in values:
                if pf.startswith("./"):
                    pf = pf[2:]
                pfs.append(pf)
            return pfs

        path_filters: List[str] = get_path_filters("path_filters")
        path_filters_exclude: List[str] = get_path_filters("path_filters_exclude")

        for f in report.files.get(ModuleId(sp_key), []):
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
                        f"{sp_key} - Skipping '{report.last_found_file}' as it matched the path_filters_exclude for "
                        f"'{self.name}'"
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
                        f"{sp_key} - Skipping '{report.last_found_file}' as it didn't match the path_filters for '"
                        f"{self.name}'"
                    )
                    continue
                else:
                    logger.debug(
                        f"{sp_key} - Selecting '{report.last_found_file}' as it matched the path_filters for '"
                        f"{self.name}'"
                    )

            # Make a sample name from the filename
            s_name = self.clean_s_name(f["fn"], f)

            if filehandles or filecontents:
                try:
                    fh: Union[io.BufferedReader, io.TextIOWrapper, None]
                    # Custom content module can now handle image files
                    (ftype, _) = mimetypes.guess_type(os.path.join(f["root"], f["fn"]))
                    if ftype is not None and ftype.startswith("image"):
                        with io.open(os.path.join(f["root"], f["fn"]), "rb") as fh:
                            # always return file handles
                            yield {**f, "s_name": s_name, "f": fh}
                    else:
                        # Everything else - should be all text files
                        with io.open(os.path.join(f["root"], f["fn"]), "r", encoding="utf-8") as fh:
                            if filehandles:
                                yield {**f, "s_name": s_name, "f": fh}
                            elif filecontents:
                                try:
                                    contents = fh.read()
                                except UnicodeDecodeError as e:
                                    logger.debug(
                                        f"Couldn't read file as utf-8: {f['fn']}, will attempt to skip non-unicode "
                                        f"characters\n{e}"
                                    )
                                    try:
                                        with io.open(
                                            os.path.join(f["root"], f["fn"]),
                                            "r",
                                            encoding="utf-8",
                                            errors="ignore",
                                        ) as fh_ignoring:
                                            yield {**f, "s_name": s_name, "f": fh_ignoring.read()}
                                    except Exception as e:
                                        logger.debug(f"Still couldn't read file: {f['fn']}\n{e}")
                                        yield {**f, "s_name": s_name, "f": None}
                                    finally:
                                        fh.close()
                                else:
                                    yield {**f, "s_name": s_name, "f": str(contents)}
                except (IOError, OSError, ValueError, UnicodeDecodeError) as e:
                    logger.debug(f"Couldn't open filehandle when returning file: {f['fn']}\n{e}")
                    yield {**f, "s_name": s_name, "f": None}
            else:
                yield {**f, "s_name": s_name, "f": None}

    def add_section(
        self,
        name: Optional[str] = None,
        anchor: Optional[Union[str, Anchor]] = None,
        id: Optional[Union[str, SectionId]] = None,
        description: str = "",
        comment: str = "",
        helptext: str = "",
        content_before_plot: str = "",
        plot: Optional[Union[Plot[Any, Any], str]] = None,
        content: str = "",
        autoformat: bool = True,
        autoformat_type: str = "markdown",
    ):
        """Add a section to the module report output"""
        if id is None and anchor is not None:
            id = str(anchor)

        if anchor is None and id is not None:
            anchor = str(id)

        if id is None:
            if name is not None:
                nid = name.lower().strip().replace(" ", "-")
                id = f"{self.anchor}-{nid}"
            else:
                sl = len(self.sections) + 1
                id = f"{self.anchor}-section-{sl}"
            if anchor is None:
                anchor = id

        assert anchor is not None
        assert id is not None

        # Prepend custom module anchor to the section if set
        cust_anchor = self.mod_cust_config.get("anchor")
        if cust_anchor:
            anchor = f"{cust_anchor}_{anchor}"
            id = f"{cust_anchor}_{id}"

        # Sanitise anchor ID and check for global duplicates
        anchor = report.save_htmlid(anchor)

        # Skip if user has a config to remove this module section
        if str(anchor) in config.remove_sections:
            logger.debug(f"Skipping section with anchor '{anchor}' because specified in user config")
            return

        # Skip if user has a config to remove this module section
        if str(id) in config.remove_sections:
            logger.debug(f"Skipping section with id '{id}' because specified in user config")
            return

        # See if we have a user comment in the config, but only if the section ID is different from the module ID
        # (otherwise it's a duplicate comment)
        if self.anchor != id and self.anchor != anchor:
            if str(id) in config.section_comments:
                comment = config.section_comments[str(id)]
            elif str(anchor) in config.section_comments:
                comment = config.section_comments[str(anchor)]

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
            name=name or "",
            anchor=Anchor(anchor),
            id=SectionId(id),
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
                section.plot_anchor = plot.anchor
                # separately keeping track of Plot objects to be rendered further
                report.plot_by_id[plot.anchor] = plot
            else:  # str
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

    def groups_for_sample(self, s_name: SampleName) -> Tuple[SampleGroup, Optional[str]]:
        """
        Takes a sample name and returns a trimmed name and groups it's assigned to.
        based on the patterns in config.sample_merge_groups.
        """
        if not config.table_sample_merge:
            return SampleGroup(s_name), None

        matched_label: Optional[str] = None
        grouping_exts: List[CleanPatternT]
        group_name = SampleGroup(s_name)
        for label, grouping_exts in config.table_sample_merge.items():
            if isinstance(grouping_exts, (str, dict)):
                grouping_exts = [grouping_exts]
            if grouping_exts:
                s_name_without_ext = SampleName(
                    self._clean_s_name(
                        s_name,
                        fn_clean_exts=grouping_exts,
                        fn_clean_trim=[],
                        prepend_dirs=False,
                    )
                )
                if s_name_without_ext != s_name:  # matched the label
                    matched_label = label
                    # Clean the rest of the name
                    group_name = SampleGroup(
                        # Leaving out fn_clean_exts and fn_clean_trim, so the default values are used, to make sure
                        # all default extentions are trimmed after we trimmed the groupping pattern.
                        self._clean_s_name(
                            s_name_without_ext,
                            prepend_dirs=False,
                        )
                    )
                    break

        return group_name, matched_label

    def group_samples_names(
        self, samples: Iterable[SampleName]
    ) -> Dict[SampleGroup, List[Tuple[Optional[str], SampleName, SampleName]]]:
        """
        Group sample name according to a named set of patterns defined in
        the config.sample_merge_groups dictionary.
        :param samples: sample names
        :return: a dict where the keys are group names, and the values are lists of tuples,
            of cleaned base names according to the cleaning rules and the original sample names
        """
        group_by_label: Dict[Optional[str], List[Tuple[SampleGroup, SampleName]]] = defaultdict(list)
        for original_name in sorted(samples):
            group_name, label = self.groups_for_sample(original_name)
            group_by_label[label].append((group_name, original_name))

        group_by_merged_name: Dict[SampleGroup, List[Tuple[Optional[str], SampleName]]] = defaultdict(list)
        for label, group in group_by_label.items():
            for group_name, original_name in group:
                group_by_merged_name[group_name].append((label, original_name))

        # Extend sample names in non-trivial groups with the group label
        return {
            group_name: [
                (
                    label,
                    SampleName(group_name) if (len(group) == 1 or not label) else SampleName(group_name + " " + label),
                    original_name,
                )
                for (label, original_name) in group
            ]
            for group_name, group in group_by_merged_name.items()
        }

    def group_samples_and_average_metrics(
        self,
        data_by_sample: Dict[Union[SampleName, str], Dict[Union[ColumnKey, str], ValueT]],
        grouping_config: SampleGroupingConfig,
    ) -> Dict[SampleGroup, List[InputRow]]:
        """
        Group samples and merges numeric metrics by averaging them, optionally normalizing using
        `normalization_metric_name`
        """

        rows_by_grouped_samples: Dict[SampleGroup, List[InputRow]] = defaultdict(list)
        for g_name, labels_s_names in self.group_samples_names([SampleName(s) for s in data_by_sample.keys()]).items():
            if len(labels_s_names) == 0:
                continue

            # We do not want "merged sample" clash with other real samples if the group is non-trivial,
            # so appending an ending to the "merged sample" name:
            if len(labels_s_names) > 1 and SampleName(g_name) in data_by_sample:
                g_name = SampleGroup(f"{g_name} (grouped)")

            # Just a single row for a trivial group
            if len(labels_s_names) == 1:
                _, s_name, original_s_name = labels_s_names[0]
                rows_by_grouped_samples[g_name] = [InputRow(sample=s_name, data=data_by_sample[original_s_name])]
                continue

            merged_row = InputRow(sample=SampleName(g_name), data={})

            # Init a dictionary of all cols that would be summed to serve as weights
            sum_by_col: Dict[ColumnKey, float] = dict()

            if grouping_config.cols_to_weighted_average:
                for _, weight_col_key in grouping_config.cols_to_weighted_average:
                    sum_by_col[weight_col_key] = 0

                # Calculate the weights
                for col in sum_by_col.keys():
                    for _, _, original_s_name in labels_s_names:
                        val = data_by_sample[original_s_name][col]
                        if isinstance(val, int) or isinstance(val, float):
                            sum_by_col[col] += float(val)

                for col, weight_col in grouping_config.cols_to_weighted_average:
                    weight = sum_by_col[weight_col]
                    if weight > 0:
                        merged_row.data[col] = (
                            sum(
                                [
                                    float(data_by_sample[original_s_name][col])
                                    * float(data_by_sample[original_s_name][weight_col])
                                    if (
                                        isinstance(data_by_sample[original_s_name][col], float)
                                        or isinstance(data_by_sample[original_s_name][col], int)
                                    )
                                    and (
                                        isinstance(data_by_sample[original_s_name][weight_col], float)
                                        or isinstance(data_by_sample[original_s_name][weight_col], int)
                                    )
                                    else 0
                                    for _, _, original_s_name in labels_s_names
                                ]
                            )
                            / weight
                        )

            if grouping_config.cols_to_average:
                for col in grouping_config.cols_to_average:
                    merged_row.data[col] = sum(
                        [
                            float(data_by_sample[original_s_name][col])
                            if (
                                isinstance(data_by_sample[original_s_name][col], float)
                                or isinstance(data_by_sample[original_s_name][col], int)
                            )
                            else 0
                            for _, _, original_s_name in labels_s_names
                        ]
                    ) / len(labels_s_names)

            if grouping_config.cols_to_sum:
                for col in grouping_config.cols_to_sum:
                    if col in sum_by_col:
                        merged_row.data[col] = sum_by_col[col]
                    else:
                        merged_row.data[col] = sum(
                            [
                                float(data_by_sample[original_s_name][col])
                                if (
                                    isinstance(data_by_sample[original_s_name][col], float)
                                    or isinstance(data_by_sample[original_s_name][col], int)
                                )
                                else 0
                                for _, _, original_s_name in labels_s_names
                            ]
                        )

            # Add count of fail statuses
            if grouping_config.extra_functions:
                for fn in grouping_config.extra_functions:
                    fn(merged_row, labels_s_names)

            rows_by_grouped_samples[g_name] = [merged_row] + [
                InputRow(sample=s_name, data=data_by_sample[original_s_name])
                for _, s_name, original_s_name in labels_s_names
            ]

        return rows_by_grouped_samples

    def clean_s_name(
        self,
        s_name: Union[str, List[str]],
        f: Union[LoadedFileDict[Any], FileDict],
        root: Optional[str] = None,
        filename: Optional[str] = None,
    ) -> str:
        """
        Helper function to take a long file name(s) and strip back to one clean sample name. Somewhat arbitrary.
        This is a user-facing version of _clean_s_name and the one that should be called in modules on raw sample names,
        because it gurantees the config options like config.prepend_dirs, config.fn_clean_exts, config.fn_clean_trim.

        search_pattern_key: the search pattern key that this file matched
        """
        return self._clean_s_name(
            s_name=s_name,
            f=f,
            root=root or f["root"],
            filename=filename or f["fn"],
            search_pattern_key=f["sp_key"],
        )

    def _clean_s_name(
        self,
        s_name: Union[str, List[str]],
        f: Optional[Union[LoadedFileDict[Any], FileDict]] = None,
        root: Optional[str] = None,
        filename: Optional[str] = None,
        search_pattern_key: Optional[str] = None,
        fn_clean_exts: Optional[List[Union[str, Dict[str, Union[str, List[str]]]]]] = None,
        fn_clean_trim: Optional[List[str]] = None,
        prepend_dirs: Optional[bool] = None,
    ) -> str:
        """
        Helper function to take a long file name(s) and strip back to one clean sample name. Somewhat arbitrary.

        search_pattern_key: the search pattern key that this file matched
        fn_clean_exts: patterns to use for cleaning (default: config.fn_clean_exts)
        fn_clean_trim: patterns to use for trimming (default: config.fn_clean_trim)
        prepend_dirs: boolean, whether to prepend dir name to s_name (default: config.prepend_dirs).
            requires `f` to be set.
        """
        if isinstance(s_name, list):
            if len(s_name) == 0:
                raise ValueError("Empty list of sample names passed to clean_s_name()")

            # Extract a sample name from a list of file names (for example, FASTQ pairs).
            # Each name is cleaned separately first:
            clean_names = [
                self._clean_s_name(
                    sn,
                    f=f,
                    root=root,
                    filename=filename,
                    search_pattern_key=search_pattern_key,
                    fn_clean_exts=fn_clean_exts,
                    fn_clean_trim=fn_clean_trim,
                    prepend_dirs=prepend_dirs,
                )
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

        sn = SampleNameMeta(original_name=SampleName(s_name))
        trimmed_name: SampleName = sn.original_name

        # Set string variables from f if it was a dict from find_log_files()
        if f is not None:
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
            trimmed_name = SampleName(filename)

        # if s_name comes from file contents, it may have a file path
        # For consistency with other modules, we keep just the basename
        trimmed_name = SampleName(os.path.basename(trimmed_name))

        if fn_clean_exts is None:
            fn_clean_exts = config.fn_clean_exts
        if fn_clean_trim is None:
            fn_clean_trim = config.fn_clean_trim
        if prepend_dirs is None:
            prepend_dirs = config.prepend_dirs

        # Prepend sample name with directory
        if prepend_dirs:
            sep = config.prepend_dirs_sep
            parts: Tuple[str, ...] = Path(root).parts if root else ()
            dirs: List[str] = [d.strip() for d in parts if d.strip() != ""]
            if config.prepend_dirs_depth != 0:
                d_idx = config.prepend_dirs_depth * -1
                if config.prepend_dirs_depth > 0:
                    dirs = dirs[d_idx:]
                else:
                    dirs = dirs[:d_idx]
            if len(dirs) > 0:
                trimmed_name = SampleName(f"{sep.join(dirs)}{sep}{trimmed_name}")

        if config.fn_clean_sample_names:
            # Split then take first section to remove everything after these matches
            _ext: Union[str, Dict[str, Union[str, List[str]]]]
            ext: Dict[str, Union[str, List[str]]]
            for _ext in fn_clean_exts:
                # Go through different filter types
                if isinstance(_ext, str):
                    ext = {"type": "truncate", "pattern": _ext}
                else:
                    ext = _ext

                # Check if this config is limited to a module
                if "module" in ext:
                    if isinstance(ext["module"], str):
                        ext["module"] = [ext["module"]]
                    if not any([m == self.anchor for m in ext["module"]]):
                        continue

                pattern = ext.get("pattern", "")
                assert isinstance(pattern, str)
                if ext.get("type") == "truncate":
                    trimmed_name = SampleName(str(trimmed_name).split(pattern, 1)[0])
                elif ext.get("type") in ("remove", "replace"):
                    if ext["type"] == "replace":
                        logger.warning(
                            "use 'config.fn_clean_sample_names.remove' instead "
                            "of 'config.fn_clean_sample_names.replace' [deprecated]"
                        )
                    trimmed_name = SampleName(str(trimmed_name).replace(pattern, ""))
                elif ext.get("type") == "regex":
                    trimmed_name = SampleName(re.sub(pattern, "", str(trimmed_name)))
                elif ext.get("type") == "regex_keep":
                    match = re.search(pattern, str(trimmed_name))
                    trimmed_name = SampleName(match.group()) if match else trimmed_name
                elif ext.get("type") is None:
                    logger.error(f'config.fn_clean_exts config was missing "type" key: {ext}')
                else:
                    logger.error(f"Unrecognised sample name cleaning pattern: {ext.get('type')}")
            # Trim off characters at the end of names
            for characters in fn_clean_trim:
                if trimmed_name.endswith(characters):
                    trimmed_name = SampleName(str(trimmed_name)[: -len(characters)])
                if trimmed_name.startswith(characters):
                    trimmed_name = SampleName(str(trimmed_name)[len(characters) :])

        # Remove trailing whitespace
        trimmed_name = SampleName(str(trimmed_name).strip())

        # If we cleaned back to an empty string, just use the original value
        if trimmed_name == "":
            trimmed_name = sn.original_name

        # Do any hard replacements that are set with --replace-names
        if config.sample_names_replace:
            for s_name_search, s_name_replace in config.sample_names_replace.items():
                try:
                    # Skip if we're looking for exact matches only
                    if config.sample_names_replace_exact:
                        # Simple strings
                        if not config.sample_names_replace_regex and str(trimmed_name) != s_name_search:
                            continue
                        # regexes
                        if config.sample_names_replace_regex and not re.fullmatch(s_name_search, trimmed_name):
                            continue
                    # Replace - regex
                    if config.sample_names_replace_regex:
                        trimmed_name = SampleName(re.sub(s_name_search, s_name_replace, str(trimmed_name)))
                    # Replace - simple string
                    else:
                        # Complete name swap
                        if config.sample_names_replace_complete:
                            if s_name_search in trimmed_name:
                                trimmed_name = SampleName(s_name_replace)
                        # Partial substring replace
                        else:
                            trimmed_name = SampleName(str(trimmed_name).replace(s_name_search, s_name_replace))
                except re.error as e:
                    logger.error(f"Error with sample name replacement regex: {e}")

        sn.trimmed_name = trimmed_name
        return trimmed_name

    def ignore_samples(
        self,
        data: Dict[SampleNameT, DataT],
        sample_names_ignore: Optional[List[str]] = None,
        sample_names_ignore_re: Optional[List[str]] = None,
    ) -> Dict[SampleNameT, DataT]:
        """Strip out samples which match `sample_names_ignore`"""
        try:
            if not isinstance(data, dict):  # type: ignore
                return data
            new_data: Dict[SampleNameT, DataT] = dict()
            for s_name, v in data.items():
                if not self.is_ignore_sample(s_name, sample_names_ignore, sample_names_ignore_re):
                    new_data[s_name] = v
            return new_data
        except (TypeError, AttributeError):
            return data

    @staticmethod
    def is_ignore_sample(
        s_name: Union[str, SampleName],
        sample_names_ignore: Optional[List[str]] = None,
        sample_names_ignore_re: Optional[List[str]] = None,
    ) -> bool:
        """Should a sample name be ignored?"""
        sample_names_ignore = sample_names_ignore or config.sample_names_ignore
        sample_names_ignore_re = sample_names_ignore_re or config.sample_names_ignore_re
        glob_match = any(fnmatch.fnmatch(s_name, sn) for sn in sample_names_ignore)
        re_match = any(re.match(sn, s_name) for sn in sample_names_ignore_re)
        return glob_match or re_match

    def general_stats_addcols(
        self,
        data_by_sample: Dict[Union[SampleName, str], Dict[Union[ColumnKey, str], ValueT]],
        headers: Optional[Dict[Union[ColumnKey, str], ColumnDict]] = None,
        namespace: Optional[str] = None,
        group_samples_config: SampleGroupingConfig = SampleGroupingConfig(),
    ):
        """Helper function to add to the General Statistics variable.
        Adds to report.general_stats and does not return anything. Fills
        in required config variables if not supplied.
        :param data_by_sample: A dict with the data. Key should be sample name, the data can be a key-value dict.
                     Or, for grouped samples, the key is the group name, and the data is a list of tuples with
                     the first element being the sample name in the group, and the second a key-value dict.
        :param headers: Dict with information for the headers,
                        such as colour scales, min and max values etc.
                        See docs/writing_python.md for more information.
        :param namespace: Append to the module name in the table column description.
                          Can be e.g. a submodule name.
        :param group_samples_config: Configuration for grouping samples.
        :return: None
        """
        if self.skip_generalstats:
            return

        rows_by_group: Dict[SampleGroup, List[InputRow]]
        if config.table_sample_merge:
            rows_by_group = self.group_samples_and_average_metrics(
                data_by_sample,
                group_samples_config,
            )
        else:
            rows_by_group = {
                SampleGroup(sname): [InputRow(sample=SampleName(sname), data=data)]
                for sname, data in data_by_sample.items()
            }

        _headers: Dict[ColumnKey, ColumnDict] = {}

        # Guess the column headers from the data if not supplied
        if headers is None or len(headers) == 0:
            column_ids: Set[ColumnKey] = set()
            for rows in rows_by_group.values():
                for row in rows:
                    column_ids.update(row.data.keys())
            for col_id in sorted(column_ids):
                _headers[col_id] = {}
        else:
            # Make a copy
            _headers = {ColumnKey(col_id): col_dict.copy() for col_id, col_dict in headers.items()}

        # Add the module name to the description if not already done
        for col_id in _headers.keys():
            # Prepend the namespace displayed in the table with the module name
            _col = _headers[col_id]
            namespace = _col["namespace"] if "namespace" in _col else namespace
            _headers[col_id]["namespace"] = self.name
            if namespace:
                _headers[col_id]["namespace"] = self.name + ": " + str(namespace)
            if "description" not in _headers[col_id]:
                _headers[col_id]["description"] = _col["title"] if "title" in _col else col_id

        # Append to report.general_stats for later assembly into table
        report.general_stats_data.append(rows_by_group)
        report.general_stats_headers.append(_headers)  # type: ignore

    def add_data_source(
        self,
        f: Optional[LoadedFileDict[Any]] = None,
        s_name: Optional[str] = None,
        path: Optional[Union[str, Path]] = None,
        module: Optional[str] = None,
        section: Optional[str] = None,
    ):
        if f is None and path is None:
            lint_error(f"add_data_source needs f or path to be set, got: {locals()}")
            return
        if module is None:
            module = self.name
        if section is None:
            section = "all_sections"
        if s_name is None and f is not None:
            s_name = f["s_name"]
        if s_name is None:
            return
        if self.is_ignore_sample(s_name):
            return
        if path is None and f is not None:
            path = os.path.abspath(os.path.join(f["root"], f["fn"]))
        report.data_sources[module][section][s_name] = str(path)

    def add_software_version(
        self,
        version: Optional[str] = None,
        sample: Optional[str] = None,
        software_name: Optional[str] = None,
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

    def write_data_file(self, data: Any, fn: str, sort_cols: bool = False, data_format: Optional[str] = None):
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
