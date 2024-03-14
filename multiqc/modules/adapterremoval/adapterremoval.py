""" MultiQC module to parse output from Adapter Removal """


import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"AdapterRemoval ver. ([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Adapter Removal",
            anchor="adapterremoval",
            target="Adapter Removal",
            href="https://github.com/MikkelSchubert/adapterremoval",
            info=" rapid adapter trimming, identification, and read merging ",
            doi=["10.1186/s13104-016-1900-2", "10.1186/1756-0500-5-337"],
        )

        self.__read_type = None
        self.__any_paired = False
        self.__collapsed = None
        self.__any_collapsed = False
        self.s_name = None
        self.adapter_removal_data = {}

        self.len_dist_plot_data = {
            "mate1": dict(),
            "mate2": dict(),
            "singleton": dict(),
            "collapsed": dict(),
            "collapsed_truncated": dict(),
            "discarded": dict(),
            "all": dict(),
        }

        parsed_data = None
        for f in self.find_log_files("adapterremoval", filehandles=True):
            self.s_name = f["s_name"]
            try:
                parsed_data = self.parse_settings_file(f)
            except ModuleNoSamplesFound:
                continue
            if parsed_data is not None:
                self.adapter_removal_data[self.s_name] = parsed_data
                self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.adapter_removal_data = self.ignore_samples(self.adapter_removal_data)

        if len(self.adapter_removal_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.adapter_removal_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.adapter_removal_data, "multiqc_adapter_removal")

        # add data to Basic Stats Table
        self.adapter_removal_stats_table()

        self.adapter_removal_retained_chart()
        self.adapter_removal_length_dist_plot()

    def parse_settings_file(self, f):
        self.result_data = {
            "total": None,
            "unaligned": None,
            "aligned": None,
            "reads_total": None,
            "retained": None,
            "percent_aligned": None,
        }

        settings_data = {"header": []}

        block_title = None
        for i, line in enumerate(f["f"]):
            line = line.rstrip("\n")
            if line == "":
                continue

            if line.startswith("AdapterRemoval"):
                version_match = re.search(VERSION_REGEX, line)
                if version_match:
                    self.add_software_version(version_match.group(1), self.s_name)

            if not block_title:
                block_title = "header"
                settings_data[block_title].append(str(line))
                continue

            if line.startswith("["):
                block_title = str(line.strip("[]"))
                settings_data[block_title] = []
                continue

            settings_data[block_title].append(str(line))

        # set data for further working
        self.set_result_data(settings_data)

        return self.result_data

    def set_result_data(self, settings_data):
        # set read and collapsed type
        self.set_ar_type(settings_data["Length distribution"])

        # set result_data
        self.set_trim_stat(settings_data["Trimming statistics"])
        self.set_len_dist(settings_data["Length distribution"])

    def set_ar_type(self, len_dist_data):
        head_line = len_dist_data[0].rstrip("\n").split("\t")

        self.__read_type = "paired" if head_line[2] == "Mate2" else "single"
        if not self.__any_paired:
            self.__any_paired = True if head_line[2] == "Mate2" else False

        self.__collapsed = True if head_line[-3] == "CollapsedTruncated" else False
        if not self.__any_collapsed:
            self.__any_collapsed = True if head_line[-3] == "CollapsedTruncated" else False

        # biological/technical relevance is not clear -> skip
        if self.__read_type == "single" and self.__collapsed:
            log.warning(f"Case single-end and collapse is not implemented -> File {self.s_name} skipped")
            raise ModuleNoSamplesFound

    def set_trim_stat(self, trim_data):
        required = [
            "total",
            "unaligned",
            "aligned",
            "discarded_m1",
            "singleton_m1",
            "retained",
            "discarded_m2",
            "singleton_m2",
            "full-length_cp",
            "truncated_cp",
        ]
        data_pattern = {"total": 0, "unaligned": 1, "aligned": 2, "discarded_m1": 3, "singleton_m1": 4, "retained": 6}

        if self.__read_type == "paired":
            data_pattern["discarded_m2"] = 5
            data_pattern["singleton_m2"] = 6
            if not self.__collapsed:
                data_pattern["retained"] = 8
            else:
                data_pattern["full-length_cp"] = 8
                data_pattern["truncated_cp"] = 9
                data_pattern["retained"] = 10

        for field in required:
            if field in data_pattern:
                tmp = trim_data[data_pattern[field]]
                value = tmp.split(": ")[1]
                self.result_data[field] = int(value)
            else:
                self.result_data[field] = 0

        reads_total = self.result_data["total"]
        aligned_total = self.result_data["aligned"]
        unaligned_total = self.result_data["unaligned"]
        if self.__read_type == "paired":
            reads_total = self.result_data["total"] * 2
            aligned_total = self.result_data["aligned"] * 2
            unaligned_total = self.result_data["unaligned"] * 2

        self.result_data["aligned_total"] = aligned_total
        self.result_data["unaligned_total"] = unaligned_total
        self.result_data["reads_total"] = reads_total
        self.result_data["discarded_total"] = reads_total - self.result_data["retained"]
        if self.__read_type == "paired":
            if not self.__collapsed:
                self.result_data["paired_reads"] = (
                    self.result_data["retained"] - self.result_data["singleton_m1"] - self.result_data["singleton_m2"]
                )
            else:
                self.result_data["paired_reads"] = (
                    self.result_data["retained"]
                    - self.result_data["full-length_cp"]
                    - self.result_data["truncated_cp"]
                    - self.result_data["singleton_m1"]
                    - self.result_data["singleton_m2"]
                )
            full_length_cp = self.result_data["full-length_cp"] * 2
            truncated_cp = self.result_data["truncated_cp"] * 2
            self.result_data["full-length_cp"] = full_length_cp
            self.result_data["truncated_cp"] = truncated_cp
        try:
            self.result_data["percent_aligned"] = (
                float(self.result_data["aligned"]) * 100.0 / float(self.result_data["total"])
            )
        except ZeroDivisionError:
            self.result_data["percent_aligned"] = 0
        if self.__collapsed:
            try:
                self.result_data["percent_collapsed"] = (
                    float(self.result_data["full-length_cp"] + self.result_data["truncated_cp"])
                    * 100.0
                    / float(self.result_data["reads_total"])
                )
            except ZeroDivisionError:
                self.result_data["percent_collapsed"] = 0
        try:
            self.result_data["percent_discarded"] = (
                float(self.result_data["discarded_m1"] + self.result_data["discarded_m2"])
                * 100.0
                / float(self.result_data["reads_total"])
            )
        except ZeroDivisionError:
            self.result_data["percent_discarded"] = 0

    def set_len_dist(self, len_dist_data):
        for line in len_dist_data[1:]:
            l_data = line.rstrip("\n").split("\t")
            l_data = list(map(int, l_data))

            # initialize file name
            if self.s_name not in self.len_dist_plot_data["mate1"]:
                self.len_dist_plot_data["mate1"][self.s_name] = dict()
                self.len_dist_plot_data["mate2"][self.s_name] = dict()
                self.len_dist_plot_data["singleton"][self.s_name] = dict()
                self.len_dist_plot_data["collapsed"][self.s_name] = dict()
                self.len_dist_plot_data["collapsed_truncated"][self.s_name] = dict()
                self.len_dist_plot_data["discarded"][self.s_name] = dict()
                self.len_dist_plot_data["all"][self.s_name] = dict()

            if self.__read_type == "single":
                if not self.__collapsed:
                    self.len_dist_plot_data["mate1"][self.s_name][l_data[0]] = l_data[1]
                    self.len_dist_plot_data["mate2"][self.s_name][l_data[0]] = None
                    self.len_dist_plot_data["singleton"][self.s_name][l_data[0]] = None
                    self.len_dist_plot_data["collapsed"][self.s_name][l_data[0]] = None
                    self.len_dist_plot_data["collapsed_truncated"][self.s_name][l_data[0]] = None
                    self.len_dist_plot_data["discarded"][self.s_name][l_data[0]] = l_data[2]
                    self.len_dist_plot_data["all"][self.s_name][l_data[0]] = l_data[3]
                else:
                    # this case should not be reached (see case at method set_ar_type())
                    pass
            else:
                if not self.__collapsed:
                    self.len_dist_plot_data["mate1"][self.s_name][l_data[0]] = l_data[1]
                    self.len_dist_plot_data["mate2"][self.s_name][l_data[0]] = l_data[2]
                    self.len_dist_plot_data["singleton"][self.s_name][l_data[0]] = l_data[3]
                    self.len_dist_plot_data["collapsed"][self.s_name][l_data[0]] = None
                    self.len_dist_plot_data["collapsed_truncated"][self.s_name][l_data[0]] = None
                    self.len_dist_plot_data["discarded"][self.s_name][l_data[0]] = l_data[4]
                    self.len_dist_plot_data["all"][self.s_name][l_data[0]] = l_data[5]
                else:
                    self.len_dist_plot_data["mate1"][self.s_name][l_data[0]] = l_data[1]
                    self.len_dist_plot_data["mate2"][self.s_name][l_data[0]] = l_data[2]
                    self.len_dist_plot_data["singleton"][self.s_name][l_data[0]] = l_data[3]
                    self.len_dist_plot_data["collapsed"][self.s_name][l_data[0]] = l_data[4]
                    self.len_dist_plot_data["collapsed_truncated"][self.s_name][l_data[0]] = l_data[5]
                    self.len_dist_plot_data["discarded"][self.s_name][l_data[0]] = l_data[6]
                    self.len_dist_plot_data["all"][self.s_name][l_data[0]] = l_data[7]

    def adapter_removal_stats_table(self):
        headers = {
            "percent_aligned": {
                "title": "% Trimmed",
                "description": "% trimmed reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn-rev",
                "shared_key": "percent_aligned",
            },
            "aligned_total": {
                "title": f"{config.read_count_prefix} Reads Trimmed",
                "description": f"Total trimmed reads ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "min": 0,
                "scale": "PuBu",
                "shared_key": "read_count",
            },
        }
        if self.__any_collapsed:
            headers["percent_collapsed"] = {
                "title": "% Collapsed",
                "description": "% collapsed reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn-rev",
                "shared_key": "percent_aligned",
            }
        headers["percent_discarded"] = {
            "title": "% Discarded",
            "description": "% discarded reads",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn-rev",
            "shared_key": "percent_discarded",
        }
        self.general_stats_addcols(self.adapter_removal_data, headers)

    def adapter_removal_retained_chart(self):
        pconfig = {
            "title": "Adapter Removal: Discarded Reads",
            "id": "ar_retained_plot",
            "ylab": "# Reads",
            "hide_zero_cats": False,
            "cpswitch_counts_label": "Number of Reads",
        }

        cats_pec = {}
        if self.__any_paired:
            cats_pec["paired_reads"] = {"name": "Uncollapsed Paired Reads"}

        cats_pec["singleton_m1"] = {"name": "Singleton R1"}

        if self.__any_paired:
            cats_pec["singleton_m2"] = {"name": "Singleton R2"}

            if self.__any_collapsed:
                cats_pec["full-length_cp"] = {"name": "Full-length Collapsed Reads"}
                cats_pec["truncated_cp"] = {"name": "Truncated Collapsed Reads"}

        cats_pec["discarded_m1"] = {"name": "Discarded R1"}

        if self.__any_paired:
            cats_pec["discarded_m2"] = {"name": "Discarded R2"}
        if self.__any_collapsed:
            retained_chart_description = "The number of input sequences that were retained, collapsed, and discarded. Be aware that the number of collapsed reads in the output FASTQ will be half of the numbers displayed in this plot, because both R1 and R2 of the collapsed sequences are counted here."
        else:
            retained_chart_description = "The number of input sequences that were retained and discarded."
        self.add_section(
            name="Retained and Discarded Reads",
            anchor="adapter_removal_retained_plot",
            description=retained_chart_description,
            plot=bargraph.plot(self.adapter_removal_data, cats_pec, pconfig),
        )

    def adapter_removal_length_dist_plot(self):
        pconfig = {
            "title": "Adapter Removal: Length Distribution",
            "id": "ar_length_count_plot",
            "ylab": "Counts",
            "xlab": "read length",
            "xDecimals": False,
            "ymin": 0,
            "tt_label": "<b>{point.x} bp trimmed</b>: {point.y:.0f}",
            "data_labels": None,
        }

        lineplot_data = [self.len_dist_plot_data["all"], self.len_dist_plot_data["mate1"]]
        data_labels = [
            {"name": "All", "ylab": "Count"},
            {"name": "Mate1", "ylab": "Count"},
        ]
        if self.__any_paired:
            lineplot_data.extend([self.len_dist_plot_data["mate2"], self.len_dist_plot_data["singleton"]])
            data_labels.extend(
                [
                    {"name": "Mate2", "ylab": "Count"},
                    {"name": "Singleton", "ylab": "Count"},
                ]
            )
            if self.__any_collapsed:
                lineplot_data.extend(
                    [self.len_dist_plot_data["collapsed"], self.len_dist_plot_data["collapsed_truncated"]]
                )
                data_labels.extend(
                    [{"name": "Collapsed", "ylab": "Count"}, {"name": "Collapsed Truncated", "ylab": "Count"}]
                )
        lineplot_data.append(self.len_dist_plot_data["discarded"])
        data_labels.append({"name": "Discarded", "ylab": "Count"})

        pconfig["data_labels"] = data_labels

        self.add_section(
            name="Length Distribution",
            anchor="ar_length_count",
            description="The length distribution of reads after processing adapter alignment.",
            plot=linegraph.plot(lineplot_data, pconfig),
        )
