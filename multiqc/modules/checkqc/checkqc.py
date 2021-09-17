"""MultiQC module to parse CheckQC JSON output"""

import logging
import json
from collections import OrderedDict
from operator import itemgetter
import re
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

handlers = (
    "ClusterPFHandler",
    "Q30Handler",
    "ErrorRateHandler",
    "ReadsPerSampleHandler",
    "UnidentifiedIndexHandler",
    "UndeterminedPercentageHandler",
)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="CheckQC",
            anchor="mymod",
            href="https://github.com/Molmed/checkQC",
            info="CheckQC is a program designed to check a set of quality criteria against an Illumina runfolder.",
        )

        self.checkqc_content = dict()

        for f in self.find_log_files("checkqc"):
            self.checkqc_content[f["f"]] = json.loads(f["f"])
            self.parse_json(f)
            self.add_data_source(f)

        if not self.checkqc_content:
            raise UserWarning

    def parse_json(self, f):
        self.add_general_stats(self.checkqc_content[f["f"]], f)
        self.add_sections(self.checkqc_content[f["f"]], f)

    def add_general_stats(self, content, f):
        """Add data to the general statistics table

        Adds only results from ReadsPerSampleHandler, as this is the
        only handler that is sample-centric, the others are lane-centric

        Args:
            content (dict): JSON dict of checkqc JSON file
            f (dict): MultiQC log file dict
        """
        general_stats = dict()

        for issue in content.get("ReadsPerSampleHandler", []):
            sample = self.clean_s_name(issue["data"]["sample_name"], f)
            if self.is_ignore_sample(sample):
                continue
            read_num = issue["data"]["sample_reads"]
            read_threshold = issue["data"]["threshold"]
            general_stats[sample] = {"read_num": read_num, "read_threshold": read_threshold}

        if general_stats:
            headers = OrderedDict()
            headers["read_num"] = {
                "title": "Too few reads",
                "description": "Too few demultiplexed reads for sample. Compare with column 'Minimal reads threshold'.",
                "suffix": "M",
            }
            headers["read_threshold"] = {
                "title": "Minimal reads threshold",
                "description": "Threshold for minimal expected number of reads for sample.",
                "suffix": "M",
            }
            self.general_stats_addcols(general_stats, headers)

    def add_reads_per_sample_section(self, issues, f):
        """Add a section for samples with too few reads

        Creates a barplot with read number and number missing to reach threshold

        Args:
            issues (dict): JSON dict from CheckQC containing ReadsPerSampleHandler results
            f (dict): MultiQC log file dict
        """
        data = {}
        error = False
        warning = False
        for issue in issues:
            is_error = issue["type"] == "error"
            error = is_error or error
            warning = (not is_error) or warning
            sample_name = issue["data"]["sample_name"]
            lane = issue["data"]["lane"]
            sample = self.clean_s_name(f"{sample_name} - {lane}", f)
            if self.is_ignore_sample(sample):
                continue

            read_num = issue["data"]["sample_reads"] * pow(10, 6)
            read_threshold = issue["data"]["threshold"] * pow(10, 6)

            data[sample] = {"read_num": read_num}
            if is_error:
                data[sample]["missing_error"] = read_threshold - read_num
            else:
                data[sample]["missing_warning"] = read_threshold - read_num

        cats = OrderedDict()
        cats["read_num"] = {
            "name": "Reads",
        }

        if warning:
            cats["missing_warning"] = {"name": "Reads missing to reach threshold for QC warning", "color": "#ffc300"}
        if error:
            cats["missing_error"] = {"name": "Reads missing to reach threshold for QC error", "color": "#ff0000"}

        pconfig = {
            "id": "checkqc_reads-per-sample-plot",
            "title": "CheckQC: Number reads too low",
            "ylab": "Number of reads",
            "xlab": "Sample - Lane",
        }

        self.add_section(
            name="Too few reads per sample",
            anchor="checkqc-readspersample",
            description="Some samples have too few reads",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def add_cluster_pf_section(self, issues, f):
        """Add a section for lanes with cluster pf too low

        Creates a barplot with cluster pf number and number missing to reach threshold

        Args:
            issues (dict): JSON dict from CheckQC containing ClusterPFHandler results
            f (dict): MultiQC log file dict
        """
        data = {}
        error = False
        warning = False
        for issue in issues:
            is_error = issue["type"] == "error"
            error = is_error or error
            warning = (not is_error) or warning
            sample = self.clean_s_name(str(issue["data"]["lane"]), f)
            lane_pf = issue["data"]["lane_pf"]
            threshold = issue["data"]["threshold"] * pow(10, 6)
            if self.is_ignore_sample(sample):
                continue

            data[sample] = {"lane_pf": lane_pf}
            if is_error:
                data[sample]["missing_error"] = threshold - lane_pf
            else:
                data[sample]["missing_warning"] = threshold - lane_pf

        cats = OrderedDict()
        cats["lane_pf"] = {
            "name": "Cluster PF",
        }

        if error:
            cats["missing_error"] = {"name": "Cluster PF missing to reach threshold for QC error", "color": "#ff0000"}
        if warning:
            cats["missing_warning"] = {
                "name": "Cluster PF missing to reach threshold for QC warning",
                "color": "#ffc300",
            }

        pconfig = {
            "id": "checkqc_cluster-pf-plot",
            "title": "CheckQC: Cluster PF too low",
            "ylab": "Number of clusters",
            "xlab": "Lanes",
        }

        self.add_section(
            name="Cluster PF too low",
            anchor="checkqc-clusterpf",
            description="Some sequencing lanes have too few clusters passing filter (PF)",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def add_q30_section(self, issues, f):
        """Add a section for lanes with Q30 too low

        Creates a barplot with %Q30 and percentage missing to reach threshold

        Args:
            issues (dict): JSON dict from CheckQC containing Q30Handler results
            f (dict): MultiQC log file dict
        """
        data = {}
        error = False
        warning = False
        for issue in issues:
            is_error = issue["type"] == "error"
            error = is_error or error
            warning = (not is_error) or warning
            lane = issue["data"]["lane"]
            read = issue["data"]["read"]
            sample = self.clean_s_name(f"{lane} - {read}", f)
            if self.is_ignore_sample(sample):
                continue
            percent_q30 = issue["data"]["percent_q30"]
            threshold = issue["data"]["threshold"]
            data[sample] = {"percent_q30": percent_q30}
            if is_error:
                data[sample]["missing_error"] = threshold - percent_q30
            else:
                data[sample]["missing_warning"] = threshold - percent_q30

        cats = OrderedDict()
        cats["percent_q30"] = {
            "name": "%Q30",
        }
        if warning:
            cats["missing_warning"] = {"name": "%Q30 missing to reach threshold for QC warning", "color": "#ffc300"}
        if error:
            cats["missing_error"] = {"name": "%Q30 missing to reach threshold for QC error", "color": "#ff0000"}

        pconfig = {
            "id": "checkqc_q30-plot",
            "title": "CheckQC: %Q30 too low",
            "ylab": "%Q30",
            "xlab": "Lane - Read",
            "cpswitch": False,
        }

        self.add_section(
            name="%Q30 too low",
            anchor="checkqc-q30",
            description="Some lanes have too low %Q30 CheckQC. %Q30 is the percentage of bases in read 1 or read 2 with base quality over 30",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def add_error_rate_section(self, issues, f):
        """Add a section for lanes with error rate too high

        Creates a barplot with error rate and missing to reach threshold

        Args:
            issues (dict): JSON dict from CheckQC containing ErrorRateHandler results
            f (dict): MultiQC log file dict
        """
        data = {}
        error = False
        warning = False
        for issue in issues:
            is_error = issue["type"] == "error"
            error = is_error or error
            warning = (not is_error) or warning
            lane = issue["data"]["lane"]
            read = issue["data"]["read"]
            sample = self.clean_s_name(f"{lane} - {read}", f)
            if self.is_ignore_sample(sample):
                continue
            error_rate = issue["data"]["error_rate"]
            threshold = issue["data"]["threshold"]
            data[sample] = {"threshold": threshold}
            if is_error:
                data[sample]["missing_error"] = error_rate - threshold
            else:
                data[sample]["missing_warning"] = error_rate - threshold

        cats = OrderedDict()
        cats["threshold"] = {
            "name": "Error rate part until threshold",
        }
        if warning:
            cats["missing_warning"] = {"name": "Error rate part above threshold for QC warning", "color": "#ffc300"}
        if error:
            cats["missing_error"] = {"name": "Error rate part above threshold for QC error", "color": "#ff0000"}

        pconfig = {
            "id": "checkqc_error-rate-plot",
            "title": "CheckQC: Error rate too high",
            "ylab": "Error rate",
            "xlab": "Lane - Read",
            "cpswitch": False,
        }

        self.add_section(
            name="Error rate too high",
            anchor="checkqc-errorrate",
            description="Some lanes have too high error rate.",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def add_undetermined_percentage_section(self, issues, f):
        """Add a section for lanes with undetermined index percentage too high

        Creates a barplot with undetermined index percentage and missing to reach threshold

        Args:
            issues (dict): JSON dict from CheckQC containing UndeterminedPercentageHandler results
            f (dict): MultiQC log file dict
        """
        data = {}
        error = False
        warning = False
        for issue in issues:
            is_error = issue["type"] == "error"
            error = is_error or error
            warning = (not is_error) or warning
            sample = self.clean_s_name(str(issue["data"]["lane"]), f)
            if self.is_ignore_sample(sample):
                continue
            p_undetermined = issue["data"]["percentage_undetermined"]
            threshold = issue["data"]["threshold"]
            computed_threshold = issue["data"]["computed_threshold"]
            phix = issue["data"]["phix_on_lane"]
            data[sample] = {"phix": phix, "threshold": threshold}
            if is_error:
                data[sample]["missing_error"] = p_undetermined - computed_threshold
            else:
                data[sample]["missing_warning"] = p_undetermined - computed_threshold

        cats = OrderedDict()
        cats["phix"] = {"name": r"% PhiX", "color": "#000000"}
        cats["threshold"] = {
            "name": r"% undetermined indexes until threshold",
        }
        if warning:
            cats["missing_warning"] = {
                "name": r"% undetermined indexes above threshold for QC warning",
                "color": "#ffc300",
            }
        if error:
            cats["missing_error"] = {"name": r"% undetermined indexes above threshold for QC error", "color": "#ff0000"}

        pconfig = {
            "id": "checkqc_undetermined-percentage-plot",
            "title": "CheckQC: Percentage undetermined indexes too high",
            "ylab": r"% undetermined indexes",
            "xlab": "Lane",
            "cpswitch": False,
        }

        self.add_section(
            name="Percentage of undetermined indexes too high",
            anchor="checkqc-undeterminedrate",
            description="Some lanes have a percentage of undetermined indexes that is too high.",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def add_unidentified_index_section(self, issues, f):
        """Add a section for lanes with unidentified index percentage too high

        Creates a barplot with undetermined index overrepresentation by lane

        Args:
            issues (dict): JSON dict from CheckQC containing UndeterminedPercentageHandler results
            f (dict): MultiQC log file dict
        """
        error = False
        warning = False
        idx_to_lane_to_rep = {}
        threshold = None

        for issue in issues:
            is_error = issue["type"] == "error"
            error = is_error or error
            warning = (not is_error) or warning
            msg = issue["data"]["msg"]
            m = re.match(
                r"Index: ([ATGC+]+) on lane: (\d+) was significantly "
                r"overrepresented \(([0-9.]+)%\) at significance "
                r"threshold of: ([0-9.]+)%\.",
                msg,
            )
            if not m:
                # TODO: Could add handling of other message types
                continue
            index, lane, overrep, threshold = m.groups()
            overrep = float(overrep)
            threshold = float(threshold)

            if index not in idx_to_lane_to_rep:
                idx_to_lane_to_rep[index] = {}
            idx_to_lane_to_rep[index][lane] = overrep

        # Sort indexes by average overrepresentation and take top 20
        sorted_idx = sorted(
            [(idx, sum(l2r.values()) / len(l2r)) for idx, l2r in idx_to_lane_to_rep.items()],
            key=itemgetter(1),
            reverse=True,
        )[:20]

        data = {}
        lanes = set()
        for (idx, _) in sorted_idx:
            sample = self.clean_s_name(idx, f)
            data[sample] = {}
            for lane, val in idx_to_lane_to_rep[idx].items():
                data[sample][f"overrep_{lane}"] = val
                lanes.add(lane)

        cats = OrderedDict()
        for lane in sorted(lanes):
            cats[f"overrep_{lane}"] = {"name": f"Lane {lane}"}

        pconfig = {
            "id": "checkqc_unidentified-percentage-plot",
            "title": f"CheckQC: Overrepresented unidentified indexes (> {threshold}%)",
            "ylab": r"% representation",
            "xlab": f"Overrepresented indexes (> {threshold}%)",
            "stacking": None,
        }

        self.add_section(
            name="Overrepresented unidentified indexes",
            anchor="checkqc-unidentifiedpercentage",
            description="Some lanes have unidentified indexes that are overrepresented.",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def add_sections(self, content, f):
        if "ReadsPerSampleHandler" in content:
            self.add_reads_per_sample_section(content["ReadsPerSampleHandler"], f)
        if "ClusterPFHandler" in content:
            self.add_cluster_pf_section(content["ClusterPFHandler"], f)
        if "Q30Handler" in content:
            self.add_q30_section(content["Q30Handler"], f)
        if "ErrorRateHandler" in content:
            self.add_error_rate_section(content["ErrorRateHandler"], f)
        if "UndeterminedPercentageHandler" in content:
            self.add_undetermined_percentage_section(content["UndeterminedPercentageHandler"], f)
        if "UnidentifiedIndexHandler" in content:
            self.add_unidentified_index_section(content["UnidentifiedIndexHandler"], f)
