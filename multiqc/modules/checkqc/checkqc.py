import json
import logging
import re
from operator import itemgetter

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

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
    """
    The module parses a CheckQC JSON file, so make sure to use CheckQC with the `--json` flag and collect the stdout in a file.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="CheckQC",
            anchor="checkqc",
            href="https://github.com/Molmed/checkQC",
            info="Checks a set of quality criteria against an Illumina runfolder.",
            comment="Samples are only shown in the report if they fail a check",
            doi="10.21105/joss.00556",
        )

        self.checkqc_data = dict()
        self.runs = set()

        self.log_files = list(self.find_log_files("checkqc"))
        self.onerun = len(self.log_files) == 1

        for f in self.log_files:
            raw_content = json.loads(f["f"])
            self.parse_checkqc_json(raw_content, f)
            self.add_data_source(f)

        if not self.checkqc_data:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.log_files)} run and {len(self.checkqc_data)} samples")

        self.write_data_file(self.checkqc_data, "multiqc_checkqc")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.add_sections()

    def parse_checkqc_json(self, content, f):
        run = self._get_unique_runname(content)
        if "ReadsPerSampleHandler" in content:
            self.get_reads_per_sample_data(content["ReadsPerSampleHandler"], run, f)
        if "ClusterPFHandler" in content:
            self.get_cluster_pf_data(content["ClusterPFHandler"], run, f)
        if "Q30Handler" in content:
            self.get_q30_data(content["Q30Handler"], run, f)
        if "ErrorRateHandler" in content:
            self.get_error_rate_data(content["ErrorRateHandler"], run, f)
        if "UndeterminedPercentageHandler" in content:
            self.get_undetermined_percentage_data(content["UndeterminedPercentageHandler"], run, f)
        if "UnidentifiedIndexHandler" in content:
            self.get_unidentified_index_data(content["UnidentifiedIndexHandler"], run, f)

    def _get_unique_runname(self, content):
        base_run_name = content["run_summary"]["instrument_and_reagent_type"]
        run_name = base_run_name
        i = 1
        while run_name in self.runs:
            run_name = f"{base_run_name}_{i}"
            i += 1
        self.runs.add(run_name)
        return run_name

    def add_sections(self):
        if "ReadsPerSampleHandler" in self.checkqc_data:
            self.add_reads_per_sample_section()
        if "ClusterPFHandler" in self.checkqc_data:
            self.add_cluster_pf_section()
        if "Q30Handler" in self.checkqc_data:
            self.add_q30_section()
        if "ErrorRateHandler" in self.checkqc_data:
            self.add_error_rate_section()
        if "UndeterminedPercentageHandler" in self.checkqc_data:
            self.add_undetermined_percentage_section()
        if "UndeterminedPercentageHandler_ZeroYield" in self.checkqc_data:
            self.add_undetermined_percentage_zero_yield_section()
        if "UnidentifiedIndexHandler" in self.checkqc_data:
            self.add_unidentified_index_section()

    def get_reads_per_sample_data(self, issues, run, f):
        """Parse data from checkQC ReadsPerSampleHandler

        Args:
            issues (dict): JSON dict from CheckQC containing ReadsPerSampleHandler results
            run (str): name of sequencing run
            f (dict): MultiQC log file dict
        """
        data = {}
        for issue in issues:
            is_error = issue["type"] == "error"
            sample_name = issue["data"]["sample_name"]
            lane = issue["data"]["lane"]
            if self.onerun:
                sample = self.clean_s_name(f"{sample_name} - Lane {lane}", f)
            else:
                sample = self.clean_s_name(f"{sample_name} - Lane {lane}, run {run}", f)
            if self.is_ignore_sample(sample):
                continue
            self.add_data_source(f, sample)

            read_num = issue["data"]["sample_reads"]
            read_threshold = issue["data"]["threshold"]

            data[sample] = {"read_num": read_num, "threshold": read_threshold}
            if is_error:
                data[sample]["missing_error"] = read_threshold - read_num
            else:
                data[sample]["missing_warning"] = read_threshold - read_num
        if data:
            if "ReadsPerSampleHandler" not in self.checkqc_data:
                self.checkqc_data["ReadsPerSampleHandler"] = data
            else:
                self.checkqc_data["ReadsPerSampleHandler"].update(data)

    def add_reads_per_sample_section(self):
        """Add a section for samples with too few reads

        Creates a barplot with read number and number missing to reach threshold

        """
        data = self.checkqc_data["ReadsPerSampleHandler"]

        warning, error = self._get_warning_error(data)

        cats = dict()
        cats["read_num"] = {
            "name": "Reads",
        }
        if warning:
            cats["missing_warning"] = {"name": "Reads missing to reach threshold for QC warning", "color": "#ffc300"}
        if error:
            cats["missing_error"] = {"name": "Reads missing to reach threshold for QC error", "color": "#f44336"}

        pconfig = {
            "id": "checkqc_reads-per-sample-plot",
            "title": "CheckQC: Number reads too low",
            "ylab": "Number of reads [M]",
            "xlab": "Sample - Lane",
        }

        self.add_section(
            name="Too few reads per sample",
            anchor="checkqc-readspersample",
            description="Checks that the number of reads assigned to a sample is high enough. See [documentation](https://checkqc.readthedocs.io/en/latest/checkQC.handlers.reads_per_sample_handler.html).",
            helptext="""The value specified in the configuration is interpreted as the number of reads demanded for a single sample,
            i.e. the number of reads per sample on a lane which has multiple samples is the threshold divided by the total number of samples on the lane.
            """,
            plot=bargraph.plot(data, cats, pconfig),
        )

    def get_cluster_pf_data(self, issues, run, f):
        """Get data from ClusterPFHandler

        Args:
            issues (dict): JSON dict from CheckQC containing ClusterPFHandler results
            run (str): name of sequencing run
            f (dict): MultiQC log file dict
        """
        data = {}
        for issue in issues:
            is_error = issue["type"] == "error"
            lane = str(issue["data"]["lane"])
            if self.onerun:
                sample = self.clean_s_name(lane, f)
            else:
                sample = self.clean_s_name(f"{lane} ({run})", f)

            lane_pf = issue["data"]["lane_pf"]
            threshold = issue["data"]["threshold"] * pow(10, 6)
            if self.is_ignore_sample(sample):
                continue
            self.add_data_source(f, sample)

            data[sample] = {"lane_pf": lane_pf}
            if is_error:
                data[sample]["missing_error"] = threshold - lane_pf
            else:
                data[sample]["missing_warning"] = threshold - lane_pf
        if data:
            if "ClusterPFHandler" not in self.checkqc_data:
                self.checkqc_data["ClusterPFHandler"] = data
            else:
                self.checkqc_data["ClusterPFHandler"].update(data)

    def add_cluster_pf_section(self):
        """Add a section for lanes with cluster pf too low

        Creates a barplot with cluster pf number and number missing to reach threshold

        """
        data = self.checkqc_data["ClusterPFHandler"]
        data = {f"Lane {k}": v for k, v in data.items()}

        warning, error = self._get_warning_error(data)

        cats = dict()
        cats["lane_pf"] = {
            "name": "Clusters passing filters",
        }

        if error:
            cats["missing_error"] = {"name": "Cluster PF missing to reach threshold for QC error", "color": "#f44336"}
        if warning:
            cats["missing_warning"] = {
                "name": "Cluster PF missing to reach threshold for QC warning",
                "color": "#ffc300",
            }

        pconfig = {
            "id": "checkqc_cluster-pf-plot",
            "title": "CheckQC: Too few clusters passing filters",
            "ylab": "Number of clusters",
        }

        self.add_section(
            name="Clusters passing filters",
            anchor="checkqc-clusterpf",
            description="Checks that the number of clusters passing filter (PF) on a sequencing lane passes the set criteria. See [docuementation](https://checkqc.readthedocs.io/en/latest/checkQC.handlers.cluster_pf_handler.html).",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def get_q30_data(self, issues, run, f):
        """Get data from Q30Handler

        Args:
            issues (dict): JSON dict from CheckQC containing Q30Handler results
            run (str): name of sequencing run
            f (dict): MultiQC log file dict
        """
        data = {}
        for issue in issues:
            is_error = issue["type"] == "error"
            lane = issue["data"]["lane"]
            read = issue["data"]["read"]
            if self.onerun:
                sample = self.clean_s_name(f"{lane} - R{read}", f)
            else:
                sample = self.clean_s_name(f"{lane} - R{read} ({run})", f)

            if self.is_ignore_sample(sample):
                continue
            self.add_data_source(f, sample)
            percent_q30 = issue["data"]["percent_q30"]
            threshold = issue["data"]["threshold"]
            data[sample] = {"percent_q30": percent_q30}
            if is_error:
                data[sample]["missing_error"] = threshold - percent_q30
            else:
                data[sample]["missing_warning"] = threshold - percent_q30
        if data:
            if "Q30Handler" not in self.checkqc_data:
                self.checkqc_data["Q30Handler"] = data
            else:
                self.checkqc_data["Q30Handler"].update(data)

    def add_q30_section(self):
        """Add a section for lanes with Q30 too low

        Creates a barplot with %Q30 and percentage missing to reach threshold

        """
        data = self.checkqc_data["Q30Handler"]
        data = {f"Lane {k}": v for k, v in data.items()}

        warning, error = self._get_warning_error(data)

        cats = dict()
        cats["percent_q30"] = {
            "name": "%Q30",
        }
        if warning:
            cats["missing_warning"] = {"name": "%Q30 missing to reach threshold for QC warning", "color": "#ffc300"}
        if error:
            cats["missing_error"] = {"name": "%Q30 missing to reach threshold for QC error", "color": "#f44336"}

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
            description="Checks that the percentage of bases in read 1 or read 2 with base quality over 30 (%Q30) passes the set criteria. See [documentation](https://checkqc.readthedocs.io/en/latest/checkQC.handlers.q30_handler.html).",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def get_error_rate_data(self, issues, run, f):
        """Get data from ErrorRateHandler

        Args:
            issues (dict): JSON dict from CheckQC containing ErrorRateHandler results
            run (str): name of sequencing run
            f (dict): MultiQC log file dict
        """
        data = {}
        for issue in issues:
            is_error = issue["type"] == "error"
            lane = issue["data"]["lane"]
            read = issue["data"]["read"]
            if self.onerun:
                sample = self.clean_s_name(f"{lane} - R{read}", f)
            else:
                sample = self.clean_s_name(f"{lane} - R{read} ({run})", f)

            if self.is_ignore_sample(sample):
                continue
            self.add_data_source(f, sample)
            error_rate = issue["data"]["error_rate"]
            threshold = issue["data"]["threshold"]
            data[sample] = {"threshold": threshold}
            if is_error:
                data[sample]["missing_error"] = error_rate - threshold
            else:
                data[sample]["missing_warning"] = error_rate - threshold
        if data:
            if "ErrorRateHandler" not in self.checkqc_data:
                self.checkqc_data["ErrorRateHandler"] = data
            else:
                self.checkqc_data["ErrorRateHandler"].update(data)

    def add_error_rate_section(self):
        """Add a section for lanes with error rate too high

        Creates a barplot with error rate and missing to reach threshold

        """
        data = self.checkqc_data["ErrorRateHandler"]
        data = {f"Lane {k}": v for k, v in data.items()}

        warning, error = self._get_warning_error(data)

        cats = dict()
        cats["threshold"] = {
            "name": "Error rate part until threshold",
        }
        if warning:
            cats["missing_warning"] = {"name": "Error rate part above threshold for QC warning", "color": "#ffc300"}
        if error:
            cats["missing_error"] = {"name": "Error rate part above threshold for QC error", "color": "#f44336"}

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
            description="Checks that the sequencing error rate per lane and read are below the specified threshold. See [documentation](https://checkqc.readthedocs.io/en/latest/checkQC.handlers.error_rate_handler.html).",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def get_undetermined_percentage_data(self, issues, run, f):
        """Get data from UndeterminedPercentageHandler

        Args:
            issues (dict): JSON dict from CheckQC containing UndeterminedPercentageHandler results
            run (str): name of sequencing run
            f (dict): MultiQC log file dict
        """
        zero_yield_data = {}
        data = {}
        for issue in issues:
            is_error = issue["type"] == "error"
            lane = str(issue["data"]["lane"])
            if self.onerun:
                sample = self.clean_s_name(lane, f)
            else:
                sample = self.clean_s_name(f"{lane} ({run})", f)

            if self.is_ignore_sample(sample):
                continue
            self.add_data_source(f, sample)
            p_undetermined = issue["data"]["percentage_undetermined"]
            if p_undetermined == "N/A":
                zero_yield_data[run] = {"lane": lane}
            else:
                threshold = issue["data"]["threshold"]
                computed_threshold = issue["data"]["computed_threshold"]
                phix = issue["data"]["phix_on_lane"]
                data[sample] = {"phix": phix, "threshold": threshold}
                if is_error:
                    data[sample]["missing_error"] = p_undetermined - computed_threshold
                else:
                    data[sample]["missing_warning"] = p_undetermined - computed_threshold
        if data:
            if "UndeterminedPercentageHandler" not in self.checkqc_data:
                self.checkqc_data["UndeterminedPercentageHandler"] = data
            else:
                self.checkqc_data["UndeterminedPercentageHandler"].update(data)
        if zero_yield_data:
            if "UndeterminedPercentageHandler_ZeroYield" not in self.checkqc_data:
                self.checkqc_data["UndeterminedPercentageHandler_ZeroYield"] = zero_yield_data
            else:
                self.checkqc_data["UndeterminedPercentageHandler_ZeroYield"].update(zero_yield_data)

    def add_undetermined_percentage_section(self):
        """Add a section for lanes with undetermined index percentage too high

        Creates a barplot with undetermined index percentage and missing to reach threshold

        Args:
            run (str): name of sequencing run
        """
        data = self.checkqc_data["UndeterminedPercentageHandler"]
        data = {f"Lane {k}": v for k, v in data.items()}

        warning, error = self._get_warning_error(data)
        cats = dict()
        cats["phix"] = {"name": r"% PhiX", "color": "#88a680"}
        cats["threshold"] = {
            "name": r"% undetermined indexes until threshold",
        }
        if warning:
            cats["missing_warning"] = {
                "name": r"% undetermined indexes above threshold for QC warning",
                "color": "#ffc300",
            }
        if error:
            cats["missing_error"] = {"name": r"% undetermined indexes above threshold for QC error", "color": "#f44336"}

        pconfig = {
            "id": "checkqc_undetermined-percentage-plot",
            "title": "CheckQC: Percentage undetermined indexes too high",
            "ylab": r"% undetermined indexes",
            "cpswitch": False,
        }

        self.add_section(
            name="% Undetermined indexes",
            anchor="checkqc-undetermined-rate",
            description="This checks that the percentage of undetermined reads on a lane is below the specified threshold. If there are no indexes specified for the lane, this will be skipped. See [documentation](https://checkqc.readthedocs.io/en/latest/checkQC.handlers.undetermined_percentage_handler.html).",
            plot=bargraph.plot(data, cats, pconfig),
        )

    def add_undetermined_percentage_zero_yield_section(self):
        """Add a section for lanes with zero yield

        Creates a table with the lane names with zero yield

        Args:
            run (str): name of sequencing run
        """
        data = self.checkqc_data["UndeterminedPercentageHandler_ZeroYield"]

        pconfig = {
            "id": "checkqc_zero-yield-table",
            "title": "CheckQC: Lanes with Yield 0",
            "scale": "Reds",
            "col1_header": "Run",
        }
        headers = {
            "lane": {
                "title": "Lane",
                "description": "Sequencing lane",
                "format": "{:,.0f}",
            },
        }

        self.add_section(
            name="Lanes with zero yield",
            anchor="checkqc-zero-yield",
            description="Some lanes have zero yield, no undetermined percentage computation possible.",
            plot=table.plot(data, headers, pconfig),
        )

    def get_unidentified_index_data(self, issues, run, f):
        """Get data for UnidentifiedIndexHandler

        Args:
            issues (dict): JSON dict from CheckQC containing UnidentifiedIndexHandler results
            run (str): name of sequencing run
            f (dict): MultiQC log file dict
        """
        idx_to_lane_to_rep = {}

        for issue in issues:
            msg = issue["data"]["msg"]
            m = re.match(
                r"Index: ([ATGC+]+) on lane: (\d+) was significantly "
                r"overrepresented \(([0-9.]+)%\) at significance "
                r"threshold of: [0-9.]+%\.",
                msg,
            )
            if not m:
                # TODO: Could add handling of other message types
                continue
            index, lane, overrep = m.groups()
            overrep = float(overrep)

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
        for idx, _ in sorted_idx:
            sample = self.clean_s_name(idx, f)

            if self.is_ignore_sample(sample):
                continue
            self.add_data_source(f, sample)
            if sample not in data:
                data[sample] = {}
            for lane, val in idx_to_lane_to_rep[idx].items():
                if self.onerun:
                    data[sample][lane] = val
                else:
                    data[sample][f"{lane} ({run})"] = val

        if data:
            if "UnidentifiedIndexHandler" not in self.checkqc_data:
                self.checkqc_data["UnidentifiedIndexHandler"] = data
            else:
                self.checkqc_data["UnidentifiedIndexHandler"].update(data)

    def add_unidentified_index_section(self):
        """Add a section for lanes with unidentified index percentage too high

        Creates a barplot with undetermined index overrepresentation by lane

        Args:
            run (str): name of sequencing run
        """
        data = self.checkqc_data["UnidentifiedIndexHandler"]

        lanes = set()
        for sample in data:
            for lane in data[sample]:
                lanes.add(lane)

        cats = dict()
        for lane in sorted(lanes):
            cats[lane] = {"name": f"Lane {lane}"}

        pconfig = {
            "id": "checkqc_unidentified-index-plot",
            "title": "CheckQC: Overrepresented unidentified indexes",
            "ylab": "% representation",
            "xlab": "Overrepresented indexes",
            "stacking": None,
        }

        self.add_section(
            name="Unidentified indexes",
            anchor="checkqc-unidentifiedpercentage",
            description="This check attempts to identify any indexes that are highly represented within unidentified reads. See the [documentation](https://checkqc.readthedocs.io/en/latest/checkQC.handlers.unidentified_index_handler.html).",
            helptext="This check ignores any indexes which have `N`’s in them. These are assumed to be read errors.",
            plot=bargraph.plot(data, cats, pconfig),
        )

    @staticmethod
    def _get_warning_error(data):
        warning = False
        error = False
        for sample in data:
            if "missing_error" in data[sample]:
                error = True
            if "missing_warning" in data[sample]:
                warning = True
            if warning and error:
                break
        return warning, error
