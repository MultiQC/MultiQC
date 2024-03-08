""" MultiQC submodule to parse output from MinIONQC summary stats """

import copy
import logging
import os
import re

import yaml

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="MinIONQC",
            anchor="minionqc",
            href="https://github.com/roblanf/minion_qc",
            info=" is a QC tool for Oxford Nanopore sequencing data",
            doi="10.1093/bioinformatics/bty654",
        )

        # Find and load any minionqc reports
        self.minionqc_raw_data = dict()  # main dataset in original YAML format
        self.minionqc_data = dict()  # main dataset. Stats from all reads
        self.qfilt_data = dict()  # Stats from quality filtered reads
        self.q_threshold_list = set()  # quality thresholds
        for f in self.find_log_files("minionqc", filehandles=True):
            # get sample name
            s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=os.path.dirname(f["root"]))

            # parses minionqc summary data
            self.parse_minionqc_report(s_name, f)

        # Filter to strip out ignored sample names
        self.minionqc_data = self.ignore_samples(self.minionqc_data)
        if len(self.minionqc_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.minionqc_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # columns to present in MultiQC summary table
        headers = self.headers_to_use()

        # Just the total number of reads in the main General Statistics table
        self.general_stats_addcols(self.minionqc_data, {"total.reads": headers["total.reads"]})

        # writes parsed data into a file (by MultiQC)
        self.write_data_file(self.minionqc_data, "multiqc_minionqc")

        # tables and plots in order
        self.table_qALL()
        self.table_qfiltered()
        self.plot_readlengths()

    def parse_minionqc_report(self, s_name, f):
        """
        Parses minionqc's 'summary.yaml' report file for results.
        Uses only the "All reads" stats. Ignores "Q>=x" part.
        """
        try:
            summary_dict = yaml.safe_load(f["f"])
        except Exception as e:
            log.error(f"Error parsing MinIONQC input file {f['f']}: {e}")
            return

        # adds files used in MultiQC report
        self.add_data_source(f, s_name)

        # Do a deep copy as dicts are immutable
        self.minionqc_raw_data[s_name] = copy.deepcopy(summary_dict)

        # get q value threshold used for reads
        q_threshold = None
        for k in summary_dict.keys():
            if k.startswith("Q>="):
                q_threshold = k

        data_dict = {"all": summary_dict["All reads"], "q_filt": summary_dict[q_threshold]}

        for q_key in ["all", "q_filt"]:
            for key_1 in ["reads", "gigabases"]:
                for key_2 in data_dict[q_key][key_1]:
                    new_key = f"{key_1} {key_2}"
                    data_dict[q_key][new_key] = data_dict[q_key][key_1][key_2]
                data_dict[q_key].pop(key_1)  # removes key after flattening

        if s_name in self.minionqc_data:
            log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

        self.minionqc_data[s_name] = data_dict["all"]  # stats for all reads
        self.qfilt_data[s_name] = data_dict["q_filt"]  # stats for q-filtered reads
        self.q_threshold_list.add(q_threshold)  # quality threshold used in this file

    @staticmethod
    def headers_to_use():
        """
        Defines features of columns to be used in multiqc table
        """
        headers = {
            "total.reads": {
                "title": "Total reads",
                "description": "Total number of reads",
                "format": "{:,.0f}",
                "scale": "Greys",
            },
            "total.gigabases": {
                "title": "Total bases (GB)",
                "description": "Total bases",
                "format": "{:,.2f}",
                "scale": "Blues",
            },
            "N50.length": {
                "title": "Reads N50",
                "description": "Minimum read length needed to cover 50% of all reads",
                "format": "{:,.0f}",
                "scale": "Purples",
            },
            "mean.q": {
                "title": "Mean Q score",
                "description": "Mean quality of reads",
                "min": 0,
                "max": 15,
                "format": "{:,.1f}",
                "hidden": True,
                "scale": "Greens",
            },
            "median.q": {
                "title": "Median Q score",
                "description": "Median quality of reads",
                "min": 0,
                "max": 15,
                "format": "{:,.1f}",
                "scale": "Greens",
            },
            "mean.length": {
                "title": "Mean length (bp)",
                "description": "Mean read length",
                "format": "{:,.0f}",
                "hidden": True,
                "scale": "Blues",
            },
            "median.length": {
                "title": "Median length (bp)",
                "description": "Median read length",
                "format": "{:,.0f}",
                "scale": "Blues",
            },
        }

        # Add row ID to avoid duplicates
        for k in headers:
            h_id = re.sub("[^0-9a-zA-Z]+", "_", headers[k]["title"])
            headers[k]["rid"] = f"rid_{h_id}"

        return headers

    def table_qALL(self):
        """Table showing stats for all reads"""

        self.add_section(
            name="Stats: All reads",
            anchor="minionqc-stats-qAll",
            description="MinIONQC statistics for all reads",
            plot=table.plot(
                self.minionqc_data,
                self.headers_to_use(),
                {
                    "namespace": "MinIONQC",
                    "id": "minionqc-stats-qAll-table",
                    "table_title": "MinIONQC Stats: All reads",
                },
            ),
        )

    def table_qfiltered(self):
        """Table showing stats for q-filtered reads"""

        description = "MinIONQC statistics for quality filtered reads. " + "Quailty threshold used: {}.".format(
            ", ".join(list(self.q_threshold_list))
        )
        if len(self.q_threshold_list) > 1:
            description += """
            <div class="alert alert-warning">
              <span class="glyphicon glyphicon-warning-sign"></span>
              <strong>Warning!</strong> More than one quality thresholds were present.
            </div>
            """
            log.warning(
                "More than one quality thresholds were present. Thresholds: {}.".format(
                    ", ".join(list(self.q_threshold_list))
                )
            )

        self.add_section(
            name="Stats: Quality filtered reads",
            anchor="minionqc-stats-qFilt",
            description=description,
            plot=table.plot(
                self.qfilt_data,
                self.headers_to_use(),
                {
                    "namespace": "MinIONQC",
                    "id": "minionqc-stats-qFilt-table",
                    "table_title": "MinIONQC Stats: Quality filtered reads",
                },
            ),
        )

    def plot_readlengths(self):
        pdata = [
            {s_name: d["All reads"]["reads"] for s_name, d in self.minionqc_raw_data.items()},
            {s_name: d["All reads"]["gigabases"] for s_name, d in self.minionqc_raw_data.items()},
        ]
        pconfig = {
            "id": "minionqc_read_lengths",
            "title": "MinIONQC: Output versus read length",
            "categories": True,
            "ylab": "# reads",
            "data_labels": [
                {"name": "All reads: Num reads", "ylab": "# reads"},
                {"name": "All reads: Num gigabases", "ylab": "# gigabases"},
            ],
        }
        for qfilt in list(self.q_threshold_list):
            try:
                pdata.extend(
                    [
                        {s_name: d[qfilt]["reads"] for s_name, d in self.minionqc_raw_data.items()},
                        {s_name: d[qfilt]["gigabases"] for s_name, d in self.minionqc_raw_data.items()},
                    ]
                )
                pconfig["data_labels"].extend(
                    [
                        {"name": f"{qfilt}: Num reads", "ylab": "# reads"},
                        {"name": f"{qfilt}: Num gigabases", "ylab": "# gigabases"},
                    ]
                )
            except KeyError:
                pass
        self.add_section(
            name="Read length output",
            anchor="minionqc-read-length-output",
            description="Number of reads / bp sequenced at given read length thresholds.",
            plot=linegraph.plot(pdata, pconfig=pconfig),
        )
