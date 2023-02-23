""" MultiQC submodule to parse output from MinIONQC summary stats """

import copy
import logging
import os
import re
from collections import OrderedDict

import yaml

from multiqc.modules.base_module import BaseMultiqcModule
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
            parsed_dict = self.parse_minionqc_report(s_name, f["f"])

            if parsed_dict is not None:
                if s_name in self.minionqc_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))

                # adds files used in MultiQC report
                self.add_data_source(f, s_name)

        # Filter to strip out ignored sample names
        self.minionqc_data = self.ignore_samples(self.minionqc_data)
        if len(self.minionqc_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.minionqc_data)))

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
            # Parsing as OrderedDict is slightly messier with YAML
            # http://stackoverflow.com/a/21048064/713980
            def dict_constructor(loader, node):
                return OrderedDict(loader.construct_pairs(node))

            yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)
            summary_dict = yaml.safe_load(f)
        except Exception as e:
            log.error("Error parsing MinIONQC input file: {}".format(f))
            return

        # Do a deep copy as dicts are immutable
        self.minionqc_raw_data[s_name] = copy.deepcopy(summary_dict)

        # get q value threshold used for reads
        q_threshold = None
        for k in summary_dict.keys():
            if k.startswith("Q>="):
                q_threshold = k

        data_dict = {}
        data_dict["all"] = summary_dict["All reads"]  # all reads
        data_dict["q_filt"] = summary_dict[q_threshold]  # quality filtered reads

        for q_key in ["all", "q_filt"]:
            for key_1 in ["reads", "gigabases"]:
                for key_2 in data_dict[q_key][key_1]:
                    new_key = "{} {}".format(key_1, key_2)
                    data_dict[q_key][new_key] = data_dict[q_key][key_1][key_2]
                data_dict[q_key].pop(key_1)  # removes key after flattening

        self.minionqc_data[s_name] = data_dict["all"]  # stats for all reads
        self.qfilt_data[s_name] = data_dict["q_filt"]  # stats for q-filtered reads
        self.q_threshold_list.add(q_threshold)  # quality threshold used in this file

    def headers_to_use(self):
        """
        Defines features of columns to be used in multiqc table
        """
        headers = OrderedDict()

        headers["total.reads"] = {
            "title": "Total reads",
            "description": "Total number of reads",
            "format": "{:,.0f}",
            "scale": "Greys",
        }
        headers["total.gigabases"] = {
            "title": "Total bases (GB)",
            "description": "Total bases",
            "format": "{:,.2f}",
            "scale": "Blues",
        }
        headers["N50.length"] = {
            "title": "Reads N50",
            "description": "Minimum read length needed to cover 50% of all reads",
            "format": "{:,.0f}",
            "scale": "Purples",
        }
        headers["mean.q"] = {
            "title": "Mean Q score",
            "description": "Mean quality of reads",
            "min": 0,
            "max": 15,
            "format": "{:,.1f}",
            "hidden": True,
            "scale": "Greens",
        }
        headers["median.q"] = {
            "title": "Median Q score",
            "description": "Median quality of reads",
            "min": 0,
            "max": 15,
            "format": "{:,.1f}",
            "scale": "Greens",
        }
        headers["mean.length"] = {
            "title": "Mean length (bp)",
            "description": "Mean read length",
            "format": "{:,.0f}",
            "hidden": True,
            "scale": "Blues",
        }
        headers["median.length"] = {
            "title": "Median length (bp)",
            "description": "Median read length",
            "format": "{:,.0f}",
            "scale": "Blues",
        }

        # Add row ID to avoid duplicates
        for k in headers:
            h_id = re.sub("[^0-9a-zA-Z]+", "_", headers[k]["title"])
            headers[k]["rid"] = "rid_{}".format(h_id)

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
                        {"name": "{}: Num reads".format(qfilt), "ylab": "# reads"},
                        {"name": "{}: Num gigabases".format(qfilt), "ylab": "# gigabases"},
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
