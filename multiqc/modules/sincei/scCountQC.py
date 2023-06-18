""" MultiQC submodule to parse output from sincei scCountQC """

import logging
from collections import OrderedDict

from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)

class scCountQCMixin:
    def parse_scCountQC(self):
        """Find scCountQC output."""
        self.sincei_scCountQC = dict()
        for f in self.find_log_files("sincei/scCountQC"):
            parsed_data = self.parsescCountQCFile(f)
            for k, v in parsed_data.items():
                if k in self.sincei_scCountQC:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.sincei_scCountQC[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section="scCountQC")

        self.sincei_scCountQC = self.ignore_samples(self.sincei_scCountQC)

        if len(self.sincei_scCountQC) > 0:
            # Write data to file
            self.write_data_file(self.sincei_scCountQC, "sincei_count_qc")

            header = OrderedDict()
#            header["SampleName"] = {
#                "title": "Sample Name",
#                "description": "Name of Sample"
#                }
            header["n_genes"] = {
                "title": "# Features",
                "description": "No. of detected features (bins or genes) with non-zero counts",
                "scale": "RdBu",
                "min": 0,
            }
            header["n_counts_log"] = {
                "title": "# Counts (log1p)",
                "description": "Total counts in features per cell (logN+1 scale)",
                "scale": "OrRd",
                "min": 0,
            }
            header["pct_50"] = {
                "title": "% Counts top 50",
                "description": "Percent of alignments in top 50 features",
                "scale": "RdYlBu_r",
                "min": 0,
                "max": 100,
            }
            header["pct_100"] = {
                "title": "% Counts top 100",
                "description": "Percent of alignments in top 100 features",
                "scale": "RdYlBu_r",
                "min": 0,
                "max": 100,
            }
            header["pct_200"] = {
                "title": "% Counts top 200",
                "description": "Percent of alignments in top 200 features",
                "scale": "RdYlBu_r",
                "min": 0,
                "max": 100,
            }
            header["pct_500"] = {
                "title": "% Counts top 500",
                "description": "Percent of alignments in top 500 features",
                "scale": "RdYlBu_r",
                "min": 0,
                "max": 100,
            }
            header["gini_coefficient"] = {
                "title": "Gini Coefficient",
                "description": "Gini coefficient of enrichment (inequality) of counts in features.",
                "scale": "OrRd",
                "min": 0,
                "max": 1,
            }
            tdata = dict()
            for k, v in self.sincei_scCountQC.items():
                tdata[k] = {
                    "SampleName": v["sample"],
                    "n_genes": v["n_genes"],
                    "n_counts_log": v["total_log1p"],
                    "pct_50": v["pct_50"],
                    "pct_100": v["pct_100"],
                    "pct_200": v["pct_200"],
                    "pct_500": v["pct_500"],
                    "gini_coefficient": v["gini"],
                }
            config = {
                "namespace": "sincei scCountQC",
                "max_table_rows": 10000
                }
            self.add_section(
                name="Counting Metrics",
                anchor="scCountQC",
                description="Statistics of distribution of counts per cells after counting using `scCountQC`",
                plot=table.plot(tdata, header, config),

            )

        return len(self.sincei_scCountQC)

    def parsescCountQCFile(self, f):
        d = {}
        firstLine = True
        for line in f["f"].splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 12:
                # This is not really the output from scCountQC!
                log.warning(
                    "{} was initially flagged as the tabular output from scCountQC, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()

            s_name = self.clean_s_name(cols[0], f)
            if s_name in d:
                log.debug("Replacing duplicate sample {}.".format(s_name))
            d[s_name] = dict()

            try:
                d[s_name]["sample"] = cols[2]
                d[s_name]["n_genes"] = float(cols[3])
                d[s_name]["total_log1p"] = float(cols[6])
                d[s_name]["pct_50"] = float(cols[7])
                d[s_name]["pct_100"] = float(cols[8])
                d[s_name]["pct_200"] = float(cols[9])
                d[s_name]["pct_500"] = float(cols[10])
                d[s_name]["gini"] = float(cols[11])
            except:
                # Obviously this isn't really the output from scCountQC
                log.warning(
                    "{} was initially flagged as the output from scCountQC, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()
        return d
