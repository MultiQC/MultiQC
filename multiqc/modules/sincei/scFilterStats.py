"""MultiQC submodule to parse output from sincei scFilterStats"""

import logging
import csv
import numpy as np
from itertools import groupby
from collections import OrderedDict

from multiqc.plots import beeswarm

# Initialise the logger
log = logging.getLogger(__name__)


class scFilterStatsMixin:
    def parse_scFilterStats(self):
        """Find scFilterStats output."""
        self.sincei_scFilterStats = dict()
        for f in self.find_log_files("sincei/scFilterStats", filehandles=True):
            parsed_data = self.parsescFilterStatsFile(f)
            for k, v in parsed_data.items():
                if k in self.sincei_scFilterStats:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.sincei_scFilterStats[k] = v
                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None)
            if len(parsed_data) > 0:
                self.add_data_source(f, section="scFilterStats")

        self.sincei_scFilterStats = self.ignore_samples(self.sincei_scFilterStats)

        if len(self.sincei_scFilterStats) > 0:
            # Write data to file
            self.write_data_file(self.sincei_scFilterStats, "sincei_read_filtering")

            header = OrderedDict()
            header["N Entries"] = {
                "title": "N entries",
                "description": "Median number of entries sampled from the file",
            }
            header["pct_Aligned"] = {
                "title": "% Aligned",
                "description": "Percent of aligned entries (Median of cells)",
                "scale": "YlGn",
                "min": 0,
                "max": 100,
            }
            header["pct_Filtered"] = {
                "title": "% Tot. Filtered",
                "suffix": "%",
                "description": "Percent of alignment that would be filtered for any reason (Median of cells)",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Blacklisted"] = {
                "title": "% Blacklisted",
                "suffix": "%",
                "description": "Percent of alignments falling (at least partially) inside a blacklisted region (Median of cells)",
                "scale": "YlOrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_MAPQ"] = {
                "title": "% MAPQ",
                "suffix": "%",
                "description": "Percent of alignments having MAPQ scores below the specified threshold (Median of cells)",
                "scale": "YlOrBn",
                "min": 0,
                "max": 100,
            }
            header["pct_Missing_Flags"] = {
                "title": "% Missing Flags",
                "suffix": "%",
                "description": "Percent of alignments lacking at least on flag specified by --samFlagInclude (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Forbidden_Flags"] = {
                "title": "% Forbidden Flags",
                "suffix": "%",
                "description": "Percent of alignments having at least one flag specified by --samFlagExclude (Median of cells)",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_sincei_Dupes"] = {
                "title": "% sincei Duplicates",
                "suffix": "%",
                "description": "Percent of alignments marked by sincei as being duplicates (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Duplication"] = {
                "title": "% Duplication",
                "suffix": "%",
                "description": "Percent of alignments originally marked as being duplicates (Median of cells)",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Singletons"] = {
                "title": "% Singletons",
                "suffix": "%",
                "description": "Percent of alignments that are singletons (i.e., paired-end reads where the mates don't align as a pair (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Excluded_Strand"] = {
                "title": "% Strand Filtered",
                "suffix": "%",
                "description": "Percent of alignments arising from the excluded strand (Median of cells)",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Excluded_Motif"] = {
                "title": "% Motif Filtered",
                "suffix": "%",
                "description": "Percent of alignments lacking the expected sequence motif (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Excluded_GC"] = {
                "title": "% GC Filtered",
                "suffix": "%",
                "description": "Percent of alignments lacking the expected GC content (Median of cells)",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Low_Aligned_Fraction"] = {
                "title": "% Low_Aligned_Fraction",
                "suffix": "%",
                "description": "Percent of alignments where the number of bases that match the reference were lower then desired (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            kys = list(list(self.sincei_scFilterStats.values())[0].keys())[1:]

            test_dict = self.getDictVal(self.sincei_scFilterStats, kys[0])
            out = {}
            for k in test_dict.keys():
                out[k] = dict.fromkeys(kys)
                for p in kys:
                    dv = self.getDictVal(self.sincei_scFilterStats, p)
                    out[k].update(dv[k])

            tdata = dict()
            for k, v in out.items():
                tdata[k] = {
                    "SampleName": k,
                    "N Entries": v["Total_sampled"],
                    "pct_Filtered": v["Filtered"],
                    "pct_Blacklisted": v["Blacklisted"],
                    "pct_Below_MAPQ": v["Low_MAPQ"],
                    "pct_Missing_Flags": v["Missing_Flags"],
                    "pct_Forbidden_Flags": v["Excluded_Flags"],
                    "pct_sincei_Dupes": v["Internal_Duplicates"],
                    "pct_Duplication": v["Marked_Duplicates"],
                    "pct_Singletons": v["Singletons"],
                    "pct_Excluded_Strand": v["Wrong_strand"],
                    "pct_Excluded_Motif": v["Wrong_motif"],
                    "pct_Excluded_GC": v["Unwanted_GC_content"],
                    "pct_Low_Aligned_Fraction": v["Low_aligned_fraction"],
                }
            config = {"namespace": "sincei scFilterStats", "max_table_rows": 10000}
            self.add_section(
                name="Filtering metrics",
                anchor="scFilterStats",
                description="Estimated percentages of alignments filtered independently for each setting in `scFilterStats`",
                plot=beeswarm.plot(tdata, header, config),
            )

        return len(self.sincei_scFilterStats)

    def parsescFilterStatsFile(self, f):
        reader = csv.DictReader(f["f"], delimiter="\t")
        if (
            len(set(["Total_sampled", "Filtered", "Blacklisted", "Wrong_motif"]).difference(set(reader.fieldnames)))
            != 0
        ):
            # This is not really the output from scFilterStats!
            log.warning(
                "{} was initially flagged as the tabular output from scFilterStats, but that seems to not be the case. Skipping...".format(
                    f["fn"]
                )
            )

        d = {}
        print(reader.fieldnames)
        for row in reader:
            s_name = self.clean_s_name(row["Cell_ID"], f)
            if s_name in d:
                log.debug("Replacing duplicate cell_id {}.".format(s_name))
            d[s_name] = dict()

            try:
                for key in reader.fieldnames:
                    d[s_name][key] = row[key]
            except:  # noqa: E722
                # Obviously this isn't really the output from scFilterStats
                log.warning(
                    "{} was initially flagged as the output from scFilterStats, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()
        return d

    def getDictVal(self, dat, val):
        dc = groupby(
            sorted(dat.items(), key=lambda x: x[1]["Cell_ID"].split("::")[0]),
            lambda x: x[1]["Cell_ID"].split("::")[0],
        )
        out = {i: {val: np.median([float(j[1][val]) for j in j])} for i, j in dc}
        return out
