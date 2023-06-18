""" MultiQC submodule to parse output from sincei scFilterStats """

import logging
from collections import OrderedDict

from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


class scFilterStatsMixin:
    def parse_scFilterStats(self):
        """Find scFilterStats output. Only the output from --table is supported."""
        self.sincei_scFilterStats = dict()
        for f in self.find_log_files("sincei/scFilterStats"):
            parsed_data = self.parsescFilterStatsFile(f)
            for k, v in parsed_data.items():
                if k in self.sincei_scFilterStats:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.sincei_scFilterStats[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section="scFilterStats")

        self.sincei_scFilterStats = self.ignore_samples(self.sincei_scFilterStats)

        if len(self.sincei_scFilterStats) > 0:
            # Write data to file
            self.write_data_file(self.sincei_scFilterStats, "sincei_read_filtering")

            header = OrderedDict()
            header["M Entries"] = {"title": "M entries", "description": "Number of entries in the file (millions)"}
            header["pct_Aligned"] = {
                "title": "% Aligned",
                "description": "Percent of aligned entries",
                "scale": "YlGn",
                "min": 0,
                "max": 100,
            }
            header["pct_Filtered"] = {
                "title": "% Tot. Filtered",
                "description": "Percent of alignment that would be filtered for any reason.",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Blacklisted"] = {
                "title": "% Blacklisted",
                "description": "Percent of alignments falling (at least partially) inside a blacklisted region",
                "scale": "YlOrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_MAPQ"] = {
                "title": "% MAPQ",
                "description": "Percent of alignments having MAPQ scores below the specified threshold",
                "scale": "YlOrBn",
                "min": 0,
                "max": 100,
            }
            header["pct_Missing_Flags"] = {
                "title": "% Missing Flags",
                "description": "Percent of alignments lacking at least on flag specified by --samFlagInclude",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Forbidden_Flags"] = {
                "title": "% Forbidden Flags",
                "description": "Percent of alignments having at least one flag specified by --samFlagExclude",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_sincei_Dupes"] = {
                "title": "% sincei Duplicates",
                "description": "Percent of alignments marked by sincei as being duplicates",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Duplication"] = {
                "title": "% Duplication",
                "description": "Percent of alignments originally marked as being duplicates",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Singletons"] = {
                "title": "% Singletons",
                "description": "Percent of alignments that are singletons (i.e., paired-end reads where the mates don't align as a pair",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Excluded_Strand"] = {
                "title": "% Strand Filtered",
                "description": "Percent of alignments arising from the excluded strand",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Excluded_Motif"] = {
                "title": "% Motif Filtered",
                "description": "Percent of alignments lacking the expected sequence motif",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Excluded_GC"] = {
                "title": "% GC Filtered",
                "description": "Percent of alignments lacking the expected GC content",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Low_Aligned_Fraction"] = {
                "title": "% Low_Aligned_Fraction",
                "description": "Percent of alignments where the number of bases that match the reference were lower then desired",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            tdata = dict()
            for k, v in self.sincei_scFilterStats.items():
                tdata[k] = {
                    "N Entries": v["total"],
                    "pct_Filtered": 100.0 * v["filtered"] / float(v["total"]),
                    "pct_Blacklisted": 100.0 * v["blacklisted"] / float(v["total"]),
                    "pct_Below_MAPQ": 100.0 * v["mapq"] / float(v["total"]),
                    "pct_Missing_Flags": 100.0 * v["required_flags"] / float(v["total"]),
                    "pct_Forbidden_Flags": 100.0 * v["excluded_flags"] / float(v["total"]),
                    "pct_sincei_Dupes": 100.0 * v["internal_dupes"] / float(v["total"]),
                    "pct_Duplication": 100.0 * v["external_dupes"] / float(v["total"]),
                    "pct_Singletons": 100.0 * v["singletons"] / float(v["total"]),
                    "pct_Excluded_Strand": 100.0 * v["wrong_strand"] / float(v["total"]),
                    "pct_Excluded_Motif": 100.0 * v["wrong_motif"] / float(v["total"]),
                    "pct_Excluded_GC": 100.0 * v["unwanted_gc"] / float(v["total"]),
                    "pct_Low_Aligned_Fraction": 100.0 * v["low_alignedfrac"] / float(v["total"]),
                }

            config = {"namespace": "sincei bamPEFragmentSize"}
            self.add_section(
                name="Filtering metrics",
                anchor="scFilterStats",
                description="Estimated percentages of alignments filtered independently for each setting in `scFilterStats`",
                plot=table.plot(tdata, header, config),
            )

        return len(self.sincei_scFilterStats)

    def parsescFilterStatsFile(self, f):
        d = {}
        firstLine = True
        for line in f["f"].splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 14:
                # This is not really the output from scFilterStats!
                log.warning(
                    "{} was initially flagged as the tabular output from scFilterStats, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()

            s_name = self.clean_s_name(cols[0], f)
            if s_name in d:
                log.debug("Replacing duplicate sample {}.".format(s_name))
            d[s_name] = dict()

            try:
                d[s_name]["total"] = self._int(cols[1])
                d[s_name]["filtered"] = float(cols[2])
                d[s_name]["blacklisted"] = self._int(cols[3])
                d[s_name]["mapq"] = float(cols[4])
                d[s_name]["required_flags"] = float(cols[5])
                d[s_name]["excluded_flags"] = float(cols[6])
                d[s_name]["internal_dupes"] = float(cols[7])
                d[s_name]["external_dupes"] = float(cols[8])
                d[s_name]["singletons"] = float(cols[9])
                d[s_name]["wrong_strand"] = float(cols[10])
                d[s_name]["wrong_motif"] = float(cols[11])
                d[s_name]["unwanted_gc"] = float(cols[12])
                d[s_name]["low_alignedfrac"] = float(cols[13])
            except:
                # Obviously this isn't really the output from scFilterStats
                log.warning(
                    "{} was initially flagged as the output from scFilterStats, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()
        return d
