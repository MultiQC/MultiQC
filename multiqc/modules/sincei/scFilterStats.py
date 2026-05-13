import logging
import csv

from multiqc.plots import table

from ._helpers import group_median_by_cell_prefix

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
                    log.warning(f"Replacing duplicate sample {k}.")
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

            header = dict()
            header["N Entries"] = {
                "title": "N entries",
                "description": "Median number of entries sampled from the file",
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
            header["pct_Below_MAPQ"] = {
                "title": "% MAPQ",
                "suffix": "%",
                "description": "Percent of alignments having MAPQ scores below the specified threshold (Median of cells)",
                "scale": "YlOrBr",
                "min": 0,
                "max": 100,
            }
            header["pct_Missing_Flags"] = {
                "title": "% Missing",
                "suffix": "%",
                "description": "Percent of alignments lacking at least on flag specified by --samFlagInclude (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Forbidden_Flags"] = {
                "title": "% Forbidden",
                "suffix": "%",
                "description": "Percent of alignments having at least one flag specified by --samFlagExclude (Median of cells)",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            }
            header["pct_sincei_Dupes"] = {
                "title": "% sincei Dups",
                "suffix": "%",
                "description": "Percent of alignments marked by sincei as being duplicates (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            header["pct_Duplication"] = {
                "title": "% Dups",
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
                "title": "% Low Aligned Fraction",
                "suffix": "%",
                "description": "Percent of alignments where the number of bases that match the reference were lower then desired (Median of cells)",
                "scale": "PuRd",
                "min": 0,
                "max": 100,
            }
            kys = [
                "Total_sampled",
                "Filtered",
                "Blacklisted",
                "Low_MAPQ",
                "Missing_Flags",
                "Excluded_Flags",
                "Internal_Duplicates",
                "Marked_Duplicates",
                "Singletons",
                "Wrong_strand",
                "Wrong_motif",
                "Unwanted_GC_content",
                "Low_aligned_fraction",
            ]

            test_dict = group_median_by_cell_prefix(self.sincei_scFilterStats, kys[0])
            out = {}
            for k in test_dict.keys():
                out[k] = dict.fromkeys(kys)
                for p in kys:
                    dv = group_median_by_cell_prefix(self.sincei_scFilterStats, p)
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
            config = {
                "id": "sincei-scFilterStats-table",
                "title": "sincei: scFilterStats filtering metrics",
                "namespace": "sincei scFilterStats",
            }
            self.add_section(
                name="Filtering metrics",
                anchor="scFilterStats",
                description=(
                    "Estimated percentages of alignments filtered for each criterion in "
                    "[`scFilterStats`](https://sincei.readthedocs.io/en/latest/content/tools/scFilterStats.html)."
                ),
                helptext="""
                `scFilterStats` samples reads from across the genome and estimates how many would
                be discarded by each filter that sincei applies during downstream analysis.

                **Filters are evaluated independently per alignment**, so the same read can be
                counted in several columns at once. As the sincei docs note, the sum of these
                percentages may exceed 100%, and they do not partition the reads into
                non-overlapping groups.

                The columns (medians across cells in each sample) are:

                - **% Tot. Filtered**: reads that would be removed for *any* reason
                - **% Blacklisted**: reads overlapping blacklisted regions
                - **% MAPQ**: alignments below the mapping-quality threshold
                - **% Missing Flags / % Forbidden Flags**: alignments failing the
                  `--samFlagInclude` / `--samFlagExclude` filters
                - **% sincei Duplicates**: duplicates detected by sincei itself
                - **% Duplication**: reads already marked as duplicates by upstream tools
                - **% Singletons**: paired-end reads whose mate did not align as a pair
                - **% Strand / % Motif / % GC / % Low Aligned Fraction**: reads removed by the
                  strand-, motif-, GC-content-, and aligned-fraction filters

                Each row aggregates one sample (the prefix of `Cell_ID` before `::`) by taking
                the median across all of its cells.
                """,
                plot=table.plot(tdata, header, config),
            )

            # General stats: N entries and % Tot. Filtered (per sample, medians across cells)
            gs_headers = {
                "sincei_n_entries": {
                    "title": "N entries",
                    "description": "sincei scFilterStats: median number of sampled reads per cell",
                    "scale": "Greens",
                    "format": "{:,.0f}",
                    "min": 0,
                },
                "sincei_pct_filtered": {
                    "title": "% Tot. Filtered",
                    "description": "sincei scFilterStats: median % of sampled reads filtered for any reason",
                    "suffix": "%",
                    "scale": "OrRd",
                    "min": 0,
                    "max": 100,
                },
            }
            gs_data = {
                k: {
                    "sincei_n_entries": v["N Entries"],
                    "sincei_pct_filtered": v["pct_Filtered"],
                }
                for k, v in tdata.items()
            }
            self.general_stats_addcols(gs_data, gs_headers)

            return len(tdata), len(self.sincei_scFilterStats)

        return 0, len(self.sincei_scFilterStats)

    def parsescFilterStatsFile(self, f):
        reader = csv.DictReader(f["f"], delimiter="\t")
        required = {"Total_sampled", "Filtered", "Blacklisted", "Wrong_motif"}
        if required.difference(set(reader.fieldnames)):
            log.warning(
                f"{f['fn']} was initially flagged as the tabular output from scFilterStats, but that seems to not be the case. Skipping..."
            )
            return dict()

        d = {}
        for row in reader:
            s_name = self.clean_s_name(row["Cell_ID"], f)
            if s_name in d:
                log.debug(f"Replacing duplicate cell_id {s_name}.")
            d[s_name] = {key: row[key] for key in reader.fieldnames}
        return d
