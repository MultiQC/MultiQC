""" MultiQC submodule to parse output from Samtools idxstats """

import logging
from collections import defaultdict

from multiqc import config
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class IdxstatsReportMixin:
    def parse_samtools_idxstats(self):
        """Find Samtools idxstats logs and parse their data"""

        self.samtools_idxstats = dict()
        for f in self.find_log_files("samtools/idxstats"):
            parsed_data = parse_single_report(f["f"])
            if len(parsed_data) > 0:
                if f["s_name"] in self.samtools_idxstats:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="idxstats")
                self.samtools_idxstats[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.samtools_idxstats = self.ignore_samples(self.samtools_idxstats)

        if len(self.samtools_idxstats) == 0:
            return 0

        # Write parsed report data to a file (restructure first)
        self.write_data_file(self.samtools_idxstats, "multiqc_samtools_idxstats")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Prep the data for the plots
        keys = list()
        pdata = dict()
        pdata_norm = dict()
        pdata_obs_exp = dict()
        xy_counts = dict()
        # Count the total mapped reads for every chromosome
        chrs_mapped = defaultdict(lambda: 0)
        sample_mapped = defaultdict(lambda: 0)
        total_mapped = 0
        # Cutoff, can be customised in config
        cutoff = float(getattr(config, "samtools_idxstats_fraction_cutoff", 0.001))
        if cutoff != 0.001:
            log.info(f"Setting idxstats cutoff to: {cutoff * 100.0}%")
        for s_name in self.samtools_idxstats:
            for chrom in self.samtools_idxstats[s_name]:
                chrs_mapped[chrom] += self.samtools_idxstats[s_name][chrom][0]
                sample_mapped[s_name] += self.samtools_idxstats[s_name][chrom][0]
                total_mapped += self.samtools_idxstats[s_name][chrom][0]
        req_reads = float(total_mapped) * cutoff
        chr_always = getattr(config, "samtools_idxstats_always", [])
        if len(chr_always) > 0:
            log.info(f"Trying to include these chromosomes in idxstats: {', '.join(chr_always)}")
        chr_ignore = getattr(config, "samtools_idxstats_ignore", [])
        if len(chr_ignore) > 0:
            log.info(f"Excluding these chromosomes from idxstats: {', '.join(chr_ignore)}")
        xchr = getattr(config, "samtools_idxstats_xchr", False)
        if xchr:
            log.info(f'Using "{xchr}" as X chromosome name')
        ychr = getattr(config, "samtools_idxstats_ychr", False)
        if ychr:
            log.info(f'Using "{ychr}" as Y chromosome name')
        # Go through again and collect all of the keys that have enough counts
        # Also get the X/Y counts if we find them
        for s_name in self.samtools_idxstats:
            x_count = False
            y_count = False
            for chrom in self.samtools_idxstats[s_name]:
                if float(chrs_mapped[chrom]) > req_reads or chrom in chr_always:
                    if chrom not in chr_ignore and chrom not in keys:
                        keys.append(chrom)
                # Collect X and Y counts if we have them
                mapped = self.samtools_idxstats[s_name][chrom][0]
                if xchr is not False:
                    if str(xchr) == str(chrom):
                        x_count = mapped
                else:
                    if chrom.lower() == "x" or chrom.lower() == "chrx":
                        x_count = mapped
                if ychr is not False:
                    if str(ychr) == str(chrom):
                        y_count = mapped
                else:
                    if chrom.lower() == "y" or chrom.lower() == "chry":
                        y_count = mapped
            # Only save these counts if we have both x and y
            if x_count and y_count:
                xy_counts[s_name] = {"x": x_count, "y": y_count}
        # Ok, one last time. We have the chromosomes that we want to plot,
        # now collect the counts
        for s_name in self.samtools_idxstats:
            pdata[s_name] = dict()
            pdata_norm[s_name] = dict()
            pdata_obs_exp[s_name] = dict()
            genome_size = float(sum([stats[1] for stats in self.samtools_idxstats[s_name].values()]))
            for k in keys:
                try:
                    pdata[s_name][k] = self.samtools_idxstats[s_name][k][0]
                    pdata_norm[s_name][k] = float(self.samtools_idxstats[s_name][k][0]) / sample_mapped[s_name]
                    chrom_size = float(self.samtools_idxstats[s_name][k][1])
                    expected_count = (chrom_size / genome_size) * float(sample_mapped[s_name])
                    pdata_obs_exp[s_name][k] = float(pdata[s_name][k]) / expected_count
                except (KeyError, ZeroDivisionError):
                    pdata[s_name][k] = 0
                    pdata_norm[s_name][k] = 0
                    pdata_obs_exp[s_name][k] = 0

        # X/Y ratio plot
        if len(xy_counts) > 0:
            xy_keys = dict()
            xy_keys["x"] = {"name": xchr if xchr else "Chromosome X"}
            xy_keys["y"] = {"name": ychr if ychr else "Chromosome Y"}
            pconfig = {
                "id": "samtools-idxstats-xy-plot",
                "title": "Samtools idxstats: chrXY mapped reads",
                "ylab": "Percent of X+Y Reads",
                "cpswitch_counts_label": "Number of Reads",
                "cpswitch_percent_label": "Percent of X+Y Reads",
                "cpswitch_c_active": False,
            }
            self.add_section(
                name="XY counts",
                anchor="samtools-idxstats-xy-counts",
                plot=bargraph.plot(xy_counts, xy_keys, pconfig),
            )

        # Mapped reads per chr line plot
        pconfig = {
            "id": "samtools-idxstats-mapped-reads-plot",
            "title": "Samtools idxstats: Mapped reads per contig",
            "ylab": "# mapped reads",
            "xlab": "Chromosome name",
            "logswitch": True,
            "categories": True,
            "tt_label": "<strong>{point.category}:</strong> {point.y:.2f}",
            "data_labels": [
                {"name": "Normalised Counts", "ylab": "Fraction of total count"},
                {"name": "Observed over Expected Counts", "ylab": "log10 ( Observed over expected counts )"},
                {"name": "Raw Counts", "ylab": "# mapped reads"},
            ],
        }
        self.add_section(
            name="Mapped reads per contig",
            anchor="samtools-idxstats",
            description="The <code>samtools idxstats</code> tool counts the number of mapped reads per chromosome / contig. "
            + f"Chromosomes with &lt; {cutoff * 100}% of the total aligned reads are omitted from this plot.",
            plot=linegraph.plot([pdata_norm, pdata_obs_exp, pdata], pconfig),
        )

        # Return the number of logs that were found
        return len(self.samtools_idxstats)


def parse_single_report(f):
    """Parse a samtools idxstats idxstats"""

    parsed_data = dict()
    for line in f.splitlines():
        s = line.split("\t")
        try:
            parsed_data[s[0]] = [int(s[2]), int(s[1])]
        except (IndexError, ValueError):
            pass
    return parsed_data
