import logging
from collections import defaultdict
from typing import Dict

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.plots.plotly.bar import BarPlotConfig

# Initialise the logger
log = logging.getLogger(__name__)


def parse_samtools_idxstats(module):
    """Find Samtools idxstats logs and parse their data"""

    module.samtools_idxstats = dict()
    for f in module.find_log_files("samtools/idxstats"):
        parsed_data = parse_single_report(f["f"])
        if len(parsed_data) > 0:
            if f["s_name"] in module.samtools_idxstats:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="idxstats")
            module.samtools_idxstats[f["s_name"]] = parsed_data

    # Filter to strip out ignored sample names
    module.samtools_idxstats = module.ignore_samples(module.samtools_idxstats)

    if len(module.samtools_idxstats) == 0:
        return 0

    # Write parsed report data to a file (restructure first)
    module.write_data_file(module.samtools_idxstats, "multiqc_samtools_idxstats")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Prep the data for the plots
    keys = list()
    pdata: Dict = dict()
    pdata_norm: Dict = dict()
    pdata_obs_exp: Dict = dict()
    xy_counts: Dict = dict()
    # Count the total mapped reads for every chromosome
    chrs_mapped: Dict[str, int] = defaultdict(lambda: 0)
    sample_mapped: Dict[str, int] = defaultdict(lambda: 0)
    total_mapped = 0
    # Cutoff, can be customised in config
    cutoff = float(getattr(config, "samtools_idxstats_fraction_cutoff", 0.001))
    if cutoff != 0.001:
        log.info(f"Setting idxstats cutoff to: {cutoff * 100.0}%")
    for s_name in module.samtools_idxstats:
        for chrom in module.samtools_idxstats[s_name]:
            chrs_mapped[chrom] += module.samtools_idxstats[s_name][chrom][0]
            sample_mapped[s_name] += module.samtools_idxstats[s_name][chrom][0]
            total_mapped += module.samtools_idxstats[s_name][chrom][0]
    req_reads = float(total_mapped) * cutoff
    chr_always = getattr(config, "samtools_idxstats_always", [])
    if len(chr_always) > 0:
        log.info(f"Trying to include these chromosomes in idxstats: {', '.join(chr_always)}")
    chr_ignore = getattr(config, "samtools_idxstats_ignore", [])
    if len(chr_ignore) > 0:
        log.info(f"Excluding these chromosomes from idxstats: {', '.join(chr_ignore)}")
    xchr: str = getattr(config, "samtools_idxstats_xchr", "")
    if xchr:
        log.info(f'Using "{xchr}" as X chromosome name')
    ychr: str = getattr(config, "samtools_idxstats_ychr", "")
    if ychr:
        log.info(f'Using "{ychr}" as Y chromosome name')
    # Go through again and collect all the keys that have enough counts
    # Also get the X/Y counts if we find them
    for s_name in module.samtools_idxstats:
        x_count = False
        y_count = False
        for chrom in module.samtools_idxstats[s_name]:
            if float(chrs_mapped[chrom]) > req_reads or chrom in chr_always:
                if chrom not in chr_ignore and chrom not in keys:
                    keys.append(chrom)
            # Collect X and Y counts if we have them
            mapped = module.samtools_idxstats[s_name][chrom][0]
            if xchr:
                if str(xchr) == str(chrom):
                    x_count = mapped
            else:
                if chrom.lower() == "x" or chrom.lower() == "chrx":
                    x_count = mapped
            if ychr:
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
    for s_name in module.samtools_idxstats:
        pdata[s_name] = dict()
        pdata_norm[s_name] = dict()
        pdata_obs_exp[s_name] = dict()
        genome_size = float(sum([stats[1] for stats in module.samtools_idxstats[s_name].values()]))
        for k in keys:
            try:
                pdata[s_name][k] = module.samtools_idxstats[s_name][k][0]
                pdata_norm[s_name][k] = float(module.samtools_idxstats[s_name][k][0]) / sample_mapped[s_name]
                chrom_size = float(module.samtools_idxstats[s_name][k][1])
                expected_count = (chrom_size / genome_size) * float(sample_mapped[s_name])
                pdata_obs_exp[s_name][k] = float(pdata[s_name][k]) / expected_count
            except (KeyError, ZeroDivisionError):
                pdata[s_name][k] = 0
                pdata_norm[s_name][k] = 0
                pdata_obs_exp[s_name][k] = 0

    # X/Y ratio plot
    if len(xy_counts) > 0:
        xy_keys = {
            "x": {"name": xchr if xchr else "Chromosome X"},
            "y": {"name": ychr if ychr else "Chromosome Y"},
        }
        module.add_section(
            name="XY counts",
            anchor="samtools-idxstats-xy-counts",
            plot=bargraph.plot(
                xy_counts,
                xy_keys,
                BarPlotConfig(
                    id="samtools-idxstats-xy-plot",
                    title="Samtools: idxstats: chrXY mapped reads",
                    ylab="Percent of X+Y Reads",
                    cpswitch_counts_label="Number of Reads",
                    cpswitch_percent_label="Percent of X+Y Reads",
                    cpswitch_c_active=False,
                ),
            ),
        )

    # Mapped reads per chr line plot
    pconfig = {
        "id": "samtools-idxstats-mapped-reads-plot",
        "title": "Samtools: idxstats: Mapped reads per contig",
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
    module.add_section(
        name="Mapped reads per contig",
        anchor="samtools-idxstats",
        description="The <code>samtools idxstats</code> tool counts the number of mapped reads per chromosome / contig. "
        + f"Chromosomes with &lt; {cutoff * 100}% of the total aligned reads are omitted from this plot.",
        plot=linegraph.plot([pdata_norm, pdata_obs_exp, pdata], pconfig),
    )

    # Return the number of logs that were found
    return len(module.samtools_idxstats)


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
