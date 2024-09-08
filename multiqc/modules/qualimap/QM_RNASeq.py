"""MultiQC Submodule to parse output from Qualimap RNASeq"""

import logging
import os
import re
from typing import Dict

from multiqc import config, BaseMultiqcModule
from multiqc.modules.qualimap import parse_numerals, get_s_name
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find Qualimap RNASeq reports and parse their data"""
    genome_results: Dict = dict()

    int_metrics = {
        "read pairs aligned": "reads_aligned",
        "reads aligned": "reads_aligned",
        "total alignments": "total_alignments",
        "non-unique alignments": "non_unique_alignments",
        "aligned to genes": "reads_aligned_genes",
        "ambiguous alignments": "ambiguous_alignments",
        "not aligned": "not_aligned",
        "exonic": "reads_aligned_exonic",
        "intronic": "reads_aligned_intronic",
        "intergenic": "reads_aligned_intergenic",
        "overlapping exon": "reads_aligned_overlapping_exon",
    }
    rate_metrics = {
        "5' bias": "5_bias",
        "3' bias": "3_bias",
        "SSP estimation (fwd/rev)": "ssp_estimation",
        "5'-3' bias": "5_3_bias",
    }

    value_regex = re.compile(r"\s+[\d,\.\xa0]+\s+")
    for f in module.find_log_files("qualimap/rnaseq/rnaseq_results"):
        preparsed_d = dict()
        for line in f["f"].splitlines():
            if "=" in line:
                key, val = line.split("=", 1)
                m = re.search(value_regex, val)
                if m:
                    val = m.group(0)
                    val = re.sub(r"\xa0", "", val)
                key = key.strip()
                preparsed_d[key] = val.strip()

        # Check we have an input filename
        if "bam file" not in preparsed_d:
            log.debug(f"Couldn't find an input filename in genome_results file {f['fn']}")
            return 0

        s_name = module.clean_s_name(preparsed_d["bam file"], f)
        if s_name in genome_results:
            log.debug(f"Duplicate genome results sample name found! Overwriting: {s_name}")

        module.add_data_source(f, s_name=s_name, section="rna_genome_results")

        d = parse_numerals(
            preparsed_d,
            float_metrics={},
            int_metrics=int_metrics,
            rate_metrics=rate_metrics,
            fpath=os.path.join(f["root"], f["fn"]),
        )

        genome_results[s_name] = d

    module.general_stats_addcols(
        genome_results,
        headers={
            "5_3_bias": {
                "title": "5'-3' bias",
                "format": "{:,.2f}",
            },
            "reads_aligned": {
                "title": f"{config.read_count_prefix} Aligned",
                "description": f"Reads Aligned ({config.read_count_desc})",
                "min": 0,
                "scale": "RdBu",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
        },
        namespace="RNASeq",
    )

    # Coverage profile
    cov_hist: Dict = dict()
    for f in module.find_log_files("qualimap/rnaseq/coverage", filehandles=True):
        s_name = get_s_name(module, f)
        # Save results
        if s_name in cov_hist:
            log.debug(f"Duplicate coverage histogram sample name found! Overwriting: {s_name}")

        module.add_data_source(f, s_name=s_name, section="rna_coverage_histogram")

        d = dict()
        for line in f["f"]:
            if line.startswith("#"):
                continue
            coverage, count = line.split(None, 1)
            coverage = int(round(float(coverage.replace(",", "."))))
            count = float(count)
            d[coverage] = count

        if len(d) == 0:
            log.debug(f"Couldn't parse contents of coverage histogram file {f['fn']}")
            return 0

        cov_hist[s_name] = d

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Filter to strip out ignored sample names
    genome_results = module.ignore_samples(genome_results)
    cov_hist = module.ignore_samples(cov_hist)

    # Plots
    # Genomic Origin Bar Graph
    # NB: Ignore 'Overlapping Exon' in report - these make the numbers add up to > 100%
    if len(genome_results) > 0:
        # Check that we have anything to plot
        if (
            sum(
                entry[key]
                for entry in genome_results.values()
                for key in ["reads_aligned_exonic", "reads_aligned_intronic", "reads_aligned_intronic"]
            )
            > 0
        ):
            # Write data to file
            module.write_data_file(genome_results, "qualimap_rnaseq_genome_results")

            gorigin_cats = {
                "reads_aligned_exonic": {"name": "Exonic"},
                "reads_aligned_intronic": {"name": "Intronic"},
                "reads_aligned_intergenic": {"name": "Intergenic"},
            }
            gorigin_pconfig = {
                "id": "qualimap_genomic_origin",
                "title": "QualiMap: RNAseq: Genomic Origin",
                "ylab": "Number of reads",
                "cpswitch_c_active": False,
            }
            genomic_origin_helptext = """
            There are currently three main approaches to map reads to transcripts in an
            RNA-seq experiment: mapping reads to a reference genome to identify expressed
            transcripts that are annotated (and discover those that are unknown), mapping
            reads to a reference transcriptome, and <i>de novo</i> assembly of transcript
            sequences (<a href="https://doi.org/10.1186/s13059-016-0881-8"
            target="_blank">Conesa et al. 2016</a>).

            For RNA-seq QC analysis, QualiMap can be used to assess alignments produced by
            the first of these approaches. For input, it requires a GTF annotation file
            along with a reference genome, which can be used to reconstruct the exon
            structure of known transcripts. This allows mapped reads to be grouped by
            whether they originate in an exonic region (for QualiMap, this may include
            5&#8242; and 3&#8242; UTR regions as well as protein-coding exons), an intron,
            or an intergenic region (see the <a href="http://qualimap.bioinfo.cipf.es/doc_html/index.html"
            target="_blank">Qualimap 2 documentation</a>).

            The inferred genomic origins of RNA-seq reads are presented here as a bar graph
            showing either the number or percentage of mapped reads in each read dataset
            that have been assigned to each type of genomic region. This graph can be used
            to assess the proportion of useful reads in an RNA-seq experiment. That
            proportion can be reduced by the presence of intron sequences, especially if
            depletion of ribosomal RNA was used during sample preparation (<a href="https://doi.org/10.1038/nrg3642"
            target="_blank">Sims et al. 2014</a>). It can also be reduced by off-target
            transcripts, which are detected in greater numbers at the sequencing depths
            needed to detect poorly-expressed transcripts (<a href="https://doi.org/10.1101/gr.124321.111"
            target="_blank">Tarazona et al. 2011</a>)."""
            module.add_section(
                name="Genomic origin of reads",
                anchor="qualimap-reads-genomic-origin",
                description="Classification of mapped reads as originating in exonic, intronic or intergenic regions. These can be displayed as either the number or percentage of mapped reads.",
                helptext=genomic_origin_helptext,
                plot=bargraph.plot(genome_results, gorigin_cats, gorigin_pconfig),
            )
        else:
            log.warning("Found zero aligned reads. Skipping 'Genomic origin of reads' plot.")

    if len(cov_hist) > 0:
        # Write data to file
        module.write_data_file(cov_hist, "qualimap_rnaseq_cov_hist")

        # Make a normalised percentage version of the coverage data
        cov_hist_percent: Dict = dict()
        for s_name in cov_hist:
            cov_hist_percent[s_name] = dict()
            total = sum(cov_hist[s_name].values())
            if total == 0:
                for k, v in cov_hist[s_name].items():
                    cov_hist_percent[s_name][k] = 0.0
            else:
                for k, v in cov_hist[s_name].items():
                    cov_hist_percent[s_name][k] = (v / total) * 100.0

        coverage_profile_helptext = """
        There are currently three main approaches to map reads to transcripts in an
        RNA-seq experiment: mapping reads to a reference genome to identify expressed
        transcripts that are annotated (and discover those that are unknown), mapping
        reads to a reference transcriptome, and <i>de novo</i> assembly of transcript
        sequences (<a href="https://doi.org/10.1186/s13059-016-0881-8"
        target="_blank">Conesa et al. 2016</a>).

        For RNA-seq QC analysis, QualiMap can be used to assess alignments produced by
        the first of these approaches. For input, it requires a GTF annotation file
        along with a reference genome, which can be used to reconstruct the exon
        structure of known transcripts. QualiMap uses this information to calculate the
        depth of coverage along the length of each annotated transcript. For a set of
        reads mapped to a transcript, the depth of coverage at a given base position is
        the number of high-quality reads that map to the transcript at that position
        (<a href="https://doi.org/10.1038/nrg3642" target="_blank">Sims et al. 2014</a>).

        QualiMap calculates coverage depth at every base position of each annotated
        transcript. To enable meaningful comparison between transcripts, base positions
        are rescaled to relative positions expressed as percentage distance along each
        transcript (*0%, 1%, &#8230;, 99%*). For the set of transcripts with at least
        one mapped read, QualiMap plots the _cumulative mapped-read depth_ (y-axis) at
        each relative transcript position (x-axis). This plot shows the gene coverage
        profile across all mapped transcripts for each read dataset. It provides a
        visual way to assess positional biases, such as an accumulation of mapped reads
        at the 3&#8242; end of transcripts, which may indicate poor RNA quality in the
        original sample (<a href="https://doi.org/10.1186/s13059-016-0881-8"
        target="_blank">Conesa et al. 2016</a>).

        The _Normalised_ plot is calculated by MultiQC to enable comparison of samples
        with varying sequencing depth. The _cumulative mapped-read depth_ at each
        position across the averaged transcript position are divided by the total for
        that sample across the entire averaged transcript.
        """
        pconfig = {
            "id": "qualimap_gene_coverage_profile",
            "title": "QualiMap: RNAseq: Coverage Profile Along Genes (total)",
            "ylab": "Cumulative mapped-read depth",
            "xlab": "Transcript Position (%)",
            "ymin": 0,
            "xmin": 0,
            "xmax": 100,
            "tt_label": "<b>{point.x}%</b>: {point.y:.2f}",
            "data_labels": [
                {"name": "Counts", "ylab": "Cumulative mapped-read depth"},
                {"name": "Normalised", "ylab": "Percentage total cumulative mapped-read depth"},
            ],
        }
        module.add_section(
            name="Gene Coverage Profile",
            anchor="qualimap-genome-fraction-coverage",
            description="Mean distribution of coverage depth across the length of all mapped transcripts.",
            helptext=coverage_profile_helptext,
            plot=linegraph.plot([cov_hist, cov_hist_percent], pconfig),
        )

    # Return the number of reports we found
    return len(genome_results.keys())
