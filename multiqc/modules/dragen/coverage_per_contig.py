import logging
import re
from collections import defaultdict
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class DragenCoveragePerContig(BaseMultiqcModule):
    def add_coverage_per_contig(self):
        perchrom_data_by_phenotype_by_sample: Dict[str, Dict] = defaultdict(dict)

        for f in self.find_log_files("dragen/wgs_contig_mean_cov"):
            perchrom_data_by_phenotype = parse_wgs_contig_mean_cov(f)
            s_name = f["s_name"]
            if s_name in perchrom_data_by_phenotype_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="wgs_contig_mean_cov")
            perchrom_data_by_phenotype_by_sample[s_name].update(perchrom_data_by_phenotype)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Filter to strip out ignored sample names:
        perchrom_data_by_phenotype_by_sample = self.ignore_samples(perchrom_data_by_phenotype_by_sample)

        # Merge tumor and normal data:
        perchrom_data_by_sample: Dict[str, Dict] = defaultdict(dict)
        for sn in perchrom_data_by_phenotype_by_sample:
            for phenotype in perchrom_data_by_phenotype_by_sample[sn]:
                new_sn = sn
                if phenotype == "normal":
                    new_sn = sn + "_normal"
                perchrom_data_by_sample[new_sn] = perchrom_data_by_phenotype_by_sample[sn][phenotype]

        if not perchrom_data_by_sample:
            return set()

        # Data is in wrong format for writing to file
        # self.write_data_file(perchrom_data_by_sample, "dragen_cov_contig")

        main_contigs_by_sample = {sn: data[0] for sn, data in perchrom_data_by_sample.items()}
        other_contigs_by_sample = {sn: data[1] for sn, data in perchrom_data_by_sample.items()}

        self.add_section(
            name="Coverage per contig",
            anchor="dragen-coverage-per-contig",
            description="""
            Average coverage per contig or chromosome.
            Calculated as the number of bases (excluding duplicate marked reads, reads
            with MAPQ=0, and clipped bases), divided by the length of the contig or
            (if a target bed is used) the total length of the target region spanning that contig.
            """,
            plot=linegraph.plot(
                main_contigs_by_sample,
                pconfig={
                    "id": "dragen_coverage_per_contig",
                    "title": "Dragen: Average coverage per contig"
                    + (" (main contigs)" if other_contigs_by_sample else ""),
                    "ylab": "Average coverage",
                    "xlab": "Region",
                    "categories": True,
                    "tt_label": "<b>{point.x}</b>: {point.y:.1f}x",
                },
            ),
        )

        if other_contigs_by_sample:
            self.add_section(
                name="Coverage per contig (non-main)",
                anchor="dragen-coverage-per-nonmain-contig",
                description="""
                Non-main contigs:
                unlocalized (*_random), unplaced (chrU_*), alts (*_alt), mitochondria (chrM), EBV, HLA.
                Zoom in to see more contigs as all labels don\'t fit the screen.
                """,
                plot=linegraph.plot(
                    other_contigs_by_sample,
                    pconfig={
                        "id": "dragen_coverage_per_non_main_contig",
                        "title": "Dragen: Average coverage of non-main contigs",
                        "ylab": "Average coverage",
                        "xlab": "Region",
                        "categories": True,
                        "tt_label": "<b>{point.x}</b>: {point.y:.1f}x",
                    },
                ),
            )

        return perchrom_data_by_sample.keys()


def parse_wgs_contig_mean_cov(f):
    """
    The Contig Mean Coverage report generates a _contig_mean_cov.csv file, which contains the estimated coverage for
    all contigs, and an autosomal estimated coverage. The file includes the following three columns

    1. Contig name
    2. Number of bases aligned to that contig, which excludes bases from duplicate marked reads, reads with MAPQ=0,
       and clipped bases.
    3. Estimated coverage, as follows: <number of bases aligned to the contig (ie, Col2)> divided by <length of the
       contig or (if a target bed is used) the total length of the target region spanning that contig>.

    T_SRR7890936_50pc.wgs_contig_mean_cov_normal.csv
    T_SRR7890936_50pc.wgs_contig_mean_cov_tumor.csv

    chr1,11292297134,48.9945
    chr10,6482885699,48.6473
    ...
    chrUn_GL000218v1,20750824,128.77
    chrX,3590295769,23.1792
    chrY,42229820,1.5987
    chrY_KI270740v1_random,0,0
    Autosomal regions ,130912665915,47.4953

    A histogram or a plot like in mosdepth, with each chrom in X axis.
    Two versions: for main contigs (chr1..chrY), and including non-main contigs (chrUn_*, *_random, HLA-*, *_alt, chrM)
    """

    main_contig_perchrom_data = dict()
    other_contig_perchrom_data = dict()

    for line in f["f"].splitlines():
        chrom, bases, depth = line.split(",")
        chrom = chrom.strip()
        depth = float(depth)
        # skipping unplaced and alternative contigs, as well as the mitochondria (might attract 100 times more coverage
        # than human chromosomes):
        if (
            chrom.startswith("chrUn_")
            or chrom.endswith("_random")
            or chrom.endswith("_alt")
            or chrom == "chrM"
            or chrom == "MT"
            or chrom == "chrEBV"
            or chrom.startswith("HLA-")
        ):
            other_contig_perchrom_data[chrom] = depth
        else:
            main_contig_perchrom_data[chrom] = depth

    def chrom_order(chrom):
        if chrom == "Autosomal regions":
            # "Autosomal regions" average coverage goes right after all the autosomal chromosomes
            return 0
        try:
            # autosomal chromosomes go first, thus getting a negative order
            return int(chrom.replace("chr", "")) - len(main_contig_perchrom_data)
        except ValueError:
            # sex and other chromosomes go in the end
            return 1

    main_contig_perchrom_data = dict(
        sorted(
            main_contig_perchrom_data.items(),
            key=lambda key_val: chrom_order(key_val[0]),
        )
    )
    other_contig_perchrom_data = dict(
        sorted(
            other_contig_perchrom_data.items(),
            key=lambda key_val: chrom_order(key_val[0]),
        )
    )

    m = re.search(r"(tumor|normal).csv", f["fn"])
    if m:
        phenotype = m.group(1)
    else:
        phenotype = "unknown"
    return {phenotype: [main_contig_perchrom_data, other_contig_perchrom_data]}
