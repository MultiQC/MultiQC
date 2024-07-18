import logging
import os
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Currently supported Longranger pipelines:

    - `wgs`
    - `targeted`
    - `align`

    Usage:

    ```bash
    longranger wgs --fastqs=/path/to/fastq --id=NA12878
    multiqc /path/to/NA12878
    ```

    This module will look for the files `_invocation` and `summary.csv` in the `NA12878` folder, i.e. the output folder of Longranger in this example. The file `summary.csv` is required. If the file `_invocation` is not found the sample will receive a generic name in the MultiQC report (`longranger#1`), instead of `NA12878` or whatever was given by the `--id` parameter.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Long Ranger",
            anchor="longranger",
            href="https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger",
            info="Sample demultiplexing, barcode processing, alignment, quality control, variant calling, phasing, "
            "and structural variant calling.",
            doi="10.1101/gr.234443.118",
        )

        def try_float_lambda(x, func, base):
            try:
                if func == "*":
                    return float(x) * base
                elif func == "/":
                    return float(x) / base
                else:
                    return x
            except Exception:
                return x

        self.headers = {
            "large_sv_calls": {
                "title": "Large SVs",
                "description": "Large structural variants called by Longranger. Not including blacklisted regions.",
                "format": "{:,.0f}",
                "scale": "PuRd",
            },
            "short_deletion_calls": {
                "title": "Short dels",
                "description": "Short deletions called by Longranger.",
                "format": "{:,.0f}",
                "scale": "PuRd",
                "hidden": True,
            },
            "genes_phased_lt_100kb": {
                "title": "genes phased < 100kb",
                "description": "Percentage of genes shorter than 100kb with >1 heterozygous SNP that are phased into a single phase block.",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "suffix": "%",
                "scale": "YlOrRd",
                "hidden": True,
            },
            "longest_phase_block": {
                "title": "Longest phased",
                "description": "Size of the longest phase block, in base pairs",
                "scale": "YlOrRd",
                "modify": lambda x: try_float_lambda(x, "/", 1000000.0),
                "suffix": "Mbp",
                "hidden": True,
            },
            "n50_phase_block": {
                "title": "N50 phased",
                "description": "N50 length of the called phase blocks, in base pairs.",
                "modify": lambda x: try_float_lambda(x, "/", 1000000.0),
                "suffix": "Mbp",
                "scale": "YlOrRd",
                "hidden": True,
            },
            "snps_phased": {
                "title": "SNPs phased",
                "description": "Percentage of called SNPs that were phased.",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "suffix": "%",
                "scale": "PuRd",
                "hidden": True,
            },
            "median_insert_size": {
                "title": "Insert size",
                "description": "Median insert size of aligned read pairs.",
                "format": "{:,.0f}",
                "suffix": "bp",
                "scale": "PuBu",
                "hidden": True,
            },
            "on_target_bases": {
                "title": "On target",
                "description": "Percentage of aligned bases mapped with the target regions in targeted mode. Only bases inside the intervals of target BED file are counted.",
                "suffix": "%",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "scale": "Greens",
            },
            "zero_coverage": {
                "title": "Zero cov",
                "description": "Percentage of non-N bases in the genome with zero coverage.",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "suffix": "%",
                "max": 100.0,
                "min": 0.0,
                "scale": "RdGy-rev",
            },
            "mean_depth": {
                "title": "Depth",
                "description": "Mean read depth, including PCR duplicate reads. In WGS mode, this is measured across the genome; in targeted mode, this is the measure inside targeted regions.",
                "suffix": "X",
                "scale": "PuBu",
            },
            "pcr_duplication": {
                "title": "PCR Dup",
                "description": "Percentage of reads marked as PCR duplicates. To be marked as PCR duplicates, reads must have the same mapping extents on the genome and the same 10x barcode.",
                "suffix": "%",
                "min": 15.0,
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "scale": "RdGy-rev",
                "hidden": True,
            },
            "mapped_reads": {
                "title": "Mapped",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "suffix": "%",
                "description": "Percentage of input reads that were mapped to the reference genome.",
                "scale": "PuBu",
                "hidden": True,
            },
            "number_reads": {
                "title": "M Reads",
                "modify": lambda x: try_float_lambda(x, "/", 1000000.0),
                "description": "Total number of reads supplied to Long Ranger. (millions)",
                "scale": "PuBu",
                "hidden": True,
            },
            "molecule_length_mean": {
                "title": "Mol size",
                "description": "The length-weighted mean input DNA length in base pairs.",
                "modify": lambda x: try_float_lambda(x, "/", 1000.0),
                "suffix": "Kbp",
                "scale": "YlGn",
            },
            "molecule_length_stddev": {
                "title": "Mol stddev",
                "description": "The length-weighted standard deviation of the input DNA length distribution in base pairs.",
                "modify": lambda x: try_float_lambda(x, "/", 1000.0),
                "suffix": "Kbp",
                "scale": "YlGn",
                "hidden": True,
            },
            "n50_linked_reads_per_molecule": {
                "title": "N50 read per mol.",
                "description": "The N50 number of read-pairs per input DNA molecule. Half of read-pairs came from molecules with this many or greater read-pairs.",
                "scale": "BuGn",
                "hidden": True,
            },
            "r1_q30_bases_fract": {
                "title": "% R1 >= Q30",
                "description": "Percentage of bases in R1 with base quality >= 30.",
                "hidden": True,
                "suffix": "%",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "scale": "Purples",
            },
            "r2_q30_bases_fract": {
                "title": "% R2 >= Q30",
                "description": "Percentage of bases in R2 with base quality >= 30.",
                "suffix": "%",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "scale": "Purples",
                "hidden": True,
            },
            "bc_on_whitelist": {
                "title": "Valid BCs",
                "description": "The Percentage of reads that carried a valid 10x barcode sequence.",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "suffix": "%",
                "scale": "BuPu",
                "hidden": True,
            },
            "bc_q30_bases_fract": {
                "title": "BC Q30",
                "description": "Percentage of bases in the barcode with base quality >= 30.",
                "suffix": "%",
                "modify": lambda x: try_float_lambda(x, "*", 100.0),
                "scale": "Purples",
                "hidden": True,
            },
            "bc_mean_qscore": {
                "title": "BC Qscore",
                "description": "The mean base quality value on the barcode bases.",
                "scale": "BuPu",
                "hidden": True,
            },
            "mean_dna_per_gem": {
                "title": "DNA per gem",
                "description": "The average number of base pairs of genomic DNA loaded into each GEM. This metric is based on the observed extents of read-pairs on each molecule.",
                "modify": lambda x: try_float_lambda(x, "/", 1000000.0),
                "suffix": "Mbp",
                "scale": "OrRd",
                "hidden": True,
            },
            "gems_detected": {
                "title": "M Gems",
                "description": "The number of Chromium GEMs that were collected and which generated a non-trivial number of read-pairs. (millions)",
                "modify": lambda x: try_float_lambda(x, "/", 1000000.0),
                "scale": "OrRd",
            },
            "corrected_loaded_mass_ng": {
                "title": "Loaded (corrected)",
                "description": "The estimated number of nanograms of DNA loaded into the input well of the Chromium chip. This metric is calculated by measuring the mean amount of DNA covered by input molecules in each GEM, then multiplying by the ratio of the chip input to the sample volume in each GEM.",
                "suffix": "ng",
                "scale": "RdYlGn",
            },
            "loaded_mass_ng": {
                "title": "Loaded",
                "description": "This metric was found to overestimate the true loading by a factor of 1.6, due primarily to denaturation of the input DNA.",
                "suffix": "ng",
                "scale": "RdYlGn",
            },
            "instrument_ids": {
                "title": "Instrument ID",
                "description": "The list of instrument IDs used to generate the input reads.",
                "scale": False,
                "hidden": True,
            },
            "longranger_version": {
                "title": "Long Ranger Version",
                "description": "The version of the Longranger software used to generate the results.",
                "scale": False,
            },
        }

        # Parse the data
        self.longranger_data = dict()
        self.paths_dict = dict()
        for f in self.find_log_files("longranger/invocation"):
            sid = self.parse_invocation(f["f"])
            self.paths_dict[os.path.basename(f["root"])] = sid

        running_name = 1
        for f in self.find_log_files("longranger/summary"):
            data = self.parse_summary(f["f"])
            updir, _ = os.path.split(f["root"])
            base_updir = os.path.basename(updir)
            sid = f"longranger#{running_name}"
            if base_updir in self.paths_dict.keys():
                sid = self.paths_dict[base_updir]
            else:
                log.debug(f"Did not find _invocation file: {f['fn']}")
                running_name += 1

            self.longranger_data[sid] = data
            self.add_software_version(data["longranger_version"], sid)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.longranger_data = self.ignore_samples(self.longranger_data)

        if len(self.longranger_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.longranger_data.keys())} reports")

        # Write parsed report data to a file
        self.write_data_file(self.longranger_data, "multiqc_longranger")

        # Add a longranger versions column if not all the same
        longranger_versions = set([d["longranger_version"] for d in self.longranger_data.values()])
        version_str = ""
        if len(longranger_versions) == 1:
            version_str = f" All samples were processed using Longranger version {list(longranger_versions)[0]}"
            del self.headers["longranger_version"]

        # Write the table
        config_table = {
            "id": "longranger_table",
            "namespace": "longranger",
            "title": "Long Ranger: Summary Statistics",
        }
        self.add_section(
            name="Run stats",
            anchor="longranger-run-stats",
            description="Statistics gathered from Longranger reports. "
            "There are more columns available but they are hidden by default." + version_str,
            helptext="""Parses the files `summary.csv` and `_invocation` found in the
                    output directory of Longranger. If `_invocation` is not found
                    the sample IDs will be missing and they will be given a running
                    number. E.g., `longranger#1` and `longranger#2`.""",
            plot=table.plot(self.longranger_data, self.headers, config_table),
        )

        # Bar plot of phasing stats
        phase_pdata = {}
        snps_phased_pct = {}
        genes_phased_pct = {}
        for s_name in self.longranger_data:
            try:
                phase_pdata[s_name] = {
                    "longest_phase_block": float(self.longranger_data[s_name]["longest_phase_block"]),
                    "n50_phase_block": float(self.longranger_data[s_name]["n50_phase_block"]),
                }
            except Exception:
                pass
            try:
                snps_phased_pct[s_name] = {
                    "snps_phased_pct": float(self.longranger_data[s_name]["snps_phased"]) * 100.0
                }
            except Exception:
                pass
            try:
                genes_phased_pct[s_name] = {
                    "genes_phased_pct": float(self.longranger_data[s_name]["genes_phased_lt_100kb"]) * 100.0
                }
            except Exception:
                pass
        phase_plot_cats = [dict(), dict(), dict()]
        phase_plot_cats[0]["longest_phase_block"] = {"name": "Longest Phase Block"}
        phase_plot_cats[0]["n50_phase_block"] = {"name": "N50 of Phase Blocks"}
        phase_plot_cats[1]["snps_phased_pct"] = {"name": "% SNPs Phased"}
        phase_plot_cats[2]["genes_phased_pct"] = {"name": "% Genes < 100kbp in a single phase block"}
        if len(phase_pdata) > 0:
            self.add_section(
                name="Phasing",
                anchor="longranger-phasing",
                description="Phasing performance from Long Ranger. Genes are only considered if &le; 100kbp in length and with at least one heterozygous SNP.",
                helptext="""
                        * Longest phased
                            * Size of the longest phase block, in base pairs
                        * N50 phased
                            * N50 length of the called phase blocks, in base pairs.
                        * % SNPs phased
                            * Percentage of called SNPs that were phased.
                        * % Genes Phased
                            * Percentage of genes shorter than 100kb with >1 heterozygous SNP that are phased into a single phase block.
                        """,
                plot=bargraph.plot(
                    [phase_pdata, snps_phased_pct, genes_phased_pct],
                    phase_plot_cats,
                    {
                        "id": "longranger-phasing-plot",
                        "title": "Long Ranger: Phasing Statistics",
                        "data_labels": [
                            {"name": "N50 Phased", "ylab": "N50 of called phase blocks (bp)"},
                            {"name": "% SNPs Phased", "ylab": "% SNPs Phased", "ymax": 100},
                            {"name": "% Genes Phased", "ylab": "% Genes Phased", "ymax": 100},
                        ],
                        "cpswitch": False,
                        "stacking": None,
                        "ylab": "N50 of called phase blocks (bp)",
                    },
                ),
            )

        # Bar plot of mapping statistics
        mapping_counts_data = {}
        for s_name in self.longranger_data:
            mapped_reads = float(self.longranger_data[s_name]["number_reads"]) * float(
                self.longranger_data[s_name]["mapped_reads"]
            )
            unmapped_reads = float(self.longranger_data[s_name]["number_reads"]) - mapped_reads
            dup_reads = mapped_reads * float(self.longranger_data[s_name]["pcr_duplication"])
            unique_reads = mapped_reads - dup_reads
            mapping_counts_data[s_name] = {
                "unique_reads": unique_reads,
                "dup_reads": dup_reads,
                "unmapped_reads": unmapped_reads,
            }
        mapping_counts_cats = dict()
        mapping_counts_cats["unique_reads"] = {"name": "Uniquely Aligned Reads", "color": "#437bb1"}
        mapping_counts_cats["dup_reads"] = {"name": "PCR Duplicate Aligned Reads", "color": "#7cb5ec"}
        mapping_counts_cats["unmapped_reads"] = {"name": "Unaligned Reads", "color": "#7f0000"}
        self.add_section(
            name="Alignment",
            anchor="longranger-alignment",
            description="Long Ranger alignment against the reference genome. To be marked as PCR duplicates, reads must have the same mapping extents on the genome and the same 10x barcode.",
            plot=bargraph.plot(
                mapping_counts_data,
                mapping_counts_cats,
                {
                    "id": "longranger-alignment-plot",
                    "title": "Long Ranger: Alignment Statistics",
                    "ylab": "Reads Counts",
                    "cpswitch_counts_label": "Read Counts",
                },
            ),
        )

    @staticmethod
    def parse_invocation(content):
        sid_pat = re.compile('    sample_id = "(.*)",')

        sid = None
        for line in content.splitlines():
            sid_m = re.match(sid_pat, line)
            if sid_m is not None:
                sid = sid_m.groups()[0]
        return sid

    @staticmethod
    def parse_summary(content):
        out_dict = dict()
        lines = content.splitlines()
        data = list(zip(lines[0].strip().split(","), lines[1].strip().split(",")))
        for i, j in data:
            if j != "":
                try:
                    out_dict[i] = float(j)
                except ValueError:
                    out_dict[i] = j

        return out_dict
