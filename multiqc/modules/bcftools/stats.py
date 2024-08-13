"""MultiQC submodule to parse output from Bcftools stats"""

import logging
import re
from typing import Dict

from multiqc import config, BaseMultiqcModule
from multiqc.plots import bargraph, linegraph, table
from multiqc.plots.plotly.bar import BarPlotConfig
from multiqc.plots.plotly.line import LinePlotConfig

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"# This file was produced by bcftools stats \(([\d\.]+)"
HTSLIB_REGEX = r"\+htslib-([\d\.]+)"


def parse_bcftools_stats(module: BaseMultiqcModule) -> int:
    """
    Find bcftools stats logs and parse their data
      Bcftools stats reports contain 'sets' of data, which can
      have multiple vcf files each (but usually don't). Here,
      we treat each 'set' as a MultiQC sample, taking the first
      input filename for each set as the name.
    """
    collapse_complementary = getattr(config, "bcftools", {}).get("collapse_complementary_changes", False)
    if collapse_complementary:
        types = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T"]
    else:
        types = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]

    bcftools_stats: Dict = dict()
    bcftools_stats_indels: Dict = dict()
    bcftools_stats_vqc_snp: Dict = dict()
    bcftools_stats_sample_variants: Dict = dict()
    bcftools_stats_sample_tstv: Dict = dict()
    bcftools_stats_sample_singletons: Dict = dict()
    bcftools_stats_sample_depth: Dict = dict()
    bcftools_stats_vqc_transi: Dict = dict()
    bcftools_stats_vqc_transv: Dict = dict()
    bcftools_stats_vqc_indels: Dict = dict()
    bcftools_stats_depth_data: Dict = dict()
    for f in module.find_log_files("bcftools/stats"):
        s_names = list()
        for line in f["f"].splitlines():
            # Get version number from file contents
            if line.startswith("# This file was produced by bcftools stats"):
                # Look for BCFtools version
                version_match = re.search(VERSION_REGEX, line)
                if version_match is None:
                    continue

                # Add BCFtools version
                bcftools_version = version_match.group(1)
                module.add_software_version(bcftools_version, f["s_name"])

                # Look for HTSlib version
                htslib_version_match = re.search(HTSLIB_REGEX, line)
                if htslib_version_match is None:
                    continue

                # Add HTSlib version if different from BCFtools version
                htslib_version = htslib_version_match.group(1)
                if htslib_version != bcftools_version:
                    module.add_software_version(htslib_version, f["s_name"], "HTSlib")

            s = line.split("\t")
            # Get the sample names - one per 'set'
            if s[0] == "ID":
                s_name = module.clean_s_name(s[2], f)
                s_names.append(s_name)
                if s_name in bcftools_stats:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                module.add_data_source(f, s_name, section="stats")
                bcftools_stats[s_name] = dict()
                bcftools_stats_indels[s_name] = dict()
                bcftools_stats_sample_variants[s_name] = dict()
                bcftools_stats_sample_tstv[s_name] = dict()
                bcftools_stats_sample_singletons[s_name] = dict()
                bcftools_stats_sample_depth[s_name] = dict()
                bcftools_stats_vqc_snp[s_name] = dict()
                bcftools_stats_vqc_transi[s_name] = dict()
                bcftools_stats_vqc_transv[s_name] = dict()
                bcftools_stats_vqc_indels[s_name] = dict()
                bcftools_stats_depth_data[s_name] = {}
                bcftools_stats_indels[s_name][0] = None  # Avoid joining line across missing 0

            # Parse key stats
            if s[0] == "SN" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                field = s[2].strip()[:-1]
                field = field.replace(" ", "_")
                value_int = int(s[3].strip())
                bcftools_stats[s_name][field] = value_int

            # Parse transitions/transversions stats
            if s[0] == "TSTV" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                fields = ["ts", "tv", "tstv", "ts_1st_ALT", "tv_1st_ALT", "tstv_1st_ALT"]
                for i, field in enumerate(fields):
                    value_str = s[i + 2].strip()
                    if "tstv" in field:
                        value = float(value_str)
                    else:
                        value = int(value_str)
                    bcftools_stats[s_name][field] = value

            # Parse substitution types
            if s[0] == "ST" and len(s_names) > 0:
                s_name = s_names[int(s[1])]

                rc = {"A": "T", "C": "G", "G": "C", "T": "A"}
                change = s[2].strip()
                if change not in types:
                    change = ">".join(rc[n] for n in change.split(">"))

                field = f"substitution_type_{change}"
                value = int(s[3].strip())
                if field not in bcftools_stats[s_name]:
                    bcftools_stats[s_name][field] = 0
                bcftools_stats[s_name][field] += value

            # Indel length distributions
            if s[0] == "IDD" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                length = int(s[2].strip())
                count = int(s[3].strip())
                bcftools_stats_indels[s_name][length] = count

            # Per-sample counts
            if s[0] == "PSC" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                fields = ["variations_hom", "variations_het"]
                for i, field in enumerate(fields):
                    bcftools_stats[s_name][field] = int(s[i + 4].strip())

            # Per-sample variant stats
            if s[0] == "PSC" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                sample = module.clean_s_name(s[2].strip(), f)
                bcftools_stats_sample_variants[s_name][sample] = dict()
                bcftools_stats_sample_variants[s_name][sample]["nSNPs"] = (
                    int(s[3].strip()) + int(s[4].strip()) + int(s[5].strip())
                )
                bcftools_stats_sample_variants[s_name][sample]["nIndels"] = int(s[8].strip())
                if len(s) >= 14:
                    nMissing = int(s[13].strip())
                else:
                    nMissing = 0
                nPresent = bcftools_stats[s_name]["number_of_records"] - nMissing
                bcftools_stats_sample_variants[s_name][sample]["nOther"] = (
                    nPresent
                    - bcftools_stats_sample_variants[s_name][sample]["nSNPs"]
                    - bcftools_stats_sample_variants[s_name][sample]["nIndels"]
                )

            # Per-sample ts/tv stats
            if s[0] == "PSC" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                sample = module.clean_s_name(s[2].strip(), f)
                bcftools_stats_sample_tstv[s_name][sample] = dict()
                if int(s[7].strip()) != 0:
                    bcftools_stats_sample_tstv[s_name][sample]["tstv"] = float(int(s[6].strip()) / int(s[7].strip()))
                else:
                    bcftools_stats_sample_tstv[s_name][sample]["tstv"] = 0

            # Per-sample singletons stats
            if s[0] == "PSC" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                sample = module.clean_s_name(s[2].strip(), f)
                bcftools_stats_sample_singletons[s_name][sample] = dict()
                bcftools_stats_sample_singletons[s_name][sample]["singletons"] = int(s[10].strip())
                bcftools_stats_sample_singletons[s_name][sample]["rest"] = (
                    int(s[3].strip()) + int(s[4].strip()) + int(s[5].strip()) - int(s[10].strip())
                )

            # Per-sample coverage stats
            if s[0] == "PSC" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                sample = module.clean_s_name(s[2].strip(), f)
                bcftools_stats_sample_depth[s_name][sample] = dict()
                bcftools_stats_sample_depth[s_name][sample]["depth"] = float(s[9].strip())

            # Depth plots
            if s[0] == "DP" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                bin_name = s[2].strip()
                percent_sites = float(s[-1].strip())
                bcftools_stats_depth_data[s_name][bin_name] = percent_sites

            # Variant Qualities
            if s[0] == "QUAL" and len(s_names) > 0:
                s_name = s_names[int(s[1])]
                try:
                    quality = float(s[2].strip())
                except ValueError:
                    quality = 0
                bcftools_stats_vqc_snp[s_name][quality] = int(s[3].strip())
                bcftools_stats_vqc_transi[s_name][quality] = int(s[4].strip())
                bcftools_stats_vqc_transv[s_name][quality] = int(s[5].strip())
                bcftools_stats_vqc_indels[s_name][quality] = int(s[6].strip())

    # Remove empty samples
    bcftools_stats = {k: v for k, v in bcftools_stats.items() if len(v) > 0}
    bcftools_stats_indels = {k: v for k, v in bcftools_stats_indels.items() if len(v) > 0}
    bcftools_stats_sample_variants = {k: v for k, v in bcftools_stats_sample_variants.items() if len(v) > 0}
    bcftools_stats_sample_tstv = {k: v for k, v in bcftools_stats_sample_tstv.items() if len(v) > 0}
    bcftools_stats_sample_singletons = {k: v for k, v in bcftools_stats_sample_singletons.items() if len(v) > 0}
    bcftools_stats_sample_depth = {k: v for k, v in bcftools_stats_sample_depth.items() if len(v) > 0}
    bcftools_stats_vqc_snp = {k: v for k, v in bcftools_stats_vqc_snp.items() if len(v) > 0}
    bcftools_stats_vqc_transi = {k: v for k, v in bcftools_stats_vqc_transi.items() if len(v) > 0}
    bcftools_stats_vqc_transv = {k: v for k, v in bcftools_stats_vqc_transv.items() if len(v) > 0}
    bcftools_stats_vqc_indels = {k: v for k, v in bcftools_stats_vqc_indels.items() if len(v) > 0}
    bcftools_stats_depth_data = {k: v for k, v in bcftools_stats_depth_data.items() if len(v) > 0}

    # Filter to strip out ignored sample names
    bcftools_stats = module.ignore_samples(bcftools_stats)
    bcftools_stats_indels = module.ignore_samples(bcftools_stats_indels)
    bcftools_stats_vqc_snp = module.ignore_samples(bcftools_stats_vqc_snp)
    bcftools_stats_sample_variants = module.ignore_samples(bcftools_stats_sample_variants)
    bcftools_stats_sample_tstv = module.ignore_samples(bcftools_stats_sample_tstv)
    bcftools_stats_sample_singletons = module.ignore_samples(bcftools_stats_sample_singletons)
    bcftools_stats_sample_depth = module.ignore_samples(bcftools_stats_sample_depth)
    bcftools_stats_vqc_transi = module.ignore_samples(bcftools_stats_vqc_transi)
    bcftools_stats_vqc_transv = module.ignore_samples(bcftools_stats_vqc_transv)
    bcftools_stats_vqc_indels = module.ignore_samples(bcftools_stats_vqc_indels)
    bcftools_stats_depth_data = module.ignore_samples(bcftools_stats_depth_data)
    if len(bcftools_stats) == 0:
        return 0

    # Write parsed report data to a file
    module.write_data_file(bcftools_stats, "multiqc_bcftools_stats")

    # Stats Table
    stats_headers = bcftools_stats_genstats_headers()
    if getattr(config, "bcftools", {}).get("write_general_stats", True):
        module.general_stats_addcols(bcftools_stats, stats_headers, "Stats")
    if getattr(config, "bcftools", {}).get("write_separate_table", False):
        module.add_section(
            name="Bcftools Stats",
            anchor="bcftools-stats_stats",
            plot=table.plot(
                bcftools_stats,
                stats_headers,
                {
                    "namespace": "Stats",
                    "id": "bcftools-stats-table",
                },
            ),
        )

    # Make bargraph plot of substitution types
    keys = {}
    for t in types:
        keys[f"substitution_type_{t}"] = {"name": t}
    module.add_section(
        name="Variant Substitution Types",
        anchor="bcftools-stats_variant_sub_types",
        plot=bargraph.plot(
            bcftools_stats,
            keys,
            BarPlotConfig(
                id="bcftools-stats-subtypes",
                title="Bcftools Stats: Substitutions",
                ylab="# Substitutions",
                cpswitch_counts_label="Number of Substitutions",
            ),
        ),
    )

    # Make histograms of variant quality
    if len(bcftools_stats_vqc_snp) > 0:
        module.add_section(
            name="Variant Quality",
            anchor="bcftools-stats_variant_quality_plot",
            plot=linegraph.plot(
                [
                    bcftools_stats_vqc_snp,
                    bcftools_stats_vqc_transi,
                    bcftools_stats_vqc_transv,
                    bcftools_stats_vqc_indels,
                ],
                LinePlotConfig(
                    id="bcftools_stats_vqc",
                    title="Bcftools Stats: Variant Quality Count",
                    ylab="Count",
                    xlab="Quality",
                    ymin=0,
                    smooth_points=600,
                    x_decimals=0,
                    data_labels=[
                        "Count SNP",
                        "Count Transitions",
                        "Count Transversions",
                        "Count Indels",
                    ],
                ),
            ),
        )

    # Make line graph of indel lengths
    if len(bcftools_stats_indels) > 0:
        module.add_section(
            name="Indel Distribution",
            anchor="bcftools-stats_indel_plot",
            plot=linegraph.plot(
                bcftools_stats_indels,
                LinePlotConfig(
                    id="bcftools_stats_indel-lengths",
                    title="Bcftools Stats: Indel Distribution",
                    ylab="Count",
                    xlab="InDel Length (bp)",
                    xsuffix=" bp",
                    ymin=0,
                ),
            ),
        )
    # Make line graph of variants per depth
    if len(bcftools_stats_depth_data) > 0:
        # Get shared list of bins and order them numerically
        all_bins = []
        for sname, val_by_bin in bcftools_stats_depth_data.items():
            all_bins.extend(list(val_by_bin.keys()))
        all_bins = sorted(all_bins, key=lambda x: int(re.sub(r"\D", "", x)))
        # Order bins in samples and fill missing bins:
        sorted_data = {
            sname: {b: bcftools_stats_depth_data[sname].get(b, 0) for b in all_bins}
            for sname in bcftools_stats_depth_data
        }
        module.add_section(
            name="Variant depths",
            anchor="bcftools-stats_depth_plot",
            description="Read depth support distribution for called variants",
            plot=linegraph.plot(
                sorted_data,
                LinePlotConfig(
                    id="bcftools_stats_variant_depths",
                    title="Bcftools Stats: Variant depths",
                    ylab="Fraction of sites (%)",
                    xlab="Variant depth",
                    ysuffix="%",
                    ymin=0,
                    ymax=100,
                    categories=True,
                    tt_decimals=1,
                ),
            ),
        )

    # Make bargraph plot of missing sites
    if len(bcftools_stats_sample_variants) > 0:
        module.add_section(
            name="Sites per sample",
            anchor="bcftools-stats_sites_per_sample",
            plot=bargraph.plot(
                list(bcftools_stats_sample_variants.values()),
                ["nSNPs", "nIndels", "nOther"],
                BarPlotConfig(
                    id="bcftools-stats-sites",
                    title="Bcftools Stats: Sites per sample",
                    ylab="# Sites",
                    cpswitch_counts_label="Number of sites",
                    data_labels=list(bcftools_stats_sample_variants),
                ),
            ),
        )

    # Make bargraph plot of ts/tv stats
    if len(bcftools_stats_sample_tstv) > 0:
        module.add_section(
            name="Ts/Tv",
            anchor="bcftools-stats_ts_tv",
            plot=bargraph.plot(
                list(bcftools_stats_sample_tstv.values()),
                ["tstv"],
                BarPlotConfig(
                    id="bcftools-stats-tstv",
                    title="Bcftools Stats: Ts/Tv",
                    ylab="Ts/Tv",
                    cpswitch_counts_label="Ts/Tv",
                    cpswitch=False,
                    data_labels=list(bcftools_stats_sample_tstv),
                ),
            ),
        )

    # Make bargraph plot of singletons stats
    if len(bcftools_stats_sample_singletons) > 0:
        module.add_section(
            name="Number of Singletons",
            anchor="bcftools-stats_singletones",
            plot=bargraph.plot(
                list(bcftools_stats_sample_singletons.values()),
                ["singletons", "rest"],
                BarPlotConfig(
                    id="bcftools-stats-singletons",
                    title="Bcftools Stats: Singletons",
                    ylab="# Singletons",
                    cpswitch_counts_label="Singletons",
                    cpswitch_c_active=False,
                    data_labels=list(bcftools_stats_sample_singletons),
                ),
            ),
        )

    # Make bargraph plot of sequencing depth stats
    if len(bcftools_stats_sample_depth) > 0:
        module.add_section(
            name="Sequencing depth",
            anchor="bcftools-stats_sequencing_depth",
            plot=bargraph.plot(
                list(bcftools_stats_sample_depth.values()),
                ["depth"],
                BarPlotConfig(
                    id="bcftools-stats-sequencing-depth",
                    title="Bcftools Stats: Sequencing depth",
                    ylab="Sequencing depth",
                    ysuffix="X",
                    cpswitch_counts_label="Sequencing depth",
                    cpswitch=False,
                    data_labels=list(bcftools_stats_sample_depth),
                ),
            ),
        )

    # Return the number of logs that were found
    return len(bcftools_stats)


def bcftools_stats_genstats_headers():
    """Add key statistics to the General Stats table"""
    stats_headers = {
        "number_of_records": {
            "title": "Vars",
            "description": "Variations total",
            "min": 0,
            "format": "{:,.0f}",
        },
        "variations_hom": {
            "title": "Hom",
            "description": "Variations homozygous",
            "min": 0,
            "format": "{:,.0f}",
        },
        "variations_het": {
            "title": "Het",
            "description": "Variations heterozygous",
            "min": 0,
            "format": "{:,.0f}",
        },
        "number_of_SNPs": {
            "title": "SNP",
            "description": "Variation SNPs",
            "min": 0,
            "format": "{:,.0f}",
        },
        "number_of_indels": {
            "title": "Indel",
            "description": "Variation Insertions/Deletions",
            "min": 0,
            "format": "{:,.0f}",
        },
        "tstv": {
            "title": "Ts/Tv",
            "description": "Variant SNP transition / transversion ratio",
            "min": 0,
            "format": "{:,.2f}",
        },
        "number_of_MNPs": {
            "title": "MNP",
            "description": "Variation multinucleotide polymorphisms",
            "min": 0,
            "format": "{:,.0f}",
            "hidden": True,
        },
        "number_of_multiallelic_sites": {
            "title": "Multiallelic",
            "description": "Variation sites with multiple alleles",
            "min": 0,
            "format": "{:,.0f}",
            "hidden": True,
        },
        "number_of_multiallelic_SNP_sites": {
            "title": "Multiallelic SNP",
            "description": "Variation sites with multiple SNPs",
            "min": 0,
            "format": "{:,.0f}",
            "hidden": True,
        },
    }
    return stats_headers
