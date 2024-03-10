""" MultiQC submodule to parse output from Bcftools stats """

import logging
import re

from multiqc import config
from multiqc.plots import bargraph, linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"# This file was produced by bcftools stats \(([\d\.]+)"
HTSLIB_REGEX = r"\+htslib-([\d\.]+)"


class StatsReportMixin:
    """Mixin loaded by the bcftools MultiqcModule class"""

    def parse_bcftools_stats(self):
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

        self.bcftools_stats = dict()
        self.bcftools_stats_indels = dict()
        self.bcftools_stats_vqc_snp = dict()
        self.bcftools_stats_sample_variants = dict()
        self.bcftools_stats_sample_tstv = dict()
        self.bcftools_stats_sample_singletons = dict()
        self.bcftools_stats_sample_depth = dict()
        self.bcftools_stats_vqc_transi = dict()
        self.bcftools_stats_vqc_transv = dict()
        self.bcftools_stats_vqc_indels = dict()
        self.bcftools_stats_depth_data = dict()
        for f in self.find_log_files("bcftools/stats"):
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
                    self.add_software_version(bcftools_version, f["s_name"])

                    # Look for HTSlib version
                    htslib_version_match = re.search(HTSLIB_REGEX, line)
                    if htslib_version_match is None:
                        continue

                    # Add HTSlib version if different from BCFtools version
                    htslib_version = htslib_version_match.group(1)
                    if htslib_version != bcftools_version:
                        self.add_software_version(htslib_version, f["s_name"], "HTSlib")

                s = line.split("\t")
                # Get the sample names - one per 'set'
                if s[0] == "ID":
                    s_name = self.clean_s_name(s[2], f)
                    s_names.append(s_name)
                    if s_name in self.bcftools_stats:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name, section="stats")
                    self.bcftools_stats[s_name] = dict()
                    self.bcftools_stats_indels[s_name] = dict()
                    self.bcftools_stats_sample_variants[s_name] = dict()
                    self.bcftools_stats_sample_tstv[s_name] = dict()
                    self.bcftools_stats_sample_singletons[s_name] = dict()
                    self.bcftools_stats_sample_depth[s_name] = dict()
                    self.bcftools_stats_vqc_snp[s_name] = dict()
                    self.bcftools_stats_vqc_transi[s_name] = dict()
                    self.bcftools_stats_vqc_transv[s_name] = dict()
                    self.bcftools_stats_vqc_indels[s_name] = dict()
                    self.bcftools_stats_depth_data[s_name] = {}
                    self.bcftools_stats_indels[s_name][0] = None  # Avoid joining line across missing 0

                # Parse key stats
                if s[0] == "SN" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    field = s[2].strip()[:-1]
                    field = field.replace(" ", "_")
                    value = int(s[3].strip())
                    self.bcftools_stats[s_name][field] = value

                # Parse transitions/transversions stats
                if s[0] == "TSTV" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    fields = ["ts", "tv", "tstv", "ts_1st_ALT", "tv_1st_ALT", "tstv_1st_ALT"]
                    for i, field in enumerate(fields):
                        value = s[i + 2].strip()
                        if "tstv" in field:
                            value = float(value)
                        else:
                            value = int(value)
                        self.bcftools_stats[s_name][field] = value

                # Parse substitution types
                if s[0] == "ST" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]

                    rc = {"A": "T", "C": "G", "G": "C", "T": "A"}
                    change = s[2].strip()
                    if change not in types:
                        change = ">".join(rc[n] for n in change.split(">"))

                    field = f"substitution_type_{change}"
                    value = int(s[3].strip())
                    if field not in self.bcftools_stats[s_name]:
                        self.bcftools_stats[s_name][field] = 0
                    self.bcftools_stats[s_name][field] += value

                # Indel length distributions
                if s[0] == "IDD" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    length = int(s[2].strip())
                    count = int(s[3].strip())
                    self.bcftools_stats_indels[s_name][length] = count

                # Per-sample counts
                if s[0] == "PSC" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    fields = ["variations_hom", "variations_het"]
                    for i, field in enumerate(fields):
                        self.bcftools_stats[s_name][field] = int(s[i + 4].strip())

                # Per-sample variant stats
                if s[0] == "PSC" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    sample = self.clean_s_name(s[2].strip(), f)
                    self.bcftools_stats_sample_variants[s_name][sample] = dict()
                    self.bcftools_stats_sample_variants[s_name][sample]["nSNPs"] = (
                        int(s[3].strip()) + int(s[4].strip()) + int(s[5].strip())
                    )
                    self.bcftools_stats_sample_variants[s_name][sample]["nIndels"] = int(s[8].strip())
                    if len(s) >= 14:
                        nMissing = int(s[13].strip())
                    else:
                        nMissing = 0
                    nPresent = self.bcftools_stats[s_name]["number_of_records"] - nMissing
                    self.bcftools_stats_sample_variants[s_name][sample]["nOther"] = (
                        nPresent
                        - self.bcftools_stats_sample_variants[s_name][sample]["nSNPs"]
                        - self.bcftools_stats_sample_variants[s_name][sample]["nIndels"]
                    )

                # Per-sample ts/tv stats
                if s[0] == "PSC" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    sample = self.clean_s_name(s[2].strip(), f)
                    self.bcftools_stats_sample_tstv[s_name][sample] = dict()
                    if int(s[7].strip()) != 0:
                        self.bcftools_stats_sample_tstv[s_name][sample]["tstv"] = float(
                            int(s[6].strip()) / int(s[7].strip())
                        )
                    else:
                        self.bcftools_stats_sample_tstv[s_name][sample]["tstv"] = 0

                # Per-sample singletons stats
                if s[0] == "PSC" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    sample = self.clean_s_name(s[2].strip(), f)
                    self.bcftools_stats_sample_singletons[s_name][sample] = dict()
                    self.bcftools_stats_sample_singletons[s_name][sample]["singletons"] = int(s[10].strip())
                    self.bcftools_stats_sample_singletons[s_name][sample]["rest"] = (
                        int(s[3].strip()) + int(s[4].strip()) + int(s[5].strip()) - int(s[10].strip())
                    )

                # Per-sample coverage stats
                if s[0] == "PSC" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    sample = self.clean_s_name(s[2].strip(), f)
                    self.bcftools_stats_sample_depth[s_name][sample] = dict()
                    self.bcftools_stats_sample_depth[s_name][sample]["depth"] = float(s[9].strip())

                # Depth plots
                if s[0] == "DP" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    bin_name = s[2].strip()
                    percent_sites = float(s[-1].strip())
                    self.bcftools_stats_depth_data[s_name][bin_name] = percent_sites

                # Variant Qualities
                if s[0] == "QUAL" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    try:
                        quality = float(s[2].strip())
                    except ValueError:
                        quality = 0
                    self.bcftools_stats_vqc_snp[s_name][quality] = int(s[3].strip())
                    self.bcftools_stats_vqc_transi[s_name][quality] = int(s[4].strip())
                    self.bcftools_stats_vqc_transv[s_name][quality] = int(s[5].strip())
                    self.bcftools_stats_vqc_indels[s_name][quality] = int(s[6].strip())

        # Remove empty samples
        self.bcftools_stats = {k: v for k, v in self.bcftools_stats.items() if len(v) > 0}
        self.bcftools_stats_indels = {k: v for k, v in self.bcftools_stats_indels.items() if len(v) > 0}
        self.bcftools_stats_sample_variants = {
            k: v for k, v in self.bcftools_stats_sample_variants.items() if len(v) > 0
        }
        self.bcftools_stats_sample_tstv = {k: v for k, v in self.bcftools_stats_sample_tstv.items() if len(v) > 0}
        self.bcftools_stats_sample_singletons = {
            k: v for k, v in self.bcftools_stats_sample_singletons.items() if len(v) > 0
        }
        self.bcftools_stats_sample_depth = {k: v for k, v in self.bcftools_stats_sample_depth.items() if len(v) > 0}
        self.bcftools_stats_vqc_snp = {k: v for k, v in self.bcftools_stats_vqc_snp.items() if len(v) > 0}
        self.bcftools_stats_vqc_transi = {k: v for k, v in self.bcftools_stats_vqc_transi.items() if len(v) > 0}
        self.bcftools_stats_vqc_transv = {k: v for k, v in self.bcftools_stats_vqc_transv.items() if len(v) > 0}
        self.bcftools_stats_vqc_indels = {k: v for k, v in self.bcftools_stats_vqc_indels.items() if len(v) > 0}
        self.bcftools_stats_depth_data = {k: v for k, v in self.bcftools_stats_depth_data.items() if len(v) > 0}

        # Filter to strip out ignored sample names
        self.bcftools_stats = self.ignore_samples(self.bcftools_stats)
        self.bcftools_stats_indels = self.ignore_samples(self.bcftools_stats_indels)
        self.bcftools_stats_vqc_snp = self.ignore_samples(self.bcftools_stats_vqc_snp)
        self.bcftools_stats_sample_variants = self.ignore_samples(self.bcftools_stats_sample_variants)
        self.bcftools_stats_sample_tstv = self.ignore_samples(self.bcftools_stats_sample_tstv)
        self.bcftools_stats_sample_singletons = self.ignore_samples(self.bcftools_stats_sample_singletons)
        self.bcftools_stats_sample_depth = self.ignore_samples(self.bcftools_stats_sample_depth)
        self.bcftools_stats_vqc_transi = self.ignore_samples(self.bcftools_stats_vqc_transi)
        self.bcftools_stats_vqc_transv = self.ignore_samples(self.bcftools_stats_vqc_transv)
        self.bcftools_stats_vqc_indels = self.ignore_samples(self.bcftools_stats_vqc_indels)
        self.bcftools_stats_depth_data = self.ignore_samples(self.bcftools_stats_depth_data)
        if len(self.bcftools_stats) == 0:
            return 0

        # Write parsed report data to a file
        self.write_data_file(self.bcftools_stats, "multiqc_bcftools_stats")

        # Stats Table
        stats_headers = self.bcftools_stats_genstats_headers()
        if getattr(config, "bcftools", {}).get("write_general_stats", True):
            self.general_stats_addcols(self.bcftools_stats, stats_headers, "Stats")
        if getattr(config, "bcftools", {}).get("write_separate_table", False):
            self.add_section(
                name="Bcftools Stats",
                anchor="bcftools-stats_stats",
                plot=table.plot(
                    self.bcftools_stats,
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
        pconfig = {
            "id": "bcftools-stats-subtypes",
            "title": "Bcftools Stats: Substitutions",
            "ylab": "# Substitutions",
            "cpswitch_counts_label": "Number of Substitutions",
        }
        self.add_section(
            name="Variant Substitution Types",
            anchor="bcftools-stats_variant_sub_types",
            plot=bargraph.plot(self.bcftools_stats, keys, pconfig),
        )

        # Make histograms of variant quality
        if len(self.bcftools_stats_vqc_snp) > 0:
            pconfig = {
                "id": "bcftools_stats_vqc",
                "title": "Bcftools Stats: Variant Quality Count",
                "ylab": "Count",
                "xlab": "Quality",
                "ymin": 0,
                "smooth_points": 600,
                "x_decimals": 0,
                "data_labels": [
                    "Count SNP",
                    "Count Transitions",
                    "Count Transversions",
                    "Count Indels",
                ],
            }
            self.add_section(
                name="Variant Quality",
                anchor="bcftools-stats_variant_quality_plot",
                plot=linegraph.plot(
                    [
                        self.bcftools_stats_vqc_snp,
                        self.bcftools_stats_vqc_transi,
                        self.bcftools_stats_vqc_transv,
                        self.bcftools_stats_vqc_indels,
                    ],
                    pconfig,
                ),
            )

        # Make line graph of indel lengths
        if len(self.bcftools_stats_indels) > 0:
            pconfig = {
                "id": "bcftools_stats_indel-lengths",
                "title": "Bcftools Stats: Indel Distribution",
                "ylab": "Count",
                "xlab": "InDel Length (bp)",
                "xsuffix": " bp",
                "ymin": 0,
            }
            self.add_section(
                name="Indel Distribution",
                anchor="bcftools-stats_indel_plot",
                plot=linegraph.plot(self.bcftools_stats_indels, pconfig),
            )
        # Make line graph of variants per depth
        if len(self.bcftools_stats_depth_data) > 0:
            # Get shared list of bins and order them numerically
            all_bins = []
            for sname, val_by_bin in self.bcftools_stats_depth_data.items():
                all_bins.extend(list(val_by_bin.keys()))
            all_bins = sorted(all_bins, key=lambda x: int(re.sub(r"\D", "", x)))
            # Order bins in samples and fill missing bins:
            sorted_data = {
                sname: {b: self.bcftools_stats_depth_data[sname].get(b, 0) for b in all_bins}
                for sname in self.bcftools_stats_depth_data
            }
            pconfig = {
                "id": "bcftools_stats_variant_depths",
                "title": "Bcftools Stats: Variant depths",
                "ylab": "Fraction of sites (%)",
                "xlab": "Variant depth",
                "ysuffix": "%",
                "ymin": 0,
                "ymax": 100,
                "categories": True,
                "tt_decimals": 1,
            }
            self.add_section(
                name="Variant depths",
                anchor="bcftools-stats_depth_plot",
                description="Read depth support distribution for called variants",
                plot=linegraph.plot(sorted_data, pconfig),
            )

        # Make bargraph plot of missing sites
        if len(self.bcftools_stats_sample_variants) > 0:
            pconfig = {
                "id": "bcftools-stats-sites",
                "title": "Bcftools Stats: Sites per sample",
                "ylab": "# Sites",
                "cpswitch_counts_label": "Number of sites",
                "data_labels": list(self.bcftools_stats_sample_variants),
            }
            self.add_section(
                name="Sites per sample",
                anchor="bcftools-stats_sites_per_sample",
                plot=bargraph.plot(
                    list(self.bcftools_stats_sample_variants.values()), ["nSNPs", "nIndels", "nOther"], pconfig
                ),
            )

        # Make bargraph plot of ts/tv stats
        if len(self.bcftools_stats_sample_tstv) > 0:
            pconfig = {
                "id": "bcftools-stats-tstv",
                "title": "Bcftools Stats: Ts/Tv",
                "ylab": "Ts/Tv",
                "cpswitch_counts_label": "Ts/Tv",
                "cpswitch": False,
                "data_labels": list(self.bcftools_stats_sample_tstv),
            }
            self.add_section(
                name="Ts/Tv",
                anchor="bcftools-stats_ts_tv",
                plot=bargraph.plot(list(self.bcftools_stats_sample_tstv.values()), ["tstv"], pconfig),
            )

        # Make bargraph plot of singletons stats
        if len(self.bcftools_stats_sample_singletons) > 0:
            pconfig = {
                "id": "bcftools-stats-singletons",
                "title": "Bcftools Stats: Singletons",
                "ylab": "# Singletons",
                "cpswitch_counts_label": "Singletons",
                "cpswitch_c_active": False,
                "data_labels": list(self.bcftools_stats_sample_singletons),
            }
            self.add_section(
                name="Number of Singletons",
                anchor="bcftools-stats_singletones",
                plot=bargraph.plot(
                    list(self.bcftools_stats_sample_singletons.values()), ["singletons", "rest"], pconfig
                ),
            )

        # Make bargraph plot of sequencing depth stats
        if len(self.bcftools_stats_sample_depth) > 0:
            pconfig = {
                "id": "bcftools-stats-sequencing-depth",
                "title": "Bcftools Stats: Sequencing depth",
                "ylab": "Sequencing depth",
                "ysuffix": "X",
                "cpswitch_counts_label": "Sequencing depth",
                "cpswitch": False,
                "data_labels": list(self.bcftools_stats_sample_depth),
            }
            self.add_section(
                name="Sequencing depth",
                anchor="bcftools-stats_sequencing_depth",
                plot=bargraph.plot(list(self.bcftools_stats_sample_depth.values()), ["depth"], pconfig),
            )

        # Return the number of logs that were found
        return len(self.bcftools_stats)

    @staticmethod
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
