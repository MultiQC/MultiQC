import logging
import os.path
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    **Note** - because this module shares sample identifiers across multiple files,
    the `--fn_as_s_name` / `config.use_filename_as_sample_name` functionality has been disabled and has no effect.

    The MultiQC module is supported since HiC-Pro v2.11.0.
    """

    """HiC-Pro module, parses log and stats files saved by HiC-Pro."""

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HiC-Pro",
            anchor="hicpro",
            href="https://github.com/nservant/HiC-Pro",
            info="Pipeline for Hi-C data processing",
            doi="10.1186/s13059-015-0831-x",
        )

        # Find and load any HiC-Pro summary reports
        self.hicpro_data: Dict = dict()
        for k in ["mmapstat", "mpairstat", "mergestat", "mRSstat", "assplit"]:
            for f in self.find_log_files(f"hicpro/{k}"):
                self.parse_hicpro_stats(f, k)

        # Update current statistics
        for s_name in self.hicpro_data:
            data = self.hicpro_data[s_name]
            try:
                if "valid_interaction" in data:
                    if "valid_interaction_rmdup" in data:
                        data["duplicates"] = data["valid_interaction"] - data["valid_interaction_rmdup"]
                        if float(data["valid_interaction"]) > 0:
                            data["percent_duplicates"] = (
                                float(data["duplicates"]) / float(data["valid_interaction"]) * 100.0
                            )
                    if "valid_pairs_on_target" in data:
                        data["valid_pairs_off_target"] = int(data["valid_interaction"]) - int(
                            data["valid_pairs_on_target"]
                        )
                    if "total_R1" in data and float(data["total_R1"]) > 0:
                        data["percent_valid"] = float(data["valid_interaction"]) / float(data["total_R1"]) * 100.0

                if "Reported_pairs" in data:
                    data["paired_reads"] = int(data["Reported_pairs"])
                    if "total_R1" in data and float(data["total_R1"]) > 0:
                        data["percent_paired_reads"] = float(data["Reported_pairs"]) / float(data["total_R1"]) * 100.0

                if "mapped_R1" in data:
                    data["Failed_To_Align_Read_R1"] = int(data["total_R1"]) - int(data["mapped_R1"])
                    if "total_R1" in data and float(data["total_R1"]) > 0:
                        data["percent_mapped_R1"] = float(data["mapped_R1"]) / float(data["total_R1"]) * 100.0

                if "mapped_R2" in data:
                    data["Failed_To_Align_Read_R2"] = int(data["total_R2"]) - int(data["mapped_R2"])
                    if "total_R2" in data and float(data["total_R2"]) > 0:
                        data["percent_mapped_R2"] = float(data["mapped_R2"]) / float(data["total_R2"]) * 100.0

            except KeyError as e:
                log.error(f"Missing expected key {e} in sample '{s_name}'")
            except (ValueError, TypeError):
                pass

        # Filter to strip out ignored sample names
        self.hicpro_data = self.ignore_samples(self.hicpro_data)

        if len(self.hicpro_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.hicpro_data)} HiC-Pro reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed data to a file
        self.write_data_file(self.hicpro_data, "multiqc_hicpro")

        # Basic Stats Table
        self.hicpro_stats_table()

        # Report sections
        self.add_section(
            name="Read Mapping",
            anchor="hicpro-mapping",
            description="Alignment of reads in single-end mode.",
            helptext="""
            HiC-Pro uses a two steps mapping strategy. End-to-end mapping is first
            performed (*Full reads Alignments*). Unmapped reads are then trimmed
            at their ligation site, and the 5\' end is re-aligned on the reference
            genome (*Trimmed Reads Alignments*)""",
            plot=self.hicpro_mapping_chart(),
        )

        self.add_section(
            name="Read Pairing",
            anchor="hicpro-pairing",
            description="Pairing of single-end mapping.",
            helptext="""
            HiC-Pro combines aligned reads as pairs. Singleton, low quality pairs or
            pairs involving multi-mapped reads are usually discarded. Note that the filtering
            at pairing level can change accrding to the parameters used.""",
            plot=self.hicpro_pairing_chart(),
        )

        self.add_section(
            name="Read Pair Filtering",
            anchor="hicpro-filtering",
            description="Filtering of read pairs.",
            helptext="""
            Aligned read pairs are filtered to select valid 3C products from two
            different restriction fragments. Read pairs coming from the same fragments,
            such as self-ligation or unligated (danging-end) fragments, are discarded.
            Ligation products involving neighboring restriction fragment (religation) are also discarded.
            Filtered pairs correspond to ligation products discarded based on the range of insert/fragment sizes defined
            in the analysis. Finally, as the ligation should be a random process, valid read pairs from all
            orientations (R=Reverse, F=forward) are expected to be observed in the same proportion.""",
            plot=self.hicpro_filtering_chart(),
        )

        self.add_section(
            name="Contact Statistics",
            anchor="hicpro-rmdup",
            description="Contacts statistics after duplicates removal.",
            helptext="""
            Description of contact frequency after duplicates removal. Intra-chromosomal
            (*cis*) interaction are expected to be more frequent than inter-chromosomal contacts (*trans*)""",
            plot=self.hicpro_contact_chart(),
        )

        allele_plot = self.hicpro_as_chart()
        if allele_plot:
            self.add_section(
                name="Allele-specific analysis",
                anchor="hicpro-asan",
                description="Assignment of valid interactions (afer duplicates removal) to parental alleles.",
                helptext="""
                Description of allele-specific contacts. Valid interactions (0-1 and 1-1, resp. 0-2, 2-2)
                are used to build the genome1 (resp. genome2) parental maps.""",
                plot=allele_plot,
            )

        capture_plot = self.hicpro_capture_chart()
        if capture_plot:
            self.add_section(
                name="Capture analysis",
                anchor="hicpro-cap",
                description="Selection of interactions overlapping the targeted region(s).",
                helptext="""
                Description of capture efficiency. Valid interactions with either two (capture-capture) or
                one (capture-reporter) interactors overlapping with the target(s) are reported.""",
                plot=capture_plot,
            )

    def parse_hicpro_stats(self, f, rsection):
        """Parse a HiC-Pro stat file"""
        s_name = self.clean_s_name(os.path.basename(f["root"]), root=os.path.dirname(f["root"]))
        if s_name not in self.hicpro_data.keys():
            self.hicpro_data[s_name] = {}

        self.add_data_source(f, s_name, section=rsection)
        for line in f["f"].splitlines():
            if not line.startswith("#"):
                s = line.split("\t")
                if s[0] in self.hicpro_data[s_name]:
                    log.debug(
                        f"Duplicated key {s[0]} found in {f['root'] + '/' + f['fn']} for sample {s_name}, overwriting"
                    )
                # Try to convert the extracted value to a number and store it in hicpro_data.
                # try-block is used to prevent program crash, because there is no
                # guarantee that the value (s[1]) can be always converted to integer.
                try:
                    self.hicpro_data[s_name][s[0]] = int(s[1])
                except ValueError:
                    # Convert to float (also works for inf and scientific (exponential) notation).
                    try:
                        self.hicpro_data[s_name][s[0]] = float(s[1])
                    # Otherwise just store the value as is.
                    except ValueError:
                        log.error(f"Could not parse value for key '{s[0]}' as a number in sample '{s_name}': '{s[1]}'")
                        self.hicpro_data[s_name][s[0]] = s[1]

    def hicpro_stats_table(self):
        """Add HiC-Pro stats to the general stats table"""
        headers = {
            "percent_duplicates": {
                "title": "% Duplicates",
                "description": "Percent of duplicated valid pairs (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlOrRd",
                "hidden": True,
            },
            "valid_interaction_rmdup": {
                "title": f"{config.read_count_prefix} Valid Pairs Unique",
                "description": f"Number of valid pairs after duplicates removal ({config.read_count_desc})",
                "min": 0,
                "scale": "RdYlBu",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
            "percent_valid": {
                "title": "% Valid Pairs",
                "description": "Percentage of valid pairs over reported ones (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": True,
            },
            "valid_interaction": {
                "title": f"{config.read_count_prefix} Valid Pairs",
                "description": f"Number of valid pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "RdYlBu",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "hidden": True,
            },
            "percent_paired_reads": {
                "title": "% Reported",
                "description": "Percentage of paired reads (%) passing the mapping filters",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": True,
            },
            "paired_reads": {
                "title": "Reported Read Pairs",
                "description": "Total number of read pairs ({}) passing the mapping filters".format(
                    config.read_count_desc
                ),
                "min": 0,
                "scale": "RdYlBu",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
            "percent_mapped_R2": {
                "title": "% Aligned [R2]",
                "description": "Percentage of aligned reads [R2] (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": True,
            },
            "mapped_R2": {
                "title": "Aligned [R2]",
                "description": f"Total number of aligned reads [R2] ({config.read_count_desc})",
                "min": 0,
                "scale": "RdYlBu",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "hidden": True,
            },
            "percent_mapped_R1": {
                "title": "% Aligned [R1]",
                "description": "Percentage of aligned reads [R1] (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": True,
            },
            "mapped_R1": {
                "title": "Aligned [R1]",
                "description": f"Total number of aligned reads [R1] ({config.read_count_desc})",
                "min": 0,
                "scale": "RdYlBu",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "hidden": True,
            },
            "total_R1": {
                "title": "Total",
                "description": "Total Number of Read Pairs",
                "min": 0,
                "scale": "RdYlBu",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }

        self.general_stats_addcols(self.hicpro_data, headers)

    def hicpro_mapping_chart(self):
        """Generate the HiC-Pro Aligned reads plot"""

        # Specify the order of the different possible categories
        keys = {
            "Full_Alignments_Read": {"color": "#005ce6", "name": "Full reads Alignments"},
            "Trimmed_Alignments_Read": {"color": "#3385ff", "name": "Trimmed reads Alignments"},
            "Failed_To_Align_Read": {"color": "#a9a2a2", "name": "Failed To Align"},
        }

        data = [{}, {}]
        for s_name in self.hicpro_data:
            for r in [1, 2]:
                n_global = self.hicpro_data[s_name].get(f"global_R{r}")
                n_local = self.hicpro_data[s_name].get(f"local_R{r}")
                n_mapped = self.hicpro_data[s_name].get(f"mapped_R{r}")
                n_total = self.hicpro_data[s_name].get(f"total_R{r}")
                d = dict()
                if n_global is not None:
                    d["Full_Alignments_Read"] = n_global
                if n_local is not None:
                    d["Trimmed_Alignments_Read"] = n_local
                if n_mapped is not None and n_total is not None:
                    try:
                        d["Failed_To_Align_Read"] = int(n_total) - int(n_mapped)
                    except ValueError:
                        log.warning(
                            f"Could not parse mapped_R{r}={n_mapped} or total_R{r}={n_total} "
                            f"as an integer number for sample {s_name}"
                        )
                data[r - 1][f"{s_name} [R{r}]"] = d

        # Config for the plot
        config = {
            "id": "hicpro_mapping_stats_plot",
            "title": "HiC-Pro: Mapping Statistics",
            "ylab": "# Reads: Read 1",
            "data_labels": [
                {"name": "Read 1", "ylab": "# Reads: Read 1"},
                {"name": "Read 2", "ylab": "# Reads: Read 2"},
            ],
        }

        return bargraph.plot(data, [keys, keys], config)

    def hicpro_pairing_chart(self):
        """Generate Pairing chart"""

        # Specify the order of the different possible categories
        keys = {
            "Unique_paired_alignments": {"color": "#005ce6", "name": "Uniquely Aligned"},
            "Low_qual_pairs": {"color": "#b97b35", "name": "Low Quality"},
            "Pairs_with_singleton": {"color": "#ff9933", "name": "Singleton"},
            "Multiple_pairs_alignments": {"color": "#e67300", "name": "Multi Aligned"},
            "Unmapped_airs": {"color": "#a9a2a2", "name": "Failed To Align"},
        }

        # Config for the plot
        config = {
            "id": "hicpro_pairing_stats_plot",
            "title": "HiC-Pro: Pairing Statistics",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        if not any([k in self.hicpro_data[s_name] for s_name in self.hicpro_data for k in keys]):
            # No data found to build this plot
            return

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_filtering_chart(self):
        """Generate the HiC-Pro filtering plot"""

        # Specify the order of the different possible categories
        keys = {
            "Valid_interaction_pairs_FF": {"color": "#ccddff", "name": "Valid Pairs FF"},
            "Valid_interaction_pairs_RR": {"color": "#6699ff", "name": "Valid Pairs RR"},
            "Valid_interaction_pairs_RF": {"color": "#0055ff", "name": "Valid Pairs RF"},
            "Valid_interaction_pairs_FR": {"color": "#003399", "name": "Valid Pairs FR"},
            "Self_Cycle_pairs": {"color": "#ffad99", "name": "Same Fragment - Self-Circle"},
            "Dangling_end_pairs": {"color": "#ff5c33", "name": "Same Fragment - Dangling Ends"},
            "Religation_pairs": {"color": "#cc2900", "name": "Re-ligation"},
            "Filtered_pairs": {"color": "#661400", "name": "Filtered pairs"},
            "Dumped_pairs": {"color": "#330a00", "name": "Dumped pairs"},
        }

        # Config for the plot
        config = {
            "id": "hicpro_filtering_plot",
            "title": "HiC-Pro: Filtering Statistics",
            "ylab": "# Read Pairs",
            "cpswitch_counts_label": "Number of Read Pairs",
        }

        if not any([k in self.hicpro_data[s_name] for s_name in self.hicpro_data for k in keys]):
            # No data found to build this plot
            return

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_contact_chart(self):
        """Generate the HiC-Pro interaction plot"""

        # Specify the order of the different possible categories
        keys = {
            "cis_shortRange": {"color": "#0039e6", "name": "Unique: cis <= 20Kbp"},
            "cis_longRange": {"color": "#809fff", "name": "Unique: cis > 20Kbp"},
            "trans_interaction": {"color": "#009933", "name": "Unique: trans"},
            "duplicates": {"color": "#a9a2a2", "name": "Duplicate read pairs"},
        }

        # Config for the plot
        config = {
            "id": "hicpro_contact_plot",
            "title": "HiC-Pro: Contact Statistics",
            "ylab": "# Pairs",
            "cpswitch_counts_label": "Number of Pairs",
        }

        if not any([k in self.hicpro_data[s_name] for s_name in self.hicpro_data for k in keys]):
            # No data found to build this plot
            return

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_as_chart(self):
        """Generate Allele-specific plot"""

        keys = {
            "Valid_pairs_from_ref_genome_(1-1)": {"color": "#e6550d", "name": "Genome1 specific read pairs (1-1)"},
            "Valid_pairs_from_ref_genome_with_one_unassigned_mate_(0-1/1-0)": {
                "color": "#fdae6b",
                "name": "Genome1 with one unassigned mate (0-1/1-0)",
            },
            "Valid_pairs_from_alt_genome_(2-2)": {"color": "#756bb1", "name": "Genome2 specific read pairs (2-2)"},
            "Valid_pairs_from_alt_genome_with_one_unassigned_mate_(0-2/2-0)": {
                "color": "#bcbddc",
                "name": "Genome2 with one unassigned mate (0-2/2-0)",
            },
            "Valid_pairs_from_alt_and_ref_genome_(1-2/2-1)": {
                "color": "#a6611a",
                "name": "Trans homologous read pairs (1-2/2/1)",
            },
            "Valid_pairs_with_both_unassigned_mated_(0-0)": {"color": "#cccccc", "name": "Unassigned read pairs"},
            "Valid_pairs_with_at_least_one_conflicting_mate_(3-)": {
                "color": "#a9a2a2",
                "name": "Conflicting read pairs",
            },
        }

        # Config for the plot
        config = {
            "id": "hicpro_asan_plot",
            "title": "HiC-Pro: Allele-specific Statistics",
            "ylab": "# Pairs",
            "cpswitch_counts_label": "Number of Pairs",
        }

        if not any([k in self.hicpro_data[s_name] for s_name in self.hicpro_data for k in keys]):
            # No data found to build this plot
            return

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_capture_chart(self):
        """Generate Capture Hi-C plot"""

        keys = {
            "valid_pairs_on_target_cap_cap": {"color": "#0039e6", "name": "Capture-Capture interactions"},
            "valid_pairs_on_target_cap_rep": {"color": "#809fff", "name": "Capture-Reporter interactions"},
            "valid_pairs_off_target": {"color": "#cccccc", "name": "Off-target valid pairs"},
        }

        # Config for the plot
        config = {
            "id": "hicpro_cap_plot",
            "title": "HiC-Pro: Capture Statistics",
            "ylab": "# Pairs",
            "cpswitch_counts_label": "Number of Pairs",
        }

        if not any([k in self.hicpro_data[s_name] for s_name in self.hicpro_data for k in keys]):
            # No data found to build this plot
            return

        return bargraph.plot(self.hicpro_data, keys, config)
