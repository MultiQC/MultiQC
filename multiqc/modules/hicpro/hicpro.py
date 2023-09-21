## Nicolas Servant
## April 2018

""" MultiQC module to parse output from HiC-Pro """


import logging
import os.path
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """HiC-Pro module, parses log and stats files saved by HiC-Pro."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="HiC-Pro",
            anchor="hicpro",
            href="https://github.com/nservant/HiC-Pro",
            info="is an efficient and flexible pipeline for Hi-C data processing. The MultiQC module is supported since HiC-Pro v2.11.0.",
            doi="10.1186/s13059-015-0831-x",
        )

        # Find and load any HiC-Pro summary reports
        self.hicpro_data = dict()
        for k in ["mmapstat", "mpairstat", "mergestat", "mRSstat", "assplit"]:
            for f in self.find_log_files("hicpro/{}".format(k)):
                self.parse_hicpro_stats(f, k)

        # Update current statistics
        for s_name in self.hicpro_data:
            data = self.hicpro_data[s_name]
            try:
                data["duplicates"] = data["valid_interaction"] - data["valid_interaction_rmdup"]
                data["percent_duplicates"] = float(data["duplicates"]) / float(data["valid_interaction"]) * 100.0
                data["percent_mapped_R1"] = float(data["mapped_R1"]) / float(data["total_R1"]) * 100.0
                data["percent_mapped_R2"] = float(data["mapped_R2"]) / float(data["total_R2"]) * 100.0
                data["paired_reads"] = int(data["Reported_pairs"])
                data["percent_paired_reads"] = float(data["Reported_pairs"]) / float(data["total_R1"]) * 100.0
                data["percent_valid"] = float(data["valid_interaction"]) / float(data["total_R1"]) * 100.0
                data["Failed_To_Align_Read_R1"] = int(data["total_R1"]) - int(data["mapped_R1"])
                data["Failed_To_Align_Read_R2"] = int(data["total_R2"]) - int(data["mapped_R2"])
            except KeyError as e:
                log.error(f"Missing expected key {e} in sample '{s_name}'")

            try:
                data["valid_pairs_off_target"] = int(data["valid_interaction"]) - int(data["valid_pairs_on_target"])
            except KeyError:
                pass

        # Filter to strip out ignored sample names
        self.hicpro_data = self.ignore_samples(self.hicpro_data)

        if len(self.hicpro_data) == 0:
            raise UserWarning

        log.info("Found {} HiC-Pro reports".format(len(self.hicpro_data)))

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
                description="Selection of interactions overlaping the targeted region(s).",
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
        for l in f["f"].splitlines():
            if not l.startswith("#"):
                s = l.split("\t")
                if s[0] in self.hicpro_data[s_name]:
                    log.debug("Duplicated keys found! Overwriting: {}".format(s[0]))
                self.hicpro_data[s_name][s[0]] = int(s[1])

    def hicpro_stats_table(self):
        """Add HiC-Pro stats to the general stats table"""
        headers = OrderedDict()

        headers["percent_duplicates"] = {
            "title": "% Duplicates",
            "description": "Percent of duplicated valid pairs (%)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlOrRd",
            "hidden": True,
        }

        headers["valid_interaction_rmdup"] = {
            "title": "{} Valid Pairs Unique".format(config.read_count_prefix),
            "description": "Number of valid pairs after duplicates removal ({})".format(config.read_count_desc),
            "min": 0,
            "scale": "RdYlBu",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        }

        headers["percent_valid"] = {
            "title": "% Valid Pairs",
            "description": "Percentage of valid pairs over reported ones (%)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGn",
            "hidden": True,
        }

        headers["valid_interaction"] = {
            "title": "{} Valid Pairs".format(config.read_count_prefix),
            "description": "Number of valid pairs ({})".format(config.read_count_desc),
            "min": 0,
            "scale": "RdYlBu",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "hidden": True,
        }

        headers["percent_paired_reads"] = {
            "title": "% Reported",
            "description": "Percentage of paired reads (%) passing the mapping filters",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGn",
            "hidden": True,
        }

        headers["paired_reads"] = {
            "title": "Reported Read Pairs",
            "description": "Total number of read pairs ({}) passing the mapping filters".format(config.read_count_desc),
            "min": "0",
            "scale": "RdYlBu",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        }

        headers["percent_mapped_R2"] = {
            "title": "% Aligned [R2]",
            "description": "Percentage of aligned reads [R2] (%)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGn",
            "hidden": True,
        }

        headers["mapped_R2"] = {
            "title": "Aligned [R2]",
            "description": "Total number of aligned reads [R2] ({})".format(config.read_count_desc),
            "min": "0",
            "scale": "RdYlBu",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "hidden": True,
        }

        headers["percent_mapped_R1"] = {
            "title": "% Aligned [R1]",
            "description": "Percentage of aligned reads [R1] (%)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGn",
            "hidden": True,
        }

        headers["mapped_R1"] = {
            "title": "Aligned [R1]",
            "description": "Total number of aligned reads [R1] ({})".format(config.read_count_desc),
            "min": "0",
            "scale": "RdYlBu",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "hidden": True,
        }

        headers["total_R1"] = {
            "title": "Total",
            "description": "Total Number of Read Pairs",
            "min": "0",
            "scale": "RdYlBu",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        }

        self.general_stats_addcols(self.hicpro_data, headers)

    def hicpro_mapping_chart(self):
        """Generate the HiC-Pro Aligned reads plot"""

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys["Full_Alignments_Read"] = {"color": "#005ce6", "name": "Full reads Alignments"}
        keys["Trimmed_Alignments_Read"] = {"color": "#3385ff", "name": "Trimmed reads Alignments"}
        keys["Failed_To_Align_Read"] = {"color": "#a9a2a2", "name": "Failed To Align"}

        data = [{}, {}]
        for s_name in self.hicpro_data:
            for r in [1, 2]:
                try:
                    data[r - 1]["{} [R{}]".format(s_name, r)] = {
                        "Full_Alignments_Read": self.hicpro_data[s_name]["global_R{}".format(r)],
                        "Trimmed_Alignments_Read": self.hicpro_data[s_name]["local_R{}".format(r)],
                        "Failed_To_Align_Read": int(self.hicpro_data[s_name]["total_R{}".format(r)])
                        - int(self.hicpro_data[s_name]["mapped_R{}".format(r)]),
                    }
                except KeyError as e:
                    log.error(f"Missing expected plot key {e} in {s_name} Read {r}")

        # Config for the plot
        config = {
            "id": "hicpro_mapping_stats_plot",
            "title": "HiC-Pro: Mapping Statistics",
            "ylab": "# Reads",
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
        keys = OrderedDict()
        keys["Unique_paired_alignments"] = {"color": "#005ce6", "name": "Uniquely Aligned"}
        keys["Low_qual_pairs"] = {"color": "#b97b35", "name": "Low Quality"}
        keys["Pairs_with_singleton"] = {"color": "#ff9933", "name": "Singleton"}
        keys["Multiple_pairs_alignments"] = {"color": "#e67300", "name": "Multi Aligned"}
        keys["Unmapped_airs"] = {"color": "#a9a2a2", "name": "Failed To Align"}

        # Config for the plot
        config = {
            "id": "hicpro_pairing_stats_plot",
            "title": "HiC-Pro: Pairing Statistics",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_filtering_chart(self):
        """Generate the HiC-Pro filtering plot"""

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys["Valid_interaction_pairs_FF"] = {"color": "#ccddff", "name": "Valid Pairs FF"}
        keys["Valid_interaction_pairs_RR"] = {"color": "#6699ff", "name": "Valid Pairs RR"}
        keys["Valid_interaction_pairs_RF"] = {"color": "#0055ff", "name": "Valid Pairs RF"}
        keys["Valid_interaction_pairs_FR"] = {"color": "#003399", "name": "Valid Pairs FR"}
        keys["Self_Cycle_pairs"] = {"color": "#ffad99", "name": "Same Fragment - Self-Circle"}
        keys["Dangling_end_pairs"] = {"color": "#ff5c33", "name": "Same Fragment - Dangling Ends"}
        keys["Religation_pairs"] = {"color": "#cc2900", "name": "Re-ligation"}
        keys["Filtered_pairs"] = {"color": "#661400", "name": "Filtered pairs"}
        keys["Dumped_pairs"] = {"color": "#330a00", "name": "Dumped pairs"}

        # Config for the plot
        config = {
            "id": "hicpro_filtering_plot",
            "title": "HiC-Pro: Filtering Statistics",
            "ylab": "# Read Pairs",
            "cpswitch_counts_label": "Number of Read Pairs",
        }

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_contact_chart(self):
        """Generate the HiC-Pro interaction plot"""

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys["cis_shortRange"] = {"color": "#0039e6", "name": "Unique: cis <= 20Kbp"}
        keys["cis_longRange"] = {"color": "#809fff", "name": "Unique: cis > 20Kbp"}
        keys["trans_interaction"] = {"color": "#009933", "name": "Unique: trans"}
        keys["duplicates"] = {"color": "#a9a2a2", "name": "Duplicate read pairs"}

        # Config for the plot
        config = {
            "id": "hicpro_contact_plot",
            "title": "HiC-Pro: Contact Statistics",
            "ylab": "# Pairs",
            "cpswitch_counts_label": "Number of Pairs",
        }

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_as_chart(self):
        """Generate Allele-specific plot"""

        keys = OrderedDict()
        keys["Valid_pairs_from_ref_genome_(1-1)"] = {"color": "#e6550d", "name": "Genome1 specific read pairs (1-1)"}
        keys["Valid_pairs_from_ref_genome_with_one_unassigned_mate_(0-1/1-0)"] = {
            "color": "#fdae6b",
            "name": "Genome1 with one unassigned mate (0-1/1-0)",
        }
        keys["Valid_pairs_from_alt_genome_(2-2)"] = {"color": "#756bb1", "name": "Genome2 specific read pairs (2-2)"}
        keys["Valid_pairs_from_alt_genome_with_one_unassigned_mate_(0-2/2-0)"] = {
            "color": "#bcbddc",
            "name": "Genome2 with one unassigned mate (0-2/2-0)",
        }
        keys["Valid_pairs_from_alt_and_ref_genome_(1-2/2-1)"] = {
            "color": "#a6611a",
            "name": "Trans homologuous read pairs (1-2/2/1)",
        }
        keys["Valid_pairs_with_both_unassigned_mated_(0-0)"] = {"color": "#cccccc", "name": "Unassigned read pairs"}
        keys["Valid_pairs_with_at_least_one_conflicting_mate_(3-)"] = {
            "color": "#a9a2a2",
            "name": "Conflicting read pairs",
        }

        # check allele-specific analysis was run
        num_samples = 0
        for s_name in self.hicpro_data:
            for k in keys:
                num_samples += sum([1 if k in self.hicpro_data[s_name] else 0])
        if num_samples == 0:
            return False

        # Config for the plot
        config = {
            "id": "hicpro_asan_plot",
            "title": "HiC-Pro: Allele-specific Statistics",
            "ylab": "# Pairs",
            "cpswitch_counts_label": "Number of Pairs",
        }

        return bargraph.plot(self.hicpro_data, keys, config)

    def hicpro_capture_chart(self):
        """Generate Capture Hi-C plot"""

        keys = OrderedDict()
        keys["valid_pairs_on_target_cap_cap"] = {"color": "#0039e6", "name": "Capture-Capture interactions"}
        keys["valid_pairs_on_target_cap_rep"] = {"color": "#809fff", "name": "Capture-Reporter interactions"}
        keys["valid_pairs_off_target"] = {"color": "#cccccc", "name": "Off-target valid pairs"}

        # Check capture info are available
        num_samples = 0
        for s_name in self.hicpro_data:
            for k in keys:
                num_samples += sum([1 if k in self.hicpro_data[s_name] else 0])
        if num_samples == 0:
            return False

        # Config for the plot
        config = {
            "id": "hicpro_cap_plot",
            "title": "HiC-Pro: Capture Statistics",
            "ylab": "# Pairs",
            "cpswitch_counts_label": "Number of Pairs",
        }

        return bargraph.plot(self.hicpro_data, keys, config)
