import logging
import re

import yaml

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Currently only the "Allele-tagging" and "Allele-sorting" reports are supported.

    The log files from the genome creation steps are not parsed and there are no plots/tables produced from the "SNP coverage" report.

    Differences between the numbers in the tagging and sorting reports are due to paired-end reads.
    For these, if only a single mate in a pair is assigned to a genome then it will "rescue" its mate and both will be "sorted" into that genome (even though only one of them was tagged).
    Conversely, if the mates in a pair are tagged as arising from different genomes, then the pair as a whole is unassignable.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="SNPsplit",
            anchor="SNPsplit",
            target="SNPsplit",
            href="https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/",
            info="Allele-specific alignment sorter. Determines allelic origin of reads that cover known SNP positions",
            doi="10.12688/f1000research.9037.2",
        )

        self.snpsplit_data = dict()

        # Parse log files generated with newer versions of SNPsplit (YAML)
        for f in self.find_log_files("snpsplit/new", filehandles=True):
            parsed_log = self.parse_new_snpsplit_log(f)
            self._save_parsed(parsed_log, f)

        # Parse log files generated with older versions of SNPsplit
        for f in self.find_log_files("snpsplit/old"):
            parsed_log = self.parse_old_snpsplit_log(f)
            self._save_parsed(parsed_log, f)

        # Filter --ignore-samples
        self.snpsplit_data = self.ignore_samples(self.snpsplit_data)

        if len(self.snpsplit_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.snpsplit_data)} reports")

        self.write_data_file(self.snpsplit_data, "multiqc_snpsplit")

        self.add_general_stats()
        self.allele_tagging_section()
        self.allele_sorting_section()

    def _save_parsed(self, parsed, f):
        s_name = self.clean_s_name(parsed[0], f)
        if s_name in self.snpsplit_data:
            log.debug(f"Replacing duplicate sample {s_name}")
        self.snpsplit_data[s_name] = parsed[1]
        self.add_data_source(f, s_name=s_name)
        self.add_software_version(parsed[1].get("version"), s_name)

    @staticmethod
    def parse_new_snpsplit_log(f):
        data = next(yaml.load_all(f["f"], Loader=yaml.SafeLoader))
        flat_data = {}
        for k in data:
            for sk in data[k]:
                key = sk
                for prefix in ["PE_", "SE_", "HiC_"]:
                    if sk.startswith(prefix):
                        key = sk[len(prefix) :]
                flat_key = f"{k.lower()}_{key}"
                flat_data[flat_key] = data[k][sk]
        input_fn = data["Meta"]["infile"]
        flat_data["version"] = data["Meta"]["version"]
        return [input_fn, flat_data]

    def parse_old_snpsplit_log(self, f):
        input_fn = None
        d = dict()

        for line in f["f"].splitlines():
            # Parse the sample name
            match = re.match(r"Input file:\W+'([^']+)'", line)
            if match:
                # Report can have two sections with two input files,
                # the second sorting the output from the tagging
                # Just take the first instance, which is the initial input file
                if input_fn is None:
                    input_fn = match.group(1)
                continue

            # Log format 1: XXX description
            regex_patterns = [
                # Allele-tagging report
                ["tagging_g1", r"(\d+) reads were specific for genome 1"],
                ["tagging_g2", r"(\d+) reads were specific for genome 2"],
                ["tagging_unassignable", r"(\d+) reads were unassignable"],
                ["tagging_bizarre", r"(\d+) contained conflicting allele-specific SNPs"],
                ["tagging_unaligned", r"Reads were unaligned and hence skipped: (\d+)"],
                ["tagging_CT_positions_skipped", r"(\d+) reads that were unassignable contained C>T SNPs"],
                # Allele-specific sorting report
                ["sorting_conflicting", r"Reads contained conflicting SNP information:\W+(\d+)"],
            ]
            for k, regex in regex_patterns:
                match = re.match(regex, line)
                if match:
                    d[k] = int(match.group(1))
                    break

            # Format 2: Decription: XXX
            sorting_patterns = [
                # Tagging meta
                ["tagging_SNP_annotation", "SNP annotation file"],
                ["tagging_SNPs_stored", "SNPs stored in total"],
                ["tagging_N_was_known_SNP", "N was present in the list of known SNPs"],
                ["tagging_N_not_known", "N was not present in the list of SNPs"],
                # Allele sorting
                ["sorting_genome1", "Reads were specific for genome 1"],
                ["sorting_genome2", "Reads were specific for genome 2"],
                ["sorting_unassignable", "Reads were unassignable"],
                ["sorting_conflicting", "Read pairs contained conflicting SNP information"],
                # Hi-C data
                ["sorting_genome1_G1_G1", "Read pairs were specific for genome 1 (G1/G1)"],
                ["sorting_genome2_G2_G2", "Read pairs were specific for genome 2 (G2/G2)"],
                ["sorting_unassignable_UA_UA", "Read pairs were unassignable (UA/UA):"],
                ["sorting_G1_UA_total", "Read pairs were a mix of G1 and UA"],
                ["sorting_G2_UA_total", "Read pairs were a mix of G2 and UA"],
                ["sorting_G1_G2_total", "Read pairs were a mix of G1 and G2"],
            ]
            for k, pattern in sorting_patterns:
                if line.startswith(pattern):
                    try:
                        d[k] = int(line.split("\t")[-1].split()[0])
                    except ValueError:
                        d[k] = line.split("\t")[-1].split()[0]
                    break

            if "tagging_N_was_known_SNP" in d and "tagging_N_not_known" in d:
                n_total = d["tagging_N_was_known_SNP"] + d["tagging_N_not_known"]
                d["tagging_percent_N_was_known_SNP"] = (d["tagging_N_was_known_SNP"] / float(n_total)) * 100.0

        return [input_fn, d]

    def add_general_stats(self):
        """Add some columns to the General Statistics table at the top of the report"""
        headers = {
            "tagging_SNP_annotation": {
                "title": "SNP annotation",
                "description": "Annotation file used for differentiating genomes",
                "scale": False,
                "modify": lambda x: f"<code>{x}</code>",
                "hidden": True,
            },
            "tagging_percent_N_was_known_SNP": {
                "title": "% Ns known SNP",
                "description": "Percentage of detected SNPs in the sample that were also present in the annotation",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
            "tagging_SNPs_stored": {
                "title": "SNPs stored",
                "description": "Total number of SNPs used for the analysis",
                "scale": "PRGn",
                "format": "{:,.0f}",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.snpsplit_data, headers)

    def allele_tagging_section(self):
        """Allele-tagging report"""
        cats = {
            "tagging_g1": {"name": "Genome 1"},
            "tagging_g2": {"name": "Genome 2"},
            "tagging_unassignable": {"name": "Not assigned"},
            "tagging_bizarre": {"name": "Conflicting SNPs"},
            "tagging_unaligned": {"name": "Unaligned reads"},
            "tagging_CT_positions_skipped": {"name": "C->T SNP"},
        }

        # Subtract C->T from unassignable
        plot_data = {}
        for s_name, ss_data in self.snpsplit_data.items():
            if "tagging_unassignable" and "tagging_CT_positions_skipped" in ss_data:
                ss_data["tagging_unassignable"] -= ss_data["tagging_CT_positions_skipped"]
            plot_data[s_name] = ss_data

        pconfig = {
            "id": "snpsplit-allele-tagging-plot",
            "title": "SNPsplit: Allele-tagging report",
            "ylab": "Reads",
            "cpswitch_counts_label": "Reads",
        }

        self.add_section(
            name="Allele-tagging report",
            description="Per-sample metrics of how many reads were assigned to each genome.",
            helptext="""
                Allele tagging works on a per-read basis. The results may therefore differ considerably
                from the Allele-specific sorting results if the samples were paired-end or Hi-C samples.

                For single hybrid genomes, Genome 1-specific reads are specific for the reference sequence
                and Genome 2-specific reads are specific for the straing specified with `--strain`.

                For dual hybrid genomes, Genome 1-specific reads are specific for the strain specified
                with `--strain`, and Genome 2-specific reads are specific for the straing specified with `--strain2`.

                Bar graph categories are:

                * `Genome 1`: Reads assigned to Genome 1
                * `Genome 2`: Reads assigned to Genome 2
                * `Not assigned`: Reads don't overlap a SNP
                * `Conflicting SNPs`: Reads contained allele-specific information for both alleles within the same read
                * `Unaligned reads`: Reads aren't aligned
                * `No match`: Reads overlap informative SNPs, but don't contain the expected nucleotide for either genome
                * `C->T SNP`: (Bisulfite sequencing data only) Read SNPs involved some form of C->T transition, rendering it non-informative for allele-assignment
                """,
            plot=bargraph.plot(plot_data, cats, pconfig),
        )

    def allele_sorting_section(self):
        """Allele-specific sorting report"""
        cats = {
            "sorting_genome1": {"name": "Genome 1"},
            "sorting_genome2": {"name": "Genome 2"},
            "sorting_unassignable": {"name": "Not assigned"},
            "sorting_conflicting": {"name": "Conflicting SNPs"},
            "sorting_genome1_G1_G1": {"name": "Genome 1 / Genome 1"},
            "sorting_genome2_G2_G2": {"name": "Genome 2 / Genome 2"},
            "sorting_unassignable_UA_UA": {"name": "Unassignable / Unassignable"},
            "sorting_G1_UA_total": {"name": "Genome 1 / unassignable"},
            "sorting_G2_UA_total": {"name": "Genome 2 / unassignable"},
            "sorting_G1_G2_total": {"name": "Different genomes"},
        }
        # HiC only

        pconfig = {
            "id": "snpsplit-sorting-plot",
            "title": "SNPsplit: Allele-specific sorting",
            "ylab": "Reads",
            "cpswitch_counts_label": "Reads",
        }

        self.add_section(
            name="Allele-specific sorting",
            description="Per-sample metrics of how reads and pairs of reads were sorted into each genome.",
            helptext="""
                Bargraph categories are:

                * `Genome 1`: Reads assigned to Genome 1
                * `Genome 2`: Reads assigned to Genome 2
                * `Not assigned`: Reads don't overlap a SNP
                * `Conflicting SNPs`: Reads contained allele-specific information for both alleles within the same read

                For HiC data, categories are:

                * `Genome 1 / Genome 1`: Pairs with both reads specific to genome 1
                * `Genome 2 / Genome 2`: Pairs with both reads specific to genome 2
                * `Unassignable / Unassignable`: Pairs with both reads unassignable (not overlapping a SNP)
                * `Genome 1 / unassignable`: One paired-end read assigned to Genome 1, one unassignable (doesn't overlap a SNP)
                * `Genome 2 / unassignable`: One paired-end read assigned to Genome 2, one unassignable (doesn't overlap a SNP)
                * `Different genomes`: Paired-end reads assigned to different genomes
                * `Conflicting SNPs`: Reads contained allele-specific information for both alleles within the same read

                Allele-specific sorting takes both reads of read pairs or Hi-C samples into account.

                Note that metrics here may differ from those in the allele-tagging report.
                This occurs when paired-end reads are used, since 'tagging' only one read in
                a pair as arising from one genome can suffice in both reads being sorted there.

                Similarly, if two reads in a pair are tagged as arising from different genomes
                then the pair becomes unassignable.
            """,
            plot=bargraph.plot(self.snpsplit_data, cats, pconfig),
        )
