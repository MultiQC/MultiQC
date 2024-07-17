import logging
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, scatter

log = logging.getLogger(__name__)

VERSION_REGEX = r"# gffcompare v([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    """
    The program `gffcompare` can be used to compare, merge, annotate and estimate accuracy
    of one or more GFF files (the "query" files), when compared with a reference annotation (also provided as GFF).

    The _Sensitivity / Precision_ values are displayed in a single plot,
    different loci levels can be switched by choosing a different dataset.

    :::warning
    Please use `gffcompare` only with single samples.
    Multi-Sample comparisons are not correctly rendered by this MultiQC module.
    :::

    Note that exported data in `multiqc_data/multiqc_gffcompare.{tsv,yaml,json}` only works when
    exporting with YAML or JSON - the default `.tsv` output will not contain any data.
    Please use `-k yaml` or `-k json` to export in a structured format.
    It is hoped to refactor this code in a future release - please submit a PR if you are interested.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="GffCompare",
            anchor="gffcompare",
            href="https://ccb.jhu.edu/software/stringtie/gffcompare.shtml",
            info="Tool to compare, merge and annotate one or more GFF files with a reference annotation in GFF format.",
            doi="10.12688/f1000research.23297.1",
        )

        # Parse stats file
        # Everything is hardcoded with linenumbers, needs to be adjusted if output should change.
        # Other versions of GffCompare with different ouptut formats might not work.

        # TODO (@ewels): It would be nice to rewrite this with regular expressions to make it more robust
        # However, for now I'm ok with it as the output is written by the tool and not stdout / stderr,
        # so the line numbering *should* be fairly portable across systems, at least for the same version of the tool.
        # If you come here because the tool has updated and the module has broken, please refactor the below code.

        self.gffcompare_data = {}
        for f in self.find_log_files("gffcompare"):
            self.add_data_source(f)
            sample = f["s_name"]
            self.gffcompare_data[sample] = {}
            lines = f["f"].splitlines()

            # Version info
            version_match = re.search(VERSION_REGEX, lines[0])
            if version_match:
                self.add_software_version(version_match.group(1), sample)

            ## Transcript and loci numbers:
            # Query
            c = [int(s) for s in lines[5].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]["counts"] = {}
            self.gffcompare_data[sample]["counts"]["query"] = {}
            self.gffcompare_data[sample]["counts"]["query"]["mRNA"] = c[0]
            self.gffcompare_data[sample]["counts"]["query"]["loci"] = c[1]
            self.gffcompare_data[sample]["counts"]["query"]["multi_exon_transcripts"] = c[2]
            c = [int(s) for s in lines[6].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]["counts"]["query"]["multi_transcript_loci"] = c[0]

            # Reference
            c = [int(s) for s in lines[7].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]["counts"]["reference"] = {}
            self.gffcompare_data[sample]["counts"]["reference"]["mRNA"] = c[0]
            self.gffcompare_data[sample]["counts"]["reference"]["loci"] = c[1]
            self.gffcompare_data[sample]["counts"]["reference"]["multi_exon_transcripts"] = c[2]
            c = [int(s) for s in lines[8].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]["counts"]["super_loci"] = c[0]

            ## Accuracy metrics (Sensitivity/Precision)
            self.gffcompare_data[sample]["accuracy"] = {}
            self.gffcompare_data[sample]["missed"] = {}
            self.gffcompare_data[sample]["novel"] = {}
            for line in lines[10:16]:
                split = line.replace("|", "").replace("Intron chain", "Intron_chain").split()
                self.gffcompare_data[sample]["accuracy"][split[0]] = {}
                self.gffcompare_data[sample]["accuracy"][split[0]]["sensitivity"] = float(split[2])
                self.gffcompare_data[sample]["accuracy"][split[0]]["precision"] = float(split[3])

            ## Additional count data
            self.gffcompare_data[sample]["counts"]["matching_intron_chains"] = [
                int(s) for s in lines[17].replace("(", " ").split() if s.isdigit()
            ][0]
            self.gffcompare_data[sample]["counts"]["matching_transcripts"] = [
                int(s) for s in lines[18].replace("(", " ").split() if s.isdigit()
            ][0]
            self.gffcompare_data[sample]["counts"]["matching_loci"] = [
                int(s) for s in lines[19].replace("(", " ").split() if s.isdigit()
            ][0]

            ## Missed/Novel genomic elements
            self.gffcompare_data[sample]["missed"]["Exons"] = [
                int(s) for s in lines[21].replace("/", " ").split() if s.isdigit()
            ]
            self.gffcompare_data[sample]["novel"]["Exons"] = [
                int(s) for s in lines[22].replace("/", " ").split() if s.isdigit()
            ]
            self.gffcompare_data[sample]["missed"]["Introns"] = [
                int(s) for s in lines[23].replace("/", " ").split() if s.isdigit()
            ]
            self.gffcompare_data[sample]["novel"]["Introns"] = [
                int(s) for s in lines[24].replace("/", " ").split() if s.isdigit()
            ]
            self.gffcompare_data[sample]["missed"]["Loci"] = [
                int(s) for s in lines[25].replace("/", " ").split() if s.isdigit()
            ]
            self.gffcompare_data[sample]["novel"]["Loci"] = [
                int(s) for s in lines[26].replace("/", " ").split() if s.isdigit()
            ]

        # Filter to strip out ignored sample names
        self.gffcompare_data = self.ignore_samples(self.gffcompare_data)

        # Raise user warning if no data found
        if len(self.gffcompare_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.gffcompare_data)} reports")

        # Add nothing to general statistics table

        # Write data file
        # TODO (@ewels) - because the data structure is hierarchical, the default tsv output format gives an empty file
        # If / when refactoring the above parsing code, restructure to a flat dict so that this works with tabular output files.
        self.write_data_file(self.gffcompare_data, "multiqc_gffcompare")

        # Report sectionsq
        self.add_section(
            name="Accuracy values",
            description=(
                """
                Displayed are the accuracy values precisiond and sensitivity for different levels of genomic features.
                The metrics are calculated for the comparison of a query and reference GTF file.
                """
            ),
            helptext=(
                """
                Accuracy metrics are calculated as described in [Burset et al. (1996)](http://dx.doi.org/10.1006/geno.1996.0298). Sensitivity is the true positive rate, Precision
                True Positives are query features that agree with features in the reference. The exact definition depends on the feature level:

                - *Base*: True positives are bases reported at the same coordinates.
                - *Exon*: Comparison units are exons that overlap in query and reference with same coordinates.
                - *Intron chain*: True positives are query transcripts for which all introns coordinates match those in the reference.
                - *Transcript*: More stringent then intron chain, all Exon coordinates need to match. Outer exon coordinates (start + end) can vary by 100 bases in default settings
                - *Locus*: Cluster of exons need to match.

                A more in depth description can be found [here](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#levels).
                """
            ),
            anchor="gffcompare_accuracy",
            plot=self.plot_accuracy(),
        )

        self.add_section(
            name="Novel features",
            description="Number of novel features, present in the query data but not found in the reference annotation.",
            anchor="gffcompare_novel",
            plot=self.plot_novel(),
        )

        self.add_section(
            name="Missing features",
            description="False negative features, which are found in the reference annotation but missed (not present) in the query data.",
            anchor="gffcompare_missed",
            plot=self.plot_missed(),
        )

    # Plotting functions
    def plot_accuracy(self):
        """Generate GffCompare accuracy plot"""

        datasets = ["Base", "Exon", "Intron", "Intron_chain", "Transcript", "Locus"]

        pconfig = {
            "id": "gffcompare_accuracy_plot",
            "title": "GffCompare: Accuracy values",
            "ylab": "Precision",
            "xlab": "Sensitivity",
            "ymin": 0,
            "ymax": 1,
            "xmin": 0,
            "xmax": 1,
            "data_labels": [{"name": x, "ylab": "Precision", "xlab": "Sensitivity"} for x in datasets],
        }

        data_classification = [
            {
                sample: {
                    "x": self.gffcompare_data[sample]["accuracy"][dataset]["sensitivity"] / 100,
                    "y": self.gffcompare_data[sample]["accuracy"][dataset]["precision"] / 100,
                    "name": dataset,
                }
                for sample in self.gffcompare_data.keys()
            }
            for dataset in datasets
        ]

        return scatter.plot(data_classification, pconfig)

    # Plot number of novel and missing transcripts
    def plot_novel(self):
        """Generate GffCompare novel elements plot"""

        cats = [
            {
                "novel": {"name": "Novel features", "color": "#f7a35c"},
                "matching": {"name": "Matching features", "color": "#8bbc21"},
            }
        ] * 3

        datasets = ["Exons", "Introns", "Loci"]
        pconfig = {
            "id": "gffcompare_novel_plot",
            "title": "Gffcompare: Novel features",
            "ylab": "Reads",
            "data_labels": datasets,
        }

        data_novel = [
            {
                sample: {
                    "novel": self.gffcompare_data[sample]["novel"][dataset][0],
                    "matching": (
                        self.gffcompare_data[sample]["novel"][dataset][1]
                        - self.gffcompare_data[sample]["novel"][dataset][0]
                    ),
                }
                for sample in self.gffcompare_data.keys()
            }
            for dataset in datasets
        ]

        return bargraph.plot(data_novel, cats, pconfig)

    def plot_missed(self):
        """Generate GffCompare missed elements plot"""

        cats = [
            {
                "missed": {"name": "Missing features", "color": "#f7a35c"},
                "matching": {"name": "Matching features", "color": "#8bbc21"},
            }
        ] * 3

        datasets = ["Exons", "Introns", "Loci"]
        pconfig = {
            "id": "gffcompare_missing_plot",
            "title": "Gffcompare: Missing features",
            "ylab": "Reads",
            "data_labels": datasets,
        }

        data_missed = [
            {
                sample: {
                    "missed": self.gffcompare_data[sample]["missed"][dataset][0],
                    "matching": (
                        self.gffcompare_data[sample]["missed"][dataset][1]
                        - self.gffcompare_data[sample]["missed"][dataset][0]
                    ),
                }
                for sample in self.gffcompare_data.keys()
            }
            for dataset in datasets
        ]

        return bargraph.plot(data_missed, cats, pconfig)
