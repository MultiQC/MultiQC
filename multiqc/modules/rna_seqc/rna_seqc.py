""" MultiQC module to parse output from RNA-SeQC """


import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, violin, heatmap, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="RNA-SeQC",
            anchor="rna_seqc",
            href="https://github.com/getzlab/rnaseqc",
            info="Fast, efficient RNA-Seq metrics for quality control and process optimization",
            doi=["10.1093/bioinformatics/btab135", "10.1093/bioinformatics/bts196"],
        )

        # Parse metrics from RNA-SeQC v1
        self.rna_seqc_metrics = dict()
        for f in self.find_log_files("rna_seqc/metrics_v1", filehandles=True):
            self.parse_metrics_rnaseqc_v1(f)
            self.add_data_source(f, section="v1")

        # Parse metrics from RNA-SeQC v2
        for f in self.find_log_files("rna_seqc/metrics_v2", filehandles=True):
            self.parse_metrics_rnaseqc_v2(f)
            self.add_data_source(f, section="v2")

        # Strip out ignored sample names from metrics
        self.rna_seqc_metrics = self.ignore_samples(self.rna_seqc_metrics)

        # Parse normalised coverage information.
        self.rna_seqc_norm_high_cov = dict()
        self.rna_seqc_norm_medium_cov = dict()
        self.rna_seqc_norm_low_cov = dict()
        for f in self.find_log_files("rna_seqc/coverage"):
            self.parse_coverage(f)
        self.rna_seqc_norm_high_cov = self.ignore_samples(self.rna_seqc_norm_high_cov)
        self.rna_seqc_norm_medium_cov = self.ignore_samples(self.rna_seqc_norm_medium_cov)
        self.rna_seqc_norm_low_cov = self.ignore_samples(self.rna_seqc_norm_low_cov)

        # Parse correlation matrices
        self.rna_seqc_pearson = None
        self.rna_seqc_spearman = None
        num_correlation_samples = 0
        for f in self.find_log_files("rna_seqc/correlation"):
            num_correlation_samples += self.parse_correlation(f)

        # Work out the maximum number of samples that we have across all reports
        num_found = max(
            len(self.rna_seqc_metrics),
            len(self.rna_seqc_norm_high_cov),
            len(self.rna_seqc_norm_medium_cov),
            len(self.rna_seqc_norm_low_cov),
            num_correlation_samples,
        )

        # Stop here if we didn't find anything
        if num_found == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {num_found} samples")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write metrics to a file
        self.write_data_file(self.rna_seqc_metrics, "multiqc_rna_seqc")

        # Generate sections for the report
        self.rnaseqc_general_stats()
        self.transcript_associated_plot()
        self.plot_correlation_heatmap()
        self.strand_barplot()
        self.coverage_lineplot()
        self.bam_statplot()

    def parse_metrics_rnaseqc_v1(self, f):
        """
        Parse the metrics.tsv file from RNA-SeQC version 1.x
        """
        headers = None
        for line in f["f"]:
            s = line.strip().split("\t")
            if headers is None:
                headers = s
            else:
                s_name = s[headers.index("Sample")]
                data = dict()
                for idx, h in enumerate(headers):
                    try:
                        data[h] = float(s[idx])
                    except ValueError:
                        data[h] = s[idx]
                if s_name in self.rna_seqc_metrics:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.rna_seqc_metrics[s_name] = data

    def parse_metrics_rnaseqc_v2(self, f):
        """
        Parse the *.metrics.tsv file from RNA-SeQC version 2.x
        """
        data = {}
        s_name = f["s_name"]
        for line in f["f"]:
            s = line.split("\t")
            if s[0] == "Sample":
                s_name = self.clean_s_name(s[1], f)
            if s[0] == "Total Reads":
                s[0] = "Total Read Number"
            if s[0] == "rRNA Rate":
                s[0] = "rRNA rate"
            try:
                data[s[0]] = float(s[1])
            except ValueError:
                data[s[0]] = s[1].strip()

        if s_name in self.rna_seqc_metrics:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.rna_seqc_metrics[s_name] = data

    def parse_coverage(self, f):
        """Parse the RNA-SeQC Normalised Coverage Files"""
        data = dict()
        s_names = None
        j = 1
        for line in f["f"].splitlines():
            s = line.strip().split("\t")
            if s_names is None:
                s_names = [self.clean_s_name(s_name, f) for s_name in s]
                for s_name in s_names:
                    data[s_name] = dict()
            else:
                for i, v in enumerate(s):
                    data[s_names[i]][j] = float(v)
                j += 1
        if f["fn"] == "meanCoverageNorm_high.txt":
            self.rna_seqc_norm_high_cov.update(data)
        elif f["fn"] == "meanCoverageNorm_medium.txt":
            self.rna_seqc_norm_medium_cov.update(data)
        elif f["fn"] == "meanCoverageNorm_low.txt":
            self.rna_seqc_norm_low_cov.update(data)

    def parse_correlation(self, f):
        """Parse RNA-SeQC correlation matrices"""
        s_names = None
        data = list()
        for line in f["f"].splitlines():
            s = line.strip().split("\t")
            if s_names is None:
                s_names = [x for x in s if x != ""]
            else:
                data.append([float(val) for val in s[1:]])

        # Filter for ignored sample names
        filtered_s_names = list()
        filtered_data = list()
        for idx, s_name in enumerate(s_names):
            s_name = self.clean_s_name(s_name, f)
            if not self.is_ignore_sample(s_name):
                filtered_s_names.append(s_name)
                filtered_data.append(data[idx])

        if f["fn"] == "corrMatrixPearson.txt":
            self.rna_seqc_pearson = (filtered_s_names, filtered_data)
        elif f["fn"] == "corrMatrixSpearman.txt":
            self.rna_seqc_spearman = (filtered_s_names, filtered_data)

        return len(filtered_s_names)

    def rnaseqc_general_stats(self):
        """
        Add key metrics to the general stats table
        """
        headers = {
            "Expression Profiling Efficiency": {
                "title": "% Expression Efficiency",
                "description": "Ratio of exonic reads to total reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "modify": lambda x: float(x) * 100.0,
            },
            "Genes Detected": {
                "title": "# Genes",
                "description": "Number of genes detected with at least 5 reads",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "rRNA rate": {
                "title": "% rRNA Alignment",
                "description": "Ribosomal RNA reads per total reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "Reds",
                "modify": lambda x: float(x) * 100.0,
            },
        }
        self.general_stats_addcols(self.rna_seqc_metrics, headers)

    def transcript_associated_plot(self):
        """Plot a bargraph showing the Transcript-associated reads"""

        # Plot bar graph of groups
        keys = {
            "Exonic Rate": {"name": "Exonic", "color": "#2f7ed8"},
            "Intronic Rate": {"name": "Intronic", "color": "#8bbc21"},
            "Intergenic Rate": {"name": "Intergenic", "color": "#0d233a"},
        }

        # Config for the plot
        pconfig = {
            "id": "rna_seqc_position_plot",
            "title": "RNA-SeQC: Transcript-associated reads",
            "ylab": "Ratio of Reads",
            "cpswitch": False,
            "ymax": 1,
            "ymin": 0,
            "tt_decimals": 3,
            "cpswitch_c_active": False,
        }
        self.add_section(
            name="Transcript-associated reads",
            anchor="rna_seqc_transcript_associated",
            helptext="""
            All of the rates given are per mapped read.

            * **Exonic Rate** is the fraction mapping within exons.
            * **Intronic Rate** is the fraction mapping within introns.
            * **Intergenic Rate** is the fraction mapping in the genomic space between genes.
            """,
            plot=bargraph.plot(self.rna_seqc_metrics, keys, pconfig),
        )

    def plot_correlation_heatmap(self):
        """Return HTML for correlation heatmap"""
        data = None
        corr_type = None
        correlation_type = getattr(config, "rna_seqc", {}).get("default_correlation", "spearman")
        if self.rna_seqc_spearman is not None and correlation_type != "pearson":
            data = self.rna_seqc_spearman
            corr_type = "Spearman"
        elif self.rna_seqc_pearson is not None:
            data = self.rna_seqc_pearson
            corr_type = "Pearson"
        if data is not None:
            pconfig = {
                "id": "rna_seqc_correlation_heatmap",
                "title": f"RNA-SeQC: {corr_type} Sample Correlation",
            }
            self.add_section(
                name=f"{corr_type} Correlation",
                anchor="rseqc-rna_seqc_correlation",
                plot=heatmap.plot(data[1], data[0], data[0], pconfig),
            )

    def strand_barplot(self):
        """Plot a bargraph showing the strandedness of alignments"""
        # Plot bar graph of groups
        keys = ["End 1 Sense", "End 1 Antisense", "End 2 Sense", "End 2 Antisense"]
        # Config for the plot
        pconfig = {
            "id": "rna_seqc_strandedness_plot",
            "title": "RNA-SeQC: Strand Specificity",
            "ylab": "% Reads",
            "cpswitch_counts_label": "# Reads",
            "cpswitch_percent_label": "% Reads",
            "ymin": 0,
            "cpswitch_c_active": False,
        }
        self.add_section(
            name="Strand Specificity",
            anchor="rna_seqc_strand_specificity",
            helptext="""
            End 1/2 Sense are the number of End 1 or 2 reads that were sequenced in the sense direction.

            End 1/2 Antisense are the number of End 1 or 2 reads that were sequenced in the antisense direction
            """,
            plot=bargraph.plot(self.rna_seqc_metrics, keys, pconfig),
        )

    def coverage_lineplot(self):
        """Make HTML for coverage line plots"""
        # Add line graph to section
        data = list()
        data_labels = list()
        if len(self.rna_seqc_norm_high_cov) > 0:
            data.append(self.rna_seqc_norm_high_cov)
            data_labels.append({"name": "High Expressed"})
        if len(self.rna_seqc_norm_medium_cov) > 0:
            data.append(self.rna_seqc_norm_medium_cov)
            data_labels.append({"name": "Medium Expressed"})
        if len(self.rna_seqc_norm_low_cov) > 0:
            data.append(self.rna_seqc_norm_low_cov)
            data_labels.append({"name": "Low Expressed"})
        pconfig = {
            "id": "rna_seqc_gene_body_coverage_plot",
            "title": "RNA-SeQC: Gene Body Coverage",
            "ylab": "% Coverage",
            "xlab": "Gene Body Percentile (5' -> 3')",
            "xmin": 0,
            "xmax": 100,
            "tt_label": "<strong>{point.x}% from 5'</strong>: {point.y:.2f}",
            "data_labels": data_labels,
        }
        if len(data) > 0:
            self.add_section(
                name="Gene Body Coverage",
                anchor="rseqc-rna_seqc_mean_coverage",
                helptext="The metrics are calculated across the transcripts with tiered expression levels.",
                plot=linegraph.plot(data, pconfig),
            )

    def bam_statplot(self):
        pconfig = {"id": "rna_seqc_bam_stat_violin", "title": "RNA-SeQC: Read metrics"}
        columns = [
            "Total Read Number",
            "Alternative Alignments",
            "Chimeric Reads",
            "Duplicate Reads",
            "End 1 Mapped Reads",
            "End 2 Mapped Reads",
            "End 1 Mismatches",
            "End 2 Mismatches",
            "End 1 Sense",
            "End 2 Sense",
            "Ambiguous Reads",
            "High Quality Reads",
            "Low Quality Reads",
            "Mapped Duplicate Reads",
            "Mapped Reads",
            "Mapped Unique Reads",
            "Non-Globin Reads",
            "Non-Globin Duplicate Reads",
            "rRNA Reads",
            "Unique Mapping, Vendor QC Passed Reads",
        ]
        keys = dict()
        for col in columns:
            keys[col] = {
                "title": col,
                "shared_key": "read_count",
                "suffix": config.read_count_prefix,
            }

        self.add_section(
            name="Read Counts",
            anchor="rna_seqc_bam_stat",
            description=f"Number of reads ({config.read_count_desc}) falling into different categories.",
            helptext="Note that many of these statistics are only available from RNA-SeQC v2.x",
            plot=violin.plot(self.rna_seqc_metrics, keys, pconfig),
        )
