import logging
import os
import re
from typing import List, Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="QoRTs",
            anchor="qorts",
            href="http://hartleys.github.io/QoRTs/",
            info="Toolkit for analysis, QC, and data management of RNA-Seq datasets.",
            extra="Aids in the detection and identification of errors, biases, and artifacts produced "
            "by paired-end high-throughput RNA-Seq technology. In addition, it can produce "
            "count data designed for use with differential expression and differential "
            "exon usage tools, as well as individual-sample and/or group-summary genome "
            "track files suitable for use with the UCSC genome browser.",
            doi="10.1186/s12859-015-0670-5",
        )

        # Parse logs
        self.qorts_data: Dict = dict()
        for f in self.find_log_files("qorts", filehandles=True):
            self.parse_qorts(f)
            self.add_data_source(f)

        # Remove empty samples
        self.qorts_data = {s: v for s, v in self.qorts_data.items() if len(v) > 0}

        # Filter to strip out ignored sample names
        self.qorts_data = self.ignore_samples(self.qorts_data)

        if len(self.qorts_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.qorts_data)} logs")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.write_data_file(self.qorts_data, "multiqc_qorts")

        # Make plots
        self.qorts_general_stats()
        self.qorts_alignment_barplot()
        self.qorts_splice_loci_barplot()
        self.qorts_splice_events_barplot()
        self.qorts_strandedness_plot()

    def parse_qorts(self, f):
        s_names: List[str] = []
        for line in f["f"]:
            s = line.split("\t")
            if not s_names:
                raw_s_names = s[1:]
                s_names = [self.clean_s_name(s_name, f) for s_name in raw_s_names]
                if len(s_names) <= 2 and raw_s_names[0].endswith("COUNT"):
                    if f["fn"] == "QC.summary.txt":
                        s_names = [self.clean_s_name(os.path.basename(os.path.normpath(f["root"])), f)]
                    else:
                        s_names = [f["s_name"]]
                for s_name in s_names:
                    if s_name in self.qorts_data:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.qorts_data[s_name] = dict()
            else:
                for i, s_name in enumerate(s_names):
                    # Hack to get around Java localisation with commas for decimal places
                    if "," in s[i + 1] and "." not in s[i + 1]:
                        s[i + 1] = s[i + 1].replace(",", ".")
                    self.qorts_data[s_name][s[0]] = float(s[i + 1])
        # Add some extra fields
        for i, s_name in enumerate(s_names):
            if "Genes_Total" in self.qorts_data[s_name] and "Genes_WithNonzeroCounts" in self.qorts_data[s_name]:
                self.qorts_data[s_name]["Genes_PercentWithNonzeroCounts"] = (
                    self.qorts_data[s_name]["Genes_WithNonzeroCounts"] / self.qorts_data[s_name]["Genes_Total"]
                ) * 100.0

    def qorts_general_stats(self):
        """Add columns to the General Statistics table"""
        headers = {
            "Genes_PercentWithNonzeroCounts": {
                "title": "% Genes with Counts",
                "description": "Percent of Genes with Non-Zero Counts",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "NumberOfChromosomesCovered": {
                "title": "Chrs Covered",
                "description": "Number of Chromosomes Covered",
                "format": "{:,.0f}",
            },
        }
        self.general_stats_addcols(self.qorts_data, headers)

    def qorts_alignment_barplot(self):
        """Alignment statistics bar plot"""
        # Specify the order of the different possible categories
        keys = [
            "ReadPairs_UniqueGene_CDS",
            "ReadPairs_UniqueGene_UTR",
            "ReadPairs_AmbigGene",
            "ReadPairs_NoGene_Intron",
            "ReadPairs_NoGene_OneKbFromGene",
            "ReadPairs_NoGene_TenKbFromGene",
            "ReadPairs_NoGene_MiddleOfNowhere",
        ]
        cats = {}
        for k in keys:
            name = k.replace("ReadPairs_", "").replace("_", ": ")
            name = re.sub("([a-z])([A-Z])", r"\g<1> \g<2>", name)
            cats[k] = {"name": name}

        # Config for the plot
        pconfig = {
            "id": "qorts_alignments",
            "title": "QoRTs: Alignment Locations",
            "ylab": "# Read Pairs",
            "cpswitch_counts_label": "Number of Read Pairs",
            "hide_empty": False,
        }

        self.add_section(
            name="Alignments",
            description="This plot displays the rate for which the sample's read-pairs are assigned to the different categories.",
            helptext="""
            The [QoRTs vignette](http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf) describes the categories in this plot as follows:

            * **Unique Gene**: The read-pair overlaps with the exonic segments of one and only one gene. For many
              downstream analyses tools, such as DESeq, DESeq2 and EdgeR, only read-pairs in this category
              are used.
            * **Ambig Gene**: The read-pair overlaps with the exons of more than one gene.
            * **No Gene: Intronic**: The read-pair does not overlap with the exons of any annotated gene, but appears
              in a region that is bridged by an annotated splice junction.
            * **No Gene: One kb From Gene**: The read-pair does not overlap with the exons of any annotated gene, but is
              within 1 kilobase from the nearest annotated gene.
            * **No Gene: Ten kb From Gene**: The read-pair does not overlap with the exons of any annotated gene, but
              is within 10 kilobases from the nearest annotated gene.
            * **No Gene: Middle Of Nowhere**: The read-pair does not overlap with the exons of any annotated gene,
              and is more than 10 kilobases from the nearest annotated gene.

            _What it means and what to look for:_

            Outliers in these plots can indicate biological variations or the presence of large mapping problems.
            They may also suggest the presence of large, highly-expressed, unannotated transcripts or genes.
            """,
            plot=bargraph.plot(self.qorts_data, cats, pconfig),
        )

    def qorts_splice_loci_barplot(self):
        # Specify the order of the different possible categories
        keys = [
            "SpliceLoci_Known_ManyReads",
            "SpliceLoci_Known_FewReads",
            "SpliceLoci_Known_NoReads",
            "SpliceLoci_Novel_ManyReads",
            "SpliceLoci_Novel_FewReads",
        ]
        cats = {}
        for k in keys:
            name = k.replace("SpliceLoci_", "").replace("_", ": ")
            name = re.sub("([a-z])([A-Z])", r"\g<1> \g<2>", name)
            cats[k] = {"name": name}

        # Config for the plot
        pconfig = {
            "id": "qorts_splice_loci",
            "title": "QoRTs: Splice Loci",
            "ylab": "# Splice Loci",
            "cpswitch_counts_label": "Number of Splice Loci",
            "hide_empty": False,
        }

        self.add_section(
            name="Splice Loci",
            description="This plot shows the number of splice junction loci of each type that appear in the sample's reads.",
            helptext="""
            The [QoRTs vignette](http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf) describes the categories in this plot as follows:

            * **Known**: The splice junction locus is found in the supplied transcript annotation gtf file.
            * **Novel**: The splice junction locus is NOT found in the supplied transcript annotation gtf file.
            * **Known: Few reads**: The locus is known, and is only covered by 1-3 read-pairs.
            * **Known: Many reads**: The locus is known, and is covered by 4 or more read-pairs.
            * **Novel: Few reads**: The locus is novel, and is only covered by 1-3 read-pairs.
            * **Novel: Many reads**: The locus is novel, and is covered by 4 or more read-pairs

            _What it means and what to look for:_

            This plot can be used to detect a number of anomalies. For example:
            whether mapping or sequencing artifacts caused a disproportionate discovery of novel splice junctions in
            one sample or batch. It can also be used as an indicator of the comprehensiveness the genome annotation.
            Replicates that are obvious outliers may have sequencing/technical issues causing false detection of splice
            junctions.

            Abnormalities in the splice junction rates are generally a symptom of larger issues which will generally be
            picked up by other metrics. Numerous factors can reduce the efficacy by which aligners map across splice
            junctions, and as such these plots become very important if the intended downstream analyses include
            transcript assembly, transcript deconvolution, differential splicing, or any other form of analysis that in
            some way involves the splice junctions themselves. These plots can be used to assess whether other minor
            abnormalities observed in the other plots are of sufficient severity to impact splice junction mapping and
            thus potentially compromise such analyses.
            """,
            plot=bargraph.plot(self.qorts_data, cats, pconfig),
        )

    def qorts_splice_events_barplot(self):
        # Specify the order of the different possible categories
        keys = [
            "SpliceEvents_KnownLociWithManyReads",
            "SpliceEvents_KnownLociWithFewReads",
            "SpliceEvents_NovelLociWithManyReads",
            "SpliceEvents_NovelLociWithFewReads",
        ]
        cats = {}
        for k in keys:
            name = k.replace("SpliceEvents_", "")
            name = re.sub("([a-z])([A-Z])", r"\g<1> \g<2>", name)
            cats[k] = {"name": name}

        # Config for the plot
        pconfig = {
            "id": "qorts_splice_events",
            "title": "QoRTs: Splice Events",
            "ylab": "# Splice Events",
            "cpswitch_counts_label": "Number of Splice Events",
            "hide_empty": False,
        }

        self.add_section(
            name="Splice Events",
            description="This plot shows the number of splice junction events falling into different junction categories.",
            helptext="""
            From the [QoRTs vignette](http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf):

            A splice junction "event" is one instance of a read-pair bridging a splice junction.
            Some reads may contain multiple splice junction events, some may contain none. If a splice junction appears
            on both reads of a read-pair, this is still only counted as a single "event".

            Note that because different samples/runs may have different total read counts and/or library sizes, this function
            is generally not the best for comparing between samples. In general, the event rates per read-pair should be
            used instead.
            This plot is used to detect whether sample-specific or batch effects have a substantial or biased effect on splice
            junction appearance, either due to differences in the original RNA, or due to artifacts that alter the rate at
            which the aligner maps across splice junctions.

            _What it means and what to look for:_

            This plot is useful for identifying mapping and/or annotation issues,
            and can indicate the comprehensiveness the genome annotation. Replicates that are obvious outliers may
            have sequencing/technical issues causing false detection of splice junctions.
            In general, abnormalities in the splice junction rates are generally a symptom of larger issues which will
            often be picked up by other metrics.
            """,
            plot=bargraph.plot(self.qorts_data, cats, pconfig),
        )

    def qorts_strandedness_plot(self):
        """Make a bar plot showing the reads assigned to each strand"""
        # Specify the order of the different possible categories
        keys = [
            "StrandTest_frFirstStrand",
            "StrandTest_frSecondStrand",
            "StrandTest_ambig_genesFountOnBothStrands",
            "StrandTest_ambig_noGenes",
            "StrandTest_ambig_other",
        ]
        cats = {}
        for k in keys:
            name = k.replace("StrandTest_", "").replace("_", " ").replace("ambig", "ambig:")
            name = re.sub("([a-z])([A-Z])", r"\g<1> \g<2>", name)
            cats[k] = {"name": name.title()}

        # Config for the plot
        pconfig = {
            "id": "qorts_strand_test",
            "title": "QoRTs: Strand Test",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
            "cpswitch_c_active": False,
        }

        self.add_section(
            name="Strandedness",
            description="This plot shows the rate at which reads appear to follow different library-type strandedness rules.",
            helptext="""
            From the [QoRTs vignette](http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf):

            This plot is used to detect whether your data is indeed stranded, and whether you are using the correct
            stranded data library type option. For unstranded libraries, one would expect close to 50-50
            First Strand - Second Strand. For stranded libraries, all points should fall closer to 99% one or the other.
            """,
            plot=bargraph.plot(self.qorts_data, cats, pconfig),
        )
