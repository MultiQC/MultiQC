import logging
import os
import re
import gzip
import colorsys
from typing import Dict, List, Optional

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
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

        # Parse version from log files
        for f in self.find_log_files("qorts/log"):
            self.parse_qorts_log(f)

        self.write_data_file(self.qorts_data, "multiqc_qorts")

        # Make plots
        self.qorts_general_stats()
        self.qorts_alignment_barplot()
        self.qorts_dropped_rates_plot()
        self.qorts_missingness_rate_plot()
        self.qorts_mapping_rates_plot()
        self.qorts_norm_factors_plot()

        self.qorts_biotype_rates_plot()
        self.qorts_chrom_type_rates_plot()
        self.qorts_strandedness_plot()

        self.qorts_splice_loci_barplot()
        self.qorts_splice_events_barplot()
        self.qorts_splice_events_per_read_pair_plot()

        self.qorts_phred_quality_profile_plot()
        self.qorts_gc_content_distribution_plot()
        self.qorts_n_rate_plot()
        self.qorts_base_composition_plot()
        self.qorts_insert_size_plot()
        self.qorts_read_length_distribution_plot()

        self.qorts_clipping_distribution_plots()
        self.qorts_leading_clipped_nucleotide_rates_plot()
        self.qorts_trailing_clipped_nucleotide_rates_plot()
        self.qorts_aligned_nucleotide_rates_plot()
        self.qorts_nvc_clip_match_by_clip_position_plot()

        self.qorts_gene_body_coverage_plot()
        self.qorts_gene_cdf_plot()

        self.qorts_cigar_length_distribution_plot()
        self.qorts_cigar_op_by_cycle_plot()
        self.qorts_deletion_profile_by_cycle_plot()
        self.qorts_insertion_profile_by_cycle_plot()
        self.qorts_splicing_profile_by_cycle_plot()

        self.qorts_on_target_plot()
        self.qorts_overlap_plot()
        self.qorts_reference_mismatch_family_plot()

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
            if "TOTAL_READ_PAIRS" in self.qorts_data[s_name] and "READ_PAIR_OK" in self.qorts_data[s_name]:
                total_pairs = self.qorts_data[s_name]["TOTAL_READ_PAIRS"]
                if total_pairs > 0:
                    self.qorts_data[s_name]["AlignmentRate"] = (
                        self.qorts_data[s_name]["READ_PAIR_OK"] / total_pairs
                    ) * 100.0
            if "READ_PAIR_OK" in self.qorts_data[s_name] and "ReadPairs_UniqueGene" in self.qorts_data[s_name]:
                aligned_pairs = self.qorts_data[s_name]["READ_PAIR_OK"]
                if aligned_pairs > 0:
                    self.qorts_data[s_name]["ReadPairs_UniqueGene_Pct"] = (
                        self.qorts_data[s_name]["ReadPairs_UniqueGene"] / aligned_pairs
                    ) * 100.0
            if "READ_PAIR_OK" in self.qorts_data[s_name] and "ReadPairs_NoGene" in self.qorts_data[s_name]:
                aligned_pairs = self.qorts_data[s_name]["READ_PAIR_OK"]
                if aligned_pairs > 0:
                    self.qorts_data[s_name]["ReadPairs_NoGene_Pct"] = (
                        self.qorts_data[s_name]["ReadPairs_NoGene"] / aligned_pairs
                    ) * 100.0

    def qorts_general_stats(self):
        """Add columns to the General Statistics table"""
        headers = {
            "TOTAL_READ_PAIRS": {
                "title": "Total Pairs",
                "description": "Total read pairs",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "READ_PAIR_OK": {
                "title": "Aligned Pairs",
                "description": "Read pairs retained after QoRTs filtering",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "AlignmentRate": {
                "title": "Alignment Rate",
                "description": "READ_PAIR_OK / TOTAL_READ_PAIRS",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": False,
            },
            "ReadPairs_UniqueGene": {
                "title": "Unique Gene",
                "description": "Read pairs assigned to unique genes",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "ReadPairs_UniqueGene_Pct": {
                "title": "% Unique Gene",
                "description": "ReadPairs_UniqueGene / READ_PAIR_OK",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": False,
            },
            "ReadPairs_NoGene_Pct": {
                "title": "% No Gene",
                "description": "ReadPairs_NoGene / READ_PAIR_OK",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn-rev",
                "hidden": False,
            },
            "Genes_PercentWithNonzeroCounts": {
                "title": "% Genes with Counts",
                "description": "Percent of Genes with Non-Zero Counts",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": False,
            },
            "NumberOfChromosomesCovered": {
                "title": "Chrs Covered",
                "description": "Number of Chromosomes Covered",
                "format": "{:,.0f}",
                "hidden": False,
            },
        }
        self.general_stats_addcols(self.qorts_data, headers)

    def qorts_compact_stats_table(self):
        rows = {}
        for sample, metrics in self.qorts_data.items():
            row = {}
            if "AlignmentRate" in metrics:
                row["alignment_rate"] = metrics["AlignmentRate"]
            if "ReadPairs_UniqueGene_Pct" in metrics:
                row["unique_gene_pct"] = metrics["ReadPairs_UniqueGene_Pct"]
            if "ReadPairs_NoGene_Pct" in metrics:
                row["no_gene_pct"] = metrics["ReadPairs_NoGene_Pct"]
            if "Genes_PercentWithNonzeroCounts" in metrics:
                row["genes_with_counts_pct"] = metrics["Genes_PercentWithNonzeroCounts"]
            if "TOTAL_READ_PAIRS" in metrics:
                row["total_pairs"] = metrics["TOTAL_READ_PAIRS"]
            if row:
                rows[sample] = row

        if not rows:
            return

        headers = {
            "alignment_rate": {
                "title": "Alignment %",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "format": "{:,.2f}",
            },
            "unique_gene_pct": {
                "title": "Unique gene %",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "format": "{:,.2f}",
            },
            "no_gene_pct": {
                "title": "No gene %",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "format": "{:,.2f}",
            },
            "genes_with_counts_pct": {
                "title": "Genes with counts %",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "format": "{:,.2f}",
            },
            "total_pairs": {
                "title": "Total pairs",
                "format": "{:,.0f}",
                "shared_key": "read_count",
            },
        }

        self.add_section(
            name="Compact stats",
            description="Tool-specific QoRTs KPI table for quick per-sample review.",
            plot=table.plot(
                rows,
                headers,
                {
                    "id": "qorts_compact_stats",
                    "title": "QoRTs: Compact Stats",
                },
            ),
        )

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
            "cpswitch_counts_label": "Count",
            "cpswitch_percent_label": "Percent",
            "hide_zero_cats": False,
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
            "cpswitch_counts_label": "Count",
            "cpswitch_percent_label": "Percent",
            "hide_zero_cats": False,
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
            "cpswitch_counts_label": "Count",
            "cpswitch_percent_label": "Percent",
            "hide_zero_cats": False,
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

    def qorts_splice_events_per_read_pair_plot(self):
        data = {}
        for sample, metrics in self.qorts_data.items():
            read_pairs = metrics.get("READ_PAIR_OK")
            if not read_pairs or read_pairs <= 0:
                continue

            row = {}
            for k in [
                "SpliceEvents",
                "SpliceEvents_KnownLoci",
                "SpliceEvents_KnownLociWithManyReads",
                "SpliceEvents_KnownLociWithFewReads",
                "SpliceEvents_NovelLoci",
                "SpliceEvents_NovelLociWithManyReads",
                "SpliceEvents_NovelLociWithFewReads",
            ]:
                if k in metrics:
                    row[k] = metrics[k] / read_pairs
            if row:
                data[sample] = row

        if not data:
            return

        cats = {
            "SpliceEvents": {"name": "All events / read-pair"},
            "SpliceEvents_KnownLoci": {"name": "Known loci / read-pair"},
            "SpliceEvents_KnownLociWithManyReads": {"name": "Known many reads / read-pair"},
            "SpliceEvents_KnownLociWithFewReads": {"name": "Known few reads / read-pair"},
            "SpliceEvents_NovelLoci": {"name": "Novel loci / read-pair"},
            "SpliceEvents_NovelLociWithManyReads": {"name": "Novel many reads / read-pair"},
            "SpliceEvents_NovelLociWithFewReads": {"name": "Novel few reads / read-pair"},
        }
        cats = {k: v for k, v in cats.items() if any(k in row for row in data.values())}

        self.add_section(
            name="Splice junction events per read-pair",
            description="QoRTs splice-junction events normalized by mapped read pairs.",
            plot=bargraph.plot(
                data,
                cats,
                {
                    "id": "qorts_splice_events_per_read_pair",
                    "title": "QoRTs: Splice Junction Events per Read-Pair",
                    "ylab": "Events per read-pair",
                    "y_decimals": 3,
                    "tt_decimals": 3,
                    "hide_zero_cats": False,
                },
            ),
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
            "cpswitch_counts_label": "Count",
            "cpswitch_percent_label": "Percent",
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

    def parse_qorts_log(self, f):
        """Parse QoRTs log files to extract version information."""
        match = re.search(r"Starting QoRTs v(\d+\.\d+[\.\d]*)", f["f"])
        if match:
            # Derive sample name from directory (same logic as parse_qorts for single-sample)
            s_name = self.clean_s_name(os.path.basename(os.path.normpath(f["root"])), f)
            # Only associate with sample if it exists in our data
            if s_name not in self.qorts_data:
                s_name = None
            self.add_software_version(match.group(1), s_name)

    def _sample_name_from_file(self, f):
        return self.clean_s_name(os.path.basename(os.path.normpath(f["root"])), f)

    def _iter_qorts_files(self, key: str):
        # Use filehandles to guarantee compressed QoRTs files are fully readable.
        return self.find_log_files(key, filehandles=True)

    def _parse_tsv_rows(self, f):
        raw = f["f"]
        if hasattr(raw, "seek"):
            raw.seek(0)
        if hasattr(raw, "read"):
            try:
                content = raw.read()
            except UnicodeDecodeError:
                file_path = os.path.join(f["root"], f["fn"]) if f.get("root") and f.get("fn") else None
                if file_path:
                    with open(file_path, "rb") as fh:
                        content = fh.read()
                else:
                    raise
        else:
            content = raw
        if isinstance(content, bytes):
            if content.startswith(b"\x1f\x8b"):
                content = gzip.decompress(content)
            content = content.decode("utf-8", errors="replace")

        lines = [line.strip() for line in str(content).splitlines() if line.strip() and not line.strip().startswith("#")]
        if len(lines) < 2:
            return None, []
        headers = lines[0].split("\t")
        rows = [line.split("\t") for line in lines[1:] if line]
        return headers, rows

    def _find_column(self, headers: List[str], patterns: List[str]) -> Optional[int]:
        for i, h in enumerate(headers):
            h_norm = h.lower()
            if any(re.search(p, h_norm) for p in patterns):
                return i
        return None

    def _safe_float(self, value: str) -> Optional[float]:
        try:
            return float(value)
        except ValueError:
            return None

    def _aggregate_by_xy(self, rows, x_idx: int, y_idx: int):
        out: Dict[float, float] = {}
        for row in rows:
            if len(row) <= max(x_idx, y_idx):
                continue
            x = self._safe_float(row[x_idx])
            y = self._safe_float(row[y_idx])
            if x is None or y is None:
                continue
            out[x] = out.get(x, 0.0) + y
        return out

    def _series_to_percent(self, series: Dict[float, float]) -> Dict[float, float]:
        total = sum(series.values())
        if total <= 0:
            return {}
        return {x: (100.0 * y / total) for x, y in series.items()}

    def _collect_summary_keys(self, regexes: List[str]) -> List[str]:
        all_keys = sorted({k for s in self.qorts_data.values() for k in s.keys()})
        return [k for k in all_keys if any(re.search(rx, k, re.IGNORECASE) for rx in regexes)]

    def _key_label(self, key: str) -> str:
        return re.sub(r"_+", " ", key).strip()

    def _summary_bar_section(
        self,
        section_name: str,
        section_desc: str,
        section_id: str,
        section_title: str,
        regexes: List[str],
        scale_to_percent: bool = False,
        counts_label: str = "Absolute value",
        lowercase_labels: bool = False,
    ):
        keys = self._collect_summary_keys(regexes)
        if not keys:
            return

        data = {}
        for sample, metrics in self.qorts_data.items():
            row = {}
            for k in keys:
                if k not in metrics:
                    continue
                v = metrics[k]
                if scale_to_percent and v <= 1.2:
                    v *= 100.0
                row[k] = v
            if row:
                data[sample] = row

        if not data:
            return

        cats = {
            k: {"name": self._key_label(k).lower() if lowercase_labels else self._key_label(k)} for k in keys
        }
        pconfig = {
            "id": section_id,
            "title": section_title,
            "ylab": "Value",
            "cpswitch_counts_label": counts_label,
            "cpswitch_percent_label": "Percent",
            "hide_zero_cats": False,
        }
        if scale_to_percent:
            pconfig["cpswitch_c_active"] = False
            pconfig["suffix"] = "%"

        self.add_section(name=section_name, description=section_desc, plot=bargraph.plot(data, cats, pconfig))

    def qorts_phred_quality_profile_plot(self):
        phred_profiles: Dict[str, Dict[float, Dict[str, float]]] = {}

        for key in ["qorts/quals_r1", "qorts/quals_r2"]:
            for f in self._iter_qorts_files(key):
                headers, rows = self._parse_tsv_rows(f)
                if not headers:
                    continue

                sample = self._sample_name_from_file(f)
                read_suffix = " R1" if key.endswith("r1") else " R2"
                sample_key = f"{sample}{read_suffix}"
                x_idx = self._find_column(headers, [r"cycle", r"readpos", r"position", r"base", r"readlen"])
                min_idx = self._find_column(headers, [r"\bmin\b"])
                q1_idx = self._find_column(headers, [r"q1", r"quartile.*1", r"lower.*quartile", r"p25"])
                med_idx = self._find_column(headers, [r"median", r"q2", r"p50"])
                q3_idx = self._find_column(headers, [r"q3", r"quartile.*3", r"upper.*quartile", r"p75"])
                max_idx = self._find_column(headers, [r"\bmax\b"])

                if x_idx is None or med_idx is None:
                    continue

                profile = phred_profiles.setdefault(sample_key, {})
                for row in rows:
                    if len(row) <= max(x_idx, med_idx):
                        continue
                    x = self._safe_float(row[x_idx])
                    med = self._safe_float(row[med_idx])
                    if x is None or med is None:
                        continue
                    profile.setdefault(x, {})
                    profile[x]["Median"] = med
                    if min_idx is not None and len(row) > min_idx:
                        val = self._safe_float(row[min_idx])
                        if val is not None:
                            profile[x]["Min"] = val
                    if q1_idx is not None and len(row) > q1_idx:
                        val = self._safe_float(row[q1_idx])
                        if val is not None:
                            profile[x]["Q1"] = val
                    if q3_idx is not None and len(row) > q3_idx:
                        val = self._safe_float(row[q3_idx])
                        if val is not None:
                            profile[x]["Q3"] = val
                    if max_idx is not None and len(row) > max_idx:
                        val = self._safe_float(row[max_idx])
                        if val is not None:
                            profile[x]["Max"] = val

        if not phred_profiles:
            return

        metrics = ["Min", "Q1", "Median", "Q3", "Max"]
        data_by_metric: List[Dict[str, Dict[float, float]]] = []
        data_labels: List[Dict[str, str]] = []
        for metric in metrics:
            metric_data: Dict[str, Dict[float, float]] = {}
            for sample, points in phred_profiles.items():
                series = {x: vals[metric] for x, vals in points.items() if metric in vals}
                if series:
                    metric_data[sample] = series
            if metric_data:
                data_by_metric.append(metric_data)
                data_labels.append({"name": metric, "ylab": "Phred score"})

        if not data_by_metric:
            return

        self.add_section(
            name="Phred quality profile",
            description="QoRTs per-cycle phred quality score profiles.",
            plot=linegraph.plot(
                data_by_metric,
                {
                    "id": "qorts_phred_quality_profile",
                    "title": "QoRTs: Phred Quality Score by Read Position",
                    "xlab": "Read position",
                    "ylab": "Phred score",
                    "ymin": 0,
                    "data_labels": data_labels,
                },
            ),
        )

    def _base_rate_data_by_base(self, keys: List[str], value_mode: str = "percent"):
        base_names = ["A", "C", "G", "T", "N"]
        data_by_base: List[Dict[str, Dict[float, float]]] = []
        data_labels: List[Dict[str, str]] = []

        for base_name in base_names:
            base_plot_data: Dict[str, Dict[float, float]] = {}
            for key in keys:
                for f in self._iter_qorts_files(key):
                    headers, rows = self._parse_tsv_rows(f)
                    if not headers:
                        continue
                    pos_idx = self._find_column(headers, [r"readpos", r"cycle", r"position"])
                    base_idx = self._find_column(headers, [r"base"])
                    ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
                    if pos_idx is None or base_idx is None or ct_idx is None:
                        continue

                    by_pos: Dict[float, Dict[str, float]] = {}
                    for row in rows:
                        if len(row) <= max(pos_idx, base_idx, ct_idx):
                            continue
                        pos = self._safe_float(row[pos_idx])
                        ct = self._safe_float(row[ct_idx])
                        base = row[base_idx].upper()
                        if pos is None or ct is None:
                            continue
                        by_pos.setdefault(pos, {})
                        by_pos[pos][base] = by_pos[pos].get(base, 0.0) + ct

                    series: Dict[float, float] = {}
                    for pos, counts in by_pos.items():
                        if value_mode == "count":
                            series[pos] = counts.get(base_name, 0.0)
                        else:
                            total = sum(counts.values())
                            if total > 0:
                                series[pos] = 100.0 * counts.get(base_name, 0.0) / total

                    if series:
                        suffix = " R1" if key.endswith("r1") else " R2"
                        sample = f"{self._sample_name_from_file(f)}{suffix}"
                        base_plot_data[sample] = series

            if base_plot_data:
                data_by_base.append(base_plot_data)
                if value_mode == "count":
                    data_labels.append({"name": f"{base_name} Count", "ylab": "Nucleotide count"})
                else:
                    data_labels.append({"name": f"{base_name} %", "ylab": "Nucleotide rate"})

        return data_by_base, data_labels

    def qorts_gc_content_distribution_plot(self):
        data_percent = {}
        for f in self._iter_qorts_files("qorts/gc_by_read"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            x_idx = self._find_column(headers, [r"gc", r"num_bases_gc"])
            y_idx = self._find_column(headers, [r"\bct\b", r"count"])
            if x_idx is None or y_idx is None:
                continue

            x_header = headers[x_idx].lower()
            x_is_gc_count = "num_bases_gc" in x_header or re.search(r"num.*gc", x_header) is not None

            series = {}
            for row in rows:
                if len(row) <= max(x_idx, y_idx):
                    continue
                x = self._safe_float(row[x_idx])
                y = self._safe_float(row[y_idx])
                if x is not None and y is not None:
                    series[x] = y
            if series:
                if x_is_gc_count:
                    read_len = max(series.keys())
                    if read_len > 0:
                        series = {100.0 * x / read_len: y for x, y in series.items()}

                sample = self._sample_name_from_file(f)
                series_pct = self._series_to_percent(series)
                if series_pct:
                    data_percent[sample] = series_pct

        if not data_percent:
            return

        self.add_section(
            name="GC content distribution",
            description="QoRTs GC content distribution by read.",
            plot=linegraph.plot(
                data_percent,
                {
                    "id": "qorts_gc_content_distribution",
                    "title": "QoRTs: GC Content Distribution",
                    "xlab": "GC content",
                    "xsuffix": "%",
                    "ylab": "Read fraction",
                    "ysuffix": "%",
                    "xmin": 0,
                    "xmax": 100,
                    "ymin": 0,
                },
            ),
        )

    def qorts_n_rate_plot(self):
        data_percent = {}
        for key in ["qorts/nvc_raw_r1", "qorts/nvc_raw_r2"]:
            for f in self._iter_qorts_files(key):
                headers, rows = self._parse_tsv_rows(f)
                if not headers:
                    continue
                pos_idx = self._find_column(headers, [r"readpos", r"cycle", r"position"])
                base_idx = self._find_column(headers, [r"base"])
                ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
                if pos_idx is None or base_idx is None or ct_idx is None:
                    continue

                by_pos: Dict[float, Dict[str, float]] = {}
                for row in rows:
                    if len(row) <= max(pos_idx, base_idx, ct_idx):
                        continue
                    pos = self._safe_float(row[pos_idx])
                    ct = self._safe_float(row[ct_idx])
                    base = row[base_idx].upper()
                    if pos is None or ct is None:
                        continue
                    by_pos.setdefault(pos, {})
                    by_pos[pos][base] = by_pos[pos].get(base, 0.0) + ct

                n_rate = {}
                for pos, counts in by_pos.items():
                    n_ct = counts.get("N", 0.0)
                    total = sum(counts.values())
                    if total > 0:
                        n_rate[pos] = 100.0 * n_ct / total

                if n_rate:
                    suffix = " R1" if key.endswith("r1") else " R2"
                    sample = self._sample_name_from_file(f) + suffix
                    data_percent[sample] = n_rate

        if not data_percent:
            return

        self.add_section(
            name="N-rate by read position",
            description="QoRTs per-cycle N-rate.",
            plot=linegraph.plot(
                data_percent,
                {
                    "id": "qorts_n_rate_by_cycle",
                    "title": "QoRTs: N-Rate by Read Position",
                    "xlab": "Read position",
                    "ylab": "N-rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                },
            ),
        )

    def qorts_base_composition_plot(self):
        data_by_base, data_labels = self._base_rate_data_by_base(
            ["qorts/nvc_raw_r1", "qorts/nvc_raw_r2"], value_mode="percent"
        )

        if not data_by_base:
            return

        self.add_section(
            name="Base composition by read position",
            description="QoRTs nucleotide composition by cycle for each read.",
            plot=linegraph.plot(
                data_by_base,
                {
                    "id": "qorts_base_composition_by_cycle",
                    "title": "QoRTs: Base Composition by Read Position",
                    "xlab": "Read position",
                    "ylab": "Nucleotide rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                    "ymax": 100,
                    "data_labels": data_labels,
                },
            ),
        )

    def qorts_leading_clipped_nucleotide_rates_plot(self):
        data_by_base, data_labels = self._base_rate_data_by_base(
            ["qorts/nvc_lead_clip_r1", "qorts/nvc_lead_clip_r2"], value_mode="percent"
        )
        if not data_by_base:
            return
        self.add_section(
            name="Leading-clipped nucleotide rates",
            description="QoRTs nucleotide rates by read position from leading-clipped reads.",
            plot=linegraph.plot(
                data_by_base,
                {
                    "id": "qorts_leading_clipped_nucleotide_rates",
                    "title": "QoRTs: Leading-Clipped Nucleotide Rates",
                    "xlab": "Read position",
                    "ylab": "Nucleotide rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                    "ymax": 100,
                    "data_labels": data_labels,
                },
            ),
        )

    def qorts_trailing_clipped_nucleotide_rates_plot(self):
        data_by_base, data_labels = self._base_rate_data_by_base(
            ["qorts/nvc_tail_clip_r1", "qorts/nvc_tail_clip_r2"], value_mode="percent"
        )
        if not data_by_base:
            return
        self.add_section(
            name="Trailing-clipped nucleotide rates",
            description="QoRTs nucleotide rates by read position from trailing-clipped reads.",
            plot=linegraph.plot(
                data_by_base,
                {
                    "id": "qorts_trailing_clipped_nucleotide_rates",
                    "title": "QoRTs: Trailing-Clipped Nucleotide Rates",
                    "xlab": "Read position",
                    "ylab": "Nucleotide rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                    "ymax": 100,
                    "data_labels": data_labels,
                },
            ),
        )

    def qorts_aligned_nucleotide_rates_plot(self):
        data_by_base, data_labels = self._base_rate_data_by_base(
            ["qorts/nvc_minus_clipping_r1", "qorts/nvc_minus_clipping_r2"], value_mode="percent"
        )
        if not data_by_base:
            return
        self.add_section(
            name="Aligned nucleotide rate by read position",
            description="QoRTs nucleotide rates by read position after clipping-aligned filtering.",
            plot=linegraph.plot(
                data_by_base,
                {
                    "id": "qorts_aligned_nucleotide_rates",
                    "title": "QoRTs: Aligned Nucleotide Rate by Read Position",
                    "xlab": "Read position",
                    "ylab": "Nucleotide rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                    "ymax": 100,
                    "data_labels": data_labels,
                },
            ),
        )

    def qorts_insert_size_plot(self):
        data_percent = {}
        for f in self._iter_qorts_files("qorts/insert_size"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            x_idx = self._find_column(headers, [r"insert", r"size", r"len"])
            y_idx = self._find_column(headers, [r"\bct\b", r"count"])
            if x_idx is None or y_idx is None:
                continue

            counts_by_insert_size = {}
            for row in rows:
                if len(row) <= max(x_idx, y_idx):
                    continue
                x = self._safe_float(row[x_idx])
                y = self._safe_float(row[y_idx])
                if x is not None and y is not None:
                    counts_by_insert_size[x] = counts_by_insert_size.get(x, 0.0) + y

            total_counts = sum(counts_by_insert_size.values())
            series_percent = {}
            if total_counts > 0:
                for insert_size, ct in counts_by_insert_size.items():
                    if 0 <= insert_size <= 1000:
                        series_percent[insert_size] = 100.0 * ct / total_counts
            sample = self._sample_name_from_file(f)
            if series_percent:
                data_percent[sample] = series_percent

        if not data_percent:
            return

        self.add_section(
            name="Insert size distribution",
            description="QoRTs insert size distribution.",
            plot=linegraph.plot(
                data_percent,
                {
                    "id": "qorts_insert_size_distribution",
                    "title": "QoRTs: Insert Size Distribution",
                    "xlab": "Insert size",
                    "ylab": "Normalized read count",
                    "ysuffix": "%",
                    "xmin": 0,
                    "xmax": 1000,
                    "ymin": 0,
                },
            ),
        )

    def qorts_read_length_distribution_plot(self):
        data_counts = {}
        for f in self._iter_qorts_files("qorts/read_len_dist"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue

            len_idx = self._find_column(headers, [r"\blen\b", r"length", r"size"])
            r1_idx = self._find_column(headers, [r"ct_r1", r"count_r1", r"\br1\b"])
            r2_idx = self._find_column(headers, [r"ct_r2", r"count_r2", r"\br2\b"])
            ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
            if len_idx is None:
                continue

            sample = self._sample_name_from_file(f)
            if r1_idx is not None:
                series_r1 = self._aggregate_by_xy(rows, len_idx, r1_idx)
                if series_r1:
                    data_counts[f"{sample} R1"] = series_r1
            if r2_idx is not None:
                series_r2 = self._aggregate_by_xy(rows, len_idx, r2_idx)
                if series_r2:
                    data_counts[f"{sample} R2"] = series_r2
            if r1_idx is None and r2_idx is None and ct_idx is not None:
                series = self._aggregate_by_xy(rows, len_idx, ct_idx)
                if series:
                    data_counts[sample] = series

        if not data_counts:
            return

        data_percent = {}
        for sample, series in data_counts.items():
            total = sum(series.values())
            if total > 0:
                data_percent[sample] = {x: 100.0 * y / total for x, y in series.items()}

        data_sets = [data_counts]
        data_labels = [{"name": "Count", "ylab": "Read count"}]
        if data_percent:
            data_sets.append(data_percent)
            data_labels.append({"name": "Percent", "ylab": "Read fraction (%)"})

        self.add_section(
            name="Read length distribution",
            description="QoRTs read length histogram.",
            plot=linegraph.plot(
                data_sets,
                {
                    "id": "qorts_read_length_distribution",
                    "title": "QoRTs: Read Length Distribution",
                    "xlab": "Read length",
                    "ylab": "Read count",
                    "xmin": 0,
                    "ymin": 0,
                    "data_labels": data_labels,
                },
            ),
        )

    def qorts_clipping_distribution_plots(self):
        clip_specs = [
            (
                "qorts/nvc_lead_clip_r1",
                "leadClipLen",
                "Leading clip distribution (R1)",
                "qorts_lead_clip_r1",
                "QoRTs: Leading Clip Length Distribution (R1)",
            ),
            (
                "qorts/nvc_lead_clip_r2",
                "leadClipLen",
                "Leading clip distribution (R2)",
                "qorts_lead_clip_r2",
                "QoRTs: Leading Clip Length Distribution (R2)",
            ),
            (
                "qorts/nvc_tail_clip_r1",
                "tailClipLen",
                "Trailing clip distribution (R1)",
                "qorts_tail_clip_r1",
                "QoRTs: Trailing Clip Length Distribution (R1)",
            ),
            (
                "qorts/nvc_tail_clip_r2",
                "tailClipLen",
                "Trailing clip distribution (R2)",
                "qorts_tail_clip_r2",
                "QoRTs: Trailing Clip Length Distribution (R2)",
            ),
        ]

        for key, clip_col_name, section_name, plot_id, title in clip_specs:
            data_counts = {}
            for f in self._iter_qorts_files(key):
                headers, rows = self._parse_tsv_rows(f)
                if not headers:
                    continue

                clip_idx = self._find_column(headers, [clip_col_name.lower(), r"cliplen", r"clip"])
                ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
                if clip_idx is None or ct_idx is None:
                    continue

                series = self._aggregate_by_xy(rows, clip_idx, ct_idx)
                if series:
                    data_counts[self._sample_name_from_file(f)] = series

            if not data_counts:
                continue

            data_percent = {}
            for sample, series in data_counts.items():
                total = sum(series.values())
                if total > 0:
                    data_percent[sample] = {x: 100.0 * y / total for x, y in series.items()}

            if plot_id in {
                "qorts_lead_clip_r1",
                "qorts_lead_clip_r2",
                "qorts_tail_clip_r1",
                "qorts_tail_clip_r2",
            }:
                if not data_percent:
                    continue
                data_sets = data_percent
                p_ylab = "Read fraction"
                p_ysuffix = "%"
                p_data_labels = None
            else:
                data_sets = [data_counts]
                p_data_labels = [{"name": "Count", "ylab": "Read count"}]
                if data_percent:
                    data_sets.append(data_percent)
                    p_data_labels.append({"name": "Percent", "ylab": "Read fraction (%)"})
                p_ylab = "Read count"
                p_ysuffix = None

            self.add_section(
                name=section_name,
                description="QoRTs clipping-length distribution.",
                plot=linegraph.plot(
                    data_sets,
                    {
                        "id": plot_id,
                        "title": title,
                        "xlab": "Clipped bases",
                        "ylab": p_ylab,
                        "ysuffix": p_ysuffix,
                        "xmin": 0,
                        "ymin": 0,
                        "data_labels": p_data_labels,
                    },
                ),
            )

    def qorts_nvc_clip_match_by_clip_position_plot(self):
        data = {}
        for key in ["qorts/nvc_minus_clipping_r1", "qorts/nvc_minus_clipping_r2"]:
            for f in self._iter_qorts_files(key):
                headers, rows = self._parse_tsv_rows(f)
                if not headers:
                    continue
                pos_idx = self._find_column(headers, [r"readpos", r"cycle", r"position"])
                base_idx = self._find_column(headers, [r"base"])
                ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
                if pos_idx is None or base_idx is None or ct_idx is None:
                    continue

                by_pos: Dict[float, Dict[str, float]] = {}
                for row in rows:
                    if len(row) <= max(pos_idx, base_idx, ct_idx):
                        continue
                    pos = self._safe_float(row[pos_idx])
                    ct = self._safe_float(row[ct_idx])
                    base = row[base_idx].upper()
                    if pos is None or ct is None:
                        continue
                    by_pos.setdefault(pos, {})
                    by_pos[pos][base] = by_pos[pos].get(base, 0.0) + ct

                series = {}
                for pos, counts in by_pos.items():
                    total = sum(counts.values())
                    if total > 0:
                        series[pos] = 100.0 * (
                            counts.get("A", 0.0)
                            + counts.get("C", 0.0)
                            + counts.get("G", 0.0)
                            + counts.get("T", 0.0)
                        ) / total

                if series:
                    suffix = " R1" if key.endswith("r1") else " R2"
                    data[f"{self._sample_name_from_file(f)}{suffix}"] = series

        if not data:
            return

        self.add_section(
            name="NVC clip match by clip position",
            description="QoRTs clip-position match profile derived from minus-clipping NVC tables.",
            plot=linegraph.plot(
                data,
                {
                    "id": "qorts_nvc_clip_match_by_clip_position",
                    "title": "QoRTs: NVC Clip Match by Clip Position",
                    "xlab": "Read position",
                    "ylab": "Match rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                    "ymax": 100,
                },
            ),
        )

    def qorts_gene_body_coverage_plot(self):
        overall = {}
        upper_mid = {}
        low_expr = {}

        for f in self._iter_qorts_files("qorts/gene_body_coverage"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            q_idx = self._find_column(headers, [r"quantile"])
            low_idx = self._find_column(headers, [r"(?:x?1\.)?bottomhalf", r"bottomhalf", r"low"])
            um_idx = self._find_column(headers, [r"(?:x?2\.)?uppermidquartile", r"uppermidquartile", r"upper.*mid"])
            high1_idx = self._find_column(headers, [r"(?:x?3\.)?75to90", r"75to90"])
            high2_idx = self._find_column(headers, [r"(?:x?4\.)?high", r"high"])

            if q_idx is None or low_idx is None or um_idx is None:
                continue

            sample = self._sample_name_from_file(f)
            d_overall = {}
            d_um = {}
            d_low = {}
            for row in rows:
                if len(row) <= max(q_idx, low_idx, um_idx):
                    continue
                q = self._safe_float(row[q_idx])
                low = self._safe_float(row[low_idx])
                um = self._safe_float(row[um_idx])
                h1 = self._safe_float(row[high1_idx]) if high1_idx is not None and len(row) > high1_idx else 0.0
                h2 = self._safe_float(row[high2_idx]) if high2_idx is not None and len(row) > high2_idx else 0.0
                if q is None or low is None or um is None:
                    continue

                # QoRTs quantiles may be 0-1 fractions; normalize to 0-100 scale.
                q_scaled = q * 100.0 if q <= 1.2 else q
                d_overall[q_scaled] = low + um + (h1 or 0.0) + (h2 or 0.0)
                d_um[q_scaled] = um
                d_low[q_scaled] = low

            # Convert raw coverage counts into relative per-bin rates.
            total_overall = sum(d_overall.values())
            total_um = sum(d_um.values())
            total_low = sum(d_low.values())
            if total_overall > 0:
                d_overall = {qv: 100.0 * v / total_overall for qv, v in d_overall.items()}
            if total_um > 0:
                d_um = {qv: 100.0 * v / total_um for qv, v in d_um.items()}
            if total_low > 0:
                d_low = {qv: 100.0 * v / total_low for qv, v in d_low.items()}

            if d_overall:
                overall[sample] = d_overall
                upper_mid[sample] = d_um
                low_expr[sample] = d_low

        if not overall:
            return

        self.add_section(
            name="Gene body coverage (overall)",
            description="QoRTs gene body coverage aggregated across all expression bins: <50 (low), 50-75 (upper-middle), 75-90 (higher), and >90 (highest expression).",
            plot=linegraph.plot(
                overall,
                {
                    "id": "qorts_gene_body_cov_overall",
                    "title": "QoRTs: Gene Body Coverage (Overall)",
                    "xlab": "Gene-body quantile",
                    "ylab": "Relative coverage",
                    "ysuffix": "%",
                    "xmin": 0,
                    "xmax": 100,
                },
            ),
        )
        self.add_section(
            name="Gene body coverage (upper-middle quartile)",
            description="QoRTs gene body coverage for the 50-75 expression bin only (upper-middle expression quartile).",
            plot=linegraph.plot(
                upper_mid,
                {
                    "id": "qorts_gene_body_cov_upper_mid",
                    "title": "QoRTs: Gene Body Coverage (Upper-Middle Quartile)",
                    "xlab": "Gene-body quantile",
                    "ylab": "Relative coverage",
                    "ysuffix": "%",
                    "xmin": 0,
                    "xmax": 100,
                },
            ),
        )
        self.add_section(
            name="Gene body coverage (low expression)",
            description="QoRTs gene body coverage for the <50 expression bin only (low-expression genes).",
            plot=linegraph.plot(
                low_expr,
                {
                    "id": "qorts_gene_body_cov_low",
                    "title": "QoRTs: Gene Body Coverage (Low Expression)",
                    "xlab": "Gene-body quantile",
                    "ylab": "Relative coverage",
                    "ysuffix": "%",
                    "xmin": 0,
                    "xmax": 100,
                },
            ),
        )

    def qorts_biotype_rates_plot(self):
        data = {}
        cats: Dict[str, Dict[str, str]] = {}
        for f in self._iter_qorts_files("qorts/biotype_counts"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            b_idx = self._find_column(headers, [r"biotype", r"type", r"category"])
            c_idx = self._find_column(headers, [r"\bct\b", r"count"])
            if b_idx is None:
                b_idx = 0
            if c_idx is None:
                continue

            sample_data = {}
            for row in rows:
                if len(row) <= max(b_idx, c_idx):
                    continue
                biotype = row[b_idx]
                ct = self._safe_float(row[c_idx])
                if ct is None:
                    continue
                sample_data[biotype] = sample_data.get(biotype, 0.0) + ct
                cats[biotype] = {"name": biotype}

            if sample_data:
                data[self._sample_name_from_file(f)] = sample_data

        if not data:
            return

        self.add_section(
            name="Biotype rates",
            description="QoRTs biotype count/rate distribution.",
            plot=bargraph.plot(
                data,
                cats,
                {
                    "id": "qorts_biotype_rates",
                    "title": "QoRTs: Biotype Rates",
                    "ylab": "Count",
                    "cpswitch_counts_label": "Count",
                    "cpswitch_percent_label": "Percent",
                    "cpswitch_c_active": False,
                    "y_decimals": 0,
                    "tt_decimals": 0,
                    "hide_zero_cats": False,
                },
            ),
        )

    def qorts_chrom_type_rates_plot(self):
        data: Dict[str, Dict[str, float]] = {}
        for f in self._iter_qorts_files("qorts/chrom_counts"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            chrom_idx = self._find_column(headers, [r"chrom", r"chr", r"type", r"category"])
            ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
            if chrom_idx is None:
                chrom_idx = 0
            if ct_idx is None:
                continue

            sample_data = {}
            for row in rows:
                if len(row) <= max(chrom_idx, ct_idx):
                    continue
                chrom = row[chrom_idx]
                ct = self._safe_float(row[ct_idx])
                if ct is None:
                    continue
                sample_data[chrom] = sample_data.get(chrom, 0.0) + ct

            if sample_data:
                data[self._sample_name_from_file(f)] = sample_data

        if not data:
            return

        chrom_totals: Dict[str, float] = {}
        for sample_data in data.values():
            for chrom, ct in sample_data.items():
                chrom_totals[chrom] = chrom_totals.get(chrom, 0.0) + ct

        total_cov = sum(chrom_totals.values())

        def _chrom_core(chrom: str) -> str:
            c = chrom.strip().lower()
            if c.startswith("chr"):
                c = c[3:]
            return c

        def _roman_to_int(s: str) -> Optional[int]:
            vals = {"i": 1, "v": 5, "x": 10, "l": 50, "c": 100, "d": 500, "m": 1000}
            s = s.lower()
            if not s or any(ch not in vals for ch in s):
                return None
            total = 0
            prev = 0
            for ch in reversed(s):
                v = vals[ch]
                if v < prev:
                    total -= v
                else:
                    total += v
                    prev = v
            return total

        def _natural_key(text: str):
            return [int(tok) if tok.isdigit() else tok for tok in re.split(r"(\d+)", text.lower())]

        def _primary_chrom_rank(chrom: str) -> Optional[tuple]:
            c = _chrom_core(chrom)

            # Numeric chromosomes, including forms like 1a / 2b seen in some assemblies.
            m_num = re.match(r"^(\d+)([a-z]*)$", c)
            if m_num:
                n = int(m_num.group(1))
                suffix = m_num.group(2)
                return (0, n, suffix)

            # Roman numeral chromosomes (common in some non-human references).
            roman_val = _roman_to_int(c)
            if roman_val is not None:
                return (1, roman_val, "")

            # Common sex / mito chromosomes across species.
            if c in {"x", "y", "z", "w"}:
                sex_order = {"x": 0, "y": 1, "z": 2, "w": 3}
                return (2, sex_order[c], "")
            if c in {"m", "mt", "mito", "mitochondria"}:
                return (3, 0, "")

            return None

        def _is_canonical(chrom: str) -> bool:
            if _primary_chrom_rank(chrom) is not None:
                return True

            # Keep explicitly named top-level chromosomes such as chrA / chrB.
            c = _chrom_core(chrom)
            if re.match(r"^[a-z]$", c):
                return True

            # Exclude common auxiliary/decoy/unlocalized patterns.
            aux_keywords = ["random", "_alt", "_fix", "_hap", "un_", "decoy", "patch"]
            if any(k in c for k in aux_keywords):
                return False

            # Names with no separators are often primary in many assemblies.
            if re.match(r"^[a-z0-9]+$", c):
                return True

            return False

            if c.isdigit():
                return True
            return _primary_chrom_rank(chrom) is not None

        # Keep canonical chromosomes; for non-canonical contigs keep only those with meaningful coverage.
        min_noncanonical_fraction = float(
            getattr(config, "qorts", {}).get("chrom_min_noncanonical_fraction", 0.001)
        )
        kept_chroms = []
        for chrom, ct in chrom_totals.items():
            frac = (ct / total_cov) if total_cov > 0 else 0.0
            if _is_canonical(chrom) or frac >= min_noncanonical_fraction:
                kept_chroms.append(chrom)

        def _chrom_sort_key(chrom: str):
            primary_rank = _primary_chrom_rank(chrom)
            if primary_rank is not None:
                return (0, *primary_rank, "")
            if _is_canonical(chrom):
                return (1, 0, "", _natural_key(chrom))
            return (2, 0, "", -chrom_totals.get(chrom, 0.0), _natural_key(chrom))

        kept_chroms = sorted(kept_chroms, key=_chrom_sort_key)

        filtered_data: Dict[str, Dict[str, float]] = {}
        for sample, sample_data in data.items():
            row = {chrom: sample_data[chrom] for chrom in kept_chroms if chrom in sample_data}
            if row:
                filtered_data[sample] = row

        if not filtered_data:
            return

        chrom_palette = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
            "#393b79",
            "#637939",
            "#8c6d31",
            "#843c39",
            "#7b4173",
            "#3182bd",
            "#31a354",
            "#756bb1",
            "#636363",
            "#e6550d",
            "#6baed6",
            "#74c476",
            "#9e9ac8",
            "#969696",
            "#fd8d3c",
            "#9ecae1",
            "#a1d99b",
            "#bcbddc",
            "#bdbdbd",
            "#fdae6b",
            "#c6dbef",
            "#c7e9c0",
        ]

        def _chrom_color(i: int) -> str:
            if i < len(chrom_palette):
                return chrom_palette[i]
            # Golden-angle hue stepping gives stable, visually separated fallback colors.
            hue = (i * 0.618033988749895) % 1.0
            r, g, b = colorsys.hsv_to_rgb(hue, 0.55, 0.95)
            return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"

        cats: Dict[str, Dict[str, str]] = {
            chrom: {"name": chrom, "color": _chrom_color(i)} for i, chrom in enumerate(kept_chroms)
        }

        self.add_section(
            name="Chrom type rates",
            description="QoRTs chromosome count/rate distribution. Low-coverage non-canonical contigs are hidden.",
            plot=bargraph.plot(
                filtered_data,
                cats,
                {
                    "id": "qorts_chrom_type_rates",
                    "title": "QoRTs: Chrom Type Rates",
                    "ylab": "Count",
                    "cpswitch_counts_label": "Count",
                    "cpswitch_percent_label": "Percent",
                    "y_decimals": 0,
                    "tt_decimals": 0,
                    "hide_zero_cats": False,
                },
            ),
        )

    def qorts_cigar_length_distribution_plot(self):
        op_data: Dict[str, Dict[str, Dict[float, float]]] = {}
        for key in ["qorts/cigar_op_lengths_r1", "qorts/cigar_op_lengths_r2"]:
            for f in self._iter_qorts_files(key):
                headers, rows = self._parse_tsv_rows(f)
                if not headers:
                    continue
                op_idx = self._find_column(headers, [r"op", r"cigar"])
                len_idx = self._find_column(headers, [r"len", r"length", r"size"])
                ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
                if op_idx is None or len_idx is None or ct_idx is None:
                    continue

                label = f"{self._sample_name_from_file(f)} {'R1' if key.endswith('r1') else 'R2'}"
                tmp: Dict[str, Dict[float, float]] = {}
                for row in rows:
                    if len(row) <= max(op_idx, len_idx, ct_idx):
                        continue
                    op = row[op_idx]
                    length = self._safe_float(row[len_idx])
                    ct = self._safe_float(row[ct_idx])
                    if length is None or ct is None:
                        continue
                    tmp.setdefault(op, {})
                    tmp[op][length] = tmp[op].get(length, 0.0) + ct

                for op, series in tmp.items():
                    series_pct = self._series_to_percent(series)
                    if series_pct:
                        op_data.setdefault(op, {})[label] = series_pct

        if not op_data:
            return

        data_sets = []
        labels = []
        for op, d in sorted(op_data.items()):
            data_sets.append(d)
            labels.append({"name": f"Op {op}"})

        self.add_section(
            name="CIGAR length distribution",
            description="QoRTs CIGAR operation length distributions (normalized to percent within each sample and op).",
            plot=linegraph.plot(
                data_sets,
                {
                    "id": "qorts_cigar_length_distribution",
                    "title": "QoRTs: CIGAR Length Distribution",
                    "xlab": "Operation length",
                    "ylab": "Distribution",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                    "data_labels": labels,
                },
            ),
        )

    def qorts_cigar_op_by_cycle_plot(self):
        op_data: Dict[str, Dict[str, Dict[float, float]]] = {}
        for key in ["qorts/cigar_op_cycle_r1", "qorts/cigar_op_cycle_r2"]:
            for f in self._iter_qorts_files(key):
                headers, rows = self._parse_tsv_rows(f)
                if not headers:
                    continue
                pos_idx = self._find_column(headers, [r"^cycle$", r"readpos", r"position", r"base"])
                if pos_idx is None:
                    continue

                op_cols = []
                for i, h in enumerate(headers):
                    if i == pos_idx:
                        continue
                    m = re.match(r"^([A-Z])_[SMEB]$", h)
                    if m:
                        op_cols.append((i, m.group(1)))
                if not op_cols:
                    continue

                label = f"{self._sample_name_from_file(f)} {'R1' if key.endswith('r1') else 'R2'}"
                tmp: Dict[str, Dict[float, float]] = {}
                for row in rows:
                    if len(row) <= pos_idx:
                        continue
                    pos = self._safe_float(row[pos_idx])
                    if pos is None:
                        continue
                    for col_idx, op in op_cols:
                        if len(row) <= col_idx:
                            continue
                        ct = self._safe_float(row[col_idx])
                        if ct is None:
                            continue
                        tmp.setdefault(op, {})
                        tmp[op][pos] = tmp[op].get(pos, 0.0) + ct

                for op, series in tmp.items():
                    series_pct = self._series_to_percent(series)
                    if series_pct:
                        op_data.setdefault(op, {})[label] = series_pct

        if not op_data:
            return

        data_sets = []
        labels = []
        for op, d in sorted(op_data.items()):
            data_sets.append(d)
            labels.append({"name": f"Op {op}"})

        self.add_section(
            name="CIGAR op by cycle",
            description="QoRTs CIGAR operation profiles across read cycles (normalized to percent within each sample and op).",
            plot=linegraph.plot(
                data_sets,
                {
                    "id": "qorts_cigar_op_by_cycle",
                    "title": "QoRTs: CIGAR Operation by Cycle",
                    "xlab": "Read position",
                    "ylab": "Distribution",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                    "data_labels": labels,
                },
            ),
        )

    def _cigar_single_op_by_cycle(self, op_letter: str) -> Dict[str, Dict[float, float]]:
        data: Dict[str, Dict[float, float]] = {}
        for key in ["qorts/cigar_op_cycle_r1", "qorts/cigar_op_cycle_r2"]:
            for f in self._iter_qorts_files(key):
                headers, rows = self._parse_tsv_rows(f)
                if not headers:
                    continue
                pos_idx = self._find_column(headers, [r"^cycle$", r"readpos", r"position", r"base"])
                if pos_idx is None:
                    continue

                op_idx = None
                for i, h in enumerate(headers):
                    if i == pos_idx:
                        continue
                    if re.match(rf"^{re.escape(op_letter)}_[SMEB]$", h):
                        op_idx = i
                        break
                if op_idx is None:
                    continue

                series = self._aggregate_by_xy(rows, pos_idx, op_idx)
                if series:
                    label = f"{self._sample_name_from_file(f)} {'R1' if key.endswith('r1') else 'R2'}"
                    series_pct = self._series_to_percent(series)
                    if series_pct:
                        data[label] = series_pct
        return data

    def qorts_deletion_profile_by_cycle_plot(self):
        data = self._cigar_single_op_by_cycle("D")
        if not data:
            return

        self.add_section(
            name="Deletion profile by read position",
            description="QoRTs deletion CIGAR operation profile by read position (normalized to percent).",
            plot=linegraph.plot(
                data,
                {
                    "id": "qorts_deletion_profile_by_cycle",
                    "title": "QoRTs: Deletion Profile by Read Position",
                    "xlab": "Read position",
                    "ylab": "Rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                },
            ),
        )

    def qorts_insertion_profile_by_cycle_plot(self):
        data = self._cigar_single_op_by_cycle("I")
        if not data:
            return

        self.add_section(
            name="Insertion profile by read position",
            description="QoRTs insertion CIGAR operation profile by read position (normalized to percent).",
            plot=linegraph.plot(
                data,
                {
                    "id": "qorts_insertion_profile_by_cycle",
                    "title": "QoRTs: Insertion Profile by Read Position",
                    "xlab": "Read position",
                    "ylab": "Rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                },
            ),
        )

    def qorts_splicing_profile_by_cycle_plot(self):
        data = self._cigar_single_op_by_cycle("N")
        if not data:
            return

        self.add_section(
            name="Splicing profile by read position",
            description="QoRTs splice-junction CIGAR operation profile by read position (normalized to percent).",
            plot=linegraph.plot(
                data,
                {
                    "id": "qorts_splicing_profile_by_cycle",
                    "title": "QoRTs: Splicing Profile by Read Position",
                    "xlab": "Read position",
                    "ylab": "Rate",
                    "ysuffix": "%",
                    "xmin": 0,
                    "ymin": 0,
                },
            ),
        )

    def qorts_gene_cdf_plot(self):
        data = {}
        for f in self._iter_qorts_files("qorts/gene_counts"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            gene_idx = self._find_column(headers, [r"gene", r"geneid", r"id"])
            ct_idx = self._find_column(headers, [r"\bct\b", r"count"])
            if gene_idx is None:
                gene_idx = 0
            if ct_idx is None:
                continue

            counts = []
            for row in rows:
                if len(row) <= max(gene_idx, ct_idx):
                    continue
                gene = row[gene_idx]
                ct = self._safe_float(row[ct_idx])
                if ct is None:
                    continue
                if gene.startswith("_"):
                    continue
                counts.append(ct)

            counts.sort(reverse=True)
            if not counts:
                continue
            total = sum(counts)
            if total <= 0:
                continue
            series = {}
            csum = 0.0
            for i, ct in enumerate(counts, start=1):
                csum += ct
                series[float(i)] = 100.0 * csum / total
            data[self._sample_name_from_file(f)] = series

        if not data:
            return

        self.add_section(
            name="Gene CDF",
            description="QoRTs cumulative distribution of reads over ranked genes.",
            plot=linegraph.plot(
                data,
                {
                    "id": "qorts_gene_cdf",
                    "title": "QoRTs: Gene Count CDF",
                    "xlab": "Gene rank",
                    "ylab": "Cumulative fraction",
                    "ysuffix": "%",
                    "xlog": True,
                    "xmin": 0,
                    "ymin": 0,
                    "ymax": 100,
                },
            ),
        )

    def qorts_mapping_rates_plot(self):
        self._summary_bar_section(
            section_name="Mapping rates",
            section_desc="QoRTs mapping-rate related summary metrics.",
            section_id="qorts_mapping_rates",
            section_title="QoRTs: Mapping Rates",
            regexes=[r"mapping.*rate", r"mapped.*rate", r"total\.aligned\.rate", r"aligned\.and\.pf\.rate", r"mm\.rate"],
            scale_to_percent=True,
        )

    def qorts_dropped_rates_plot(self):
        self._summary_bar_section(
            section_name="Dropped rates",
            section_desc="QoRTs dropped-read related summary metrics.",
            section_id="qorts_dropped_rates",
            section_title="QoRTs: Dropped Rates",
            regexes=[r"^dropped", r"dropped_", r"drop.*rate", r"discard"],
            scale_to_percent=False,
            counts_label="Count",
            lowercase_labels=True,
        )

    def qorts_missingness_rate_plot(self):
        self._summary_bar_section(
            section_name="Missingness rate",
            section_desc="QoRTs missingness-related summary metrics.",
            section_id="qorts_missingness_rate",
            section_title="QoRTs: Missingness Rate",
            regexes=[r"missing", r"no[_\. ]?gene", r"middleofnowhere"],
            scale_to_percent=False,
        )

    def qorts_norm_factors_plot(self):
        self._summary_bar_section(
            section_name="Norm factors",
            section_desc="QoRTs normalization-factor summary metrics.",
            section_id="qorts_norm_factors",
            section_title="QoRTs: Norm Factors",
            regexes=[
                r"^norm[_\.].*",
                r".*[_\.]norm[_\.].*",
                r"^deseq.*norm.*",
                r"^edger.*norm.*",
                r"^normalization[_\.]factor.*",
            ],
        )

    def qorts_on_target_plot(self):
        data = {}
        for f in self._iter_qorts_files("qorts/on_target"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            start_idx = self._find_column(headers, [r"start"])
            end_idx = self._find_column(headers, [r"end"])
            read_cov_idx = self._find_column(headers, [r"read_?coverage", r"readcoverage"])
            pair_cov_idx = self._find_column(headers, [r"readpair_?coverage", r"paircoverage"])
            if start_idx is None or end_idx is None or (read_cov_idx is None and pair_cov_idx is None):
                continue

            depths = []
            for row in rows:
                if len(row) <= max(start_idx, end_idx):
                    continue
                start = self._safe_float(row[start_idx])
                end = self._safe_float(row[end_idx])
                if start is None or end is None:
                    continue
                span = end - start
                if span <= 0:
                    continue
                cov = None
                if read_cov_idx is not None and len(row) > read_cov_idx:
                    cov = self._safe_float(row[read_cov_idx])
                elif pair_cov_idx is not None and len(row) > pair_cov_idx:
                    cov = self._safe_float(row[pair_cov_idx])
                if cov is None:
                    continue
                depths.append(cov / span)

            depths.sort()
            if not depths:
                continue
            n = len(depths)
            series = {100.0 * (i + 1) / n: depth for i, depth in enumerate(depths)}
            data[self._sample_name_from_file(f)] = series

        if not data:
            return

        self.add_section(
            name="On-target",
            description="QoRTs on-target depth quantile profile.",
            plot=linegraph.plot(
                data,
                {
                    "id": "qorts_on_target",
                    "title": "QoRTs: On-Target Depth Profile",
                    "xlab": "Target depth quantile",
                    "ylab": "Depth",
                    "xmin": 0,
                    "xmax": 100,
                    "ymin": 0,
                },
            ),
        )

    def qorts_overlap_plot(self):
        data_cov = {}
        for f in self._iter_qorts_files("qorts/overlap_coverage"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            pos_idx = self._find_column(headers, [r"\bpos\b", r"readpos", r"cycle"])
            ct_r1_idx = self._find_column(headers, [r"ct_r1", r"count_r1", r"r1"])
            ct_r2_idx = self._find_column(headers, [r"ct_r2", r"count_r2", r"r2"])
            if pos_idx is None:
                continue

            sample = self._sample_name_from_file(f)
            if ct_r1_idx is not None:
                series = self._aggregate_by_xy(rows, pos_idx, ct_r1_idx)
                if series:
                    data_cov[f"{sample} R1"] = series
            if ct_r2_idx is not None:
                series = self._aggregate_by_xy(rows, pos_idx, ct_r2_idx)
                if series:
                    data_cov[f"{sample} R2"] = series

        if data_cov:
            data_cov_pct = {sample: self._series_to_percent(series) for sample, series in data_cov.items()}
            data_cov_pct = {sample: series for sample, series in data_cov_pct.items() if series}
        else:
            data_cov_pct = {}

        if data_cov_pct:
            self.add_section(
                name="Overlap coverage",
                description="QoRTs overlap coverage across read positions (normalized to percent).",
                plot=linegraph.plot(
                    data_cov_pct,
                    {
                        "id": "qorts_overlap_coverage",
                        "title": "QoRTs: Overlap Coverage",
                        "xlab": "Read position",
                        "ylab": "Coverage rate",
                        "ysuffix": "%",
                        "xmin": 0,
                        "ymin": 0,
                    },
                ),
            )

        data_mm = {}
        for f in self._iter_qorts_files("qorts/overlap_mismatch_by_read"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            pos_idx = self._find_column(headers, [r"\bpos\b", r"readpos", r"cycle"])
            mm_r1_idx = self._find_column(headers, [r"ct_r1", r"mismatch.*r1"])
            mm_r2_idx = self._find_column(headers, [r"ct_r2", r"mismatch.*r2"])
            if pos_idx is None:
                continue

            sample = self._sample_name_from_file(f)
            if mm_r1_idx is not None:
                series = self._aggregate_by_xy(rows, pos_idx, mm_r1_idx)
                if series:
                    data_mm[f"{sample} R1"] = series
            if mm_r2_idx is not None:
                series = self._aggregate_by_xy(rows, pos_idx, mm_r2_idx)
                if series:
                    data_mm[f"{sample} R2"] = series

        if data_mm:
            data_mm_pct = {sample: self._series_to_percent(series) for sample, series in data_mm.items()}
            data_mm_pct = {sample: series for sample, series in data_mm_pct.items() if series}
        else:
            data_mm_pct = {}

        if data_mm_pct:
            self.add_section(
                name="Overlap mismatch",
                description="QoRTs overlap mismatch profile across read positions (normalized to percent).",
                plot=linegraph.plot(
                    data_mm_pct,
                    {
                        "id": "qorts_overlap_mismatch",
                        "title": "QoRTs: Overlap Mismatch by Read Position",
                        "xlab": "Read position",
                        "ylab": "Mismatch rate",
                        "ysuffix": "%",
                        "xmin": 0,
                        "ymin": 0,
                    },
                ),
            )

    def qorts_reference_mismatch_family_plot(self):
        data_counts = {}
        for f in self._iter_qorts_files("qorts/reference_mismatch_counts"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            ref_idx = self._find_column(headers, [r"ref", r"refbase"])
            read_idx = self._find_column(headers, [r"read", r"readbase"])
            ct_idx = self._find_column(headers, [r"ct_r1", r"\bct\b", r"count"])
            if ref_idx is None or read_idx is None or ct_idx is None:
                continue

            sample_data = {}
            for row in rows:
                if len(row) <= max(ref_idx, read_idx, ct_idx):
                    continue
                k = f"{row[ref_idx]}->{row[read_idx]}"
                ct = self._safe_float(row[ct_idx])
                if ct is None:
                    continue
                sample_data[k] = sample_data.get(k, 0.0) + ct
            if sample_data:
                data_counts[self._sample_name_from_file(f)] = sample_data

        if data_counts:
            cats = {k: {"name": k} for k in sorted({x for s in data_counts.values() for x in s.keys()})}
            self.add_section(
                name="Reference mismatch counts",
                description="QoRTs reference mismatch substitution counts.",
                plot=bargraph.plot(
                    data_counts,
                    cats,
                    {
                        "id": "qorts_reference_mismatch_counts",
                        "title": "QoRTs: Reference Mismatch Counts",
                        "ylab": "Count",
                        "cpswitch_counts_label": "Count",
                        "cpswitch_percent_label": "Percent",
                        "hide_zero_cats": False,
                    },
                ),
            )

        data_score = {}
        for f in self._iter_qorts_files("qorts/reference_mismatch_by_score"):
            headers, rows = self._parse_tsv_rows(f)
            if not headers:
                continue
            score_idx = self._find_column(headers, [r"score", r"qual"])
            mm_idx = self._find_column(headers, [r"mismatch", r"ct", r"count"])
            cov_idx = self._find_column(headers, [r"coverage", r"qualcoverage"])
            if score_idx is None or mm_idx is None:
                continue

            series = {}
            for row in rows:
                if len(row) <= max(score_idx, mm_idx):
                    continue
                score = self._safe_float(row[score_idx])
                mm = self._safe_float(row[mm_idx])
                if score is None or mm is None:
                    continue
                if cov_idx is not None and len(row) > cov_idx:
                    cov = self._safe_float(row[cov_idx])
                    if cov and cov > 0:
                        series[score] = 100.0 * mm / cov
                    else:
                        series[score] = mm
                else:
                    series[score] = mm
            if series:
                data_score[self._sample_name_from_file(f)] = series

        if data_score:
            self.add_section(
                name="Reference mismatch by score",
                description="QoRTs mismatch profile by base quality score.",
                plot=linegraph.plot(
                    data_score,
                    {
                        "id": "qorts_reference_mismatch_by_score",
                        "title": "QoRTs: Reference Mismatch by Score",
                        "xlab": "Quality score",
                        "ylab": "Mismatch rate",
                        "ysuffix": "%",
                        "xmin": 0,
                        "ymin": 0,
                    },
                ),
            )

    def qorts_runtime_performance_plot(self):
        self._summary_bar_section(
            section_name="Run time performance",
            section_desc="QoRTs runtime benchmark metrics.",
            section_id="qorts_runtime_performance",
            section_title="QoRTs: Runtime Performance",
            regexes=[r"^benchmark", r"runtime", r"minutes", r"seconds", r"time"],
        )
