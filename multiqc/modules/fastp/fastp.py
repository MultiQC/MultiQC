import json
import logging
import re
from collections import Counter, defaultdict
from typing import Any, Dict, Optional, Tuple

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table
from multiqc.types import ColumnKey

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    By default, the module generates the sample names based on the `--report_title` / `-R`
    option in the fastp command line (if present), or the input FastQ file names if not.

    If you prefer, you can tell the module to use the filenames as sample names instead.
    To do so, use the following config option:

    ```yaml
    use_filename_as_sample_name:
        - fastp
    ```

    See [Using log filenames as sample names](https://docs.seqera.io/multiqc/getting_started/config#using-log-filenames-as-sample-names)
    for more details.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="fastp",
            anchor="fastp",
            href="https://github.com/OpenGene/fastp",
            info="All-in-one FASTQ preprocessor (QC, adapters, trimming, filtering, splitting...)",
            extra="""
            Fastp goes through fastq files in a folder and perform a series of quality control and filtering.
            Quality control and reporting are displayed both before and after filtering, allowing for a clear
            depiction of the consequences of the filtering process. Notably, the latter can be conducted on a
            variety of parameters including quality scores, length, as well as the presence of adapters, polyG,
            or polyX tailing.""",
            doi="10.1093/bioinformatics/bty560",
        )

        data_by_sample = dict()
        for f in self.find_log_files("fastp", filehandles=True):
            s_name, parsed_json = self.parse_fastp_log(f)
            if not s_name:
                continue
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            data_by_sample[s_name] = parsed_json

        # Filter to strip out ignored sample names
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")

        # Find and load any fastp reports
        self.fastp_data = dict()
        self.fastp_duplication_plotdata = dict()
        self.fastp_insert_size_data = dict()
        self.fastp_all_data = dict()
        self.fastp_qual_plotdata = dict()
        self.fastp_gc_content_data = dict()
        self.fastp_n_content_data = dict()
        self.fastp_overrepresented_sequences = dict()
        for k in [
            "read1_before_filtering",
            "read2_before_filtering",
            "read1_after_filtering",
            "read2_after_filtering",
            "merged_and_filtered",
        ]:
            self.fastp_qual_plotdata[k] = dict()
            self.fastp_gc_content_data[k] = dict()
            self.fastp_n_content_data[k] = dict()
            self.fastp_overrepresented_sequences[k] = dict()
        for s_name, parsed_json in data_by_sample.items():
            self.process_parsed_data(parsed_json, s_name)

        # Save entire original parsed JSON
        self.write_data_file(self.fastp_all_data, "multiqc_fastp")

        # General Stats Table
        self.fastp_general_stats_table(self.fastp_data)

        # Filtering statistics bar plot
        self.add_section(
            name="Filtered Reads",
            anchor="fastp-filtered-reads-chart",
            description="Filtering statistics of sampled reads.",
            plot=self.fastp_filtered_reads_chart(),
        )

        # Duplication rate plot
        self.add_section(
            name="Duplication Rates",
            anchor="fastp-duprates",
            description="Duplication rates of sampled reads.",
            plot=linegraph.plot(
                self.fastp_duplication_plotdata,
                {
                    "id": "fastp-duprates-plot",
                    "title": "Fastp: Duplication Rate",
                    "xlab": "Duplication level",
                    "ylab": "Read percent",
                    "y_clipmax": 100,
                    "ymin": 0,
                    "tt_label": "{point.x}: {point.y:.2f}%",
                },
            ),
        )

        self.add_section(
            name="Insert Sizes",
            anchor="fastp-insert-size",
            description="Insert size estimation of sampled reads.",
            plot=linegraph.plot(
                self.fastp_insert_size_data,
                {
                    "id": "fastp-insert-size-plot",
                    "title": "Fastp: Insert Size Distribution",
                    "xlab": "Insert size",
                    "ylab": "Read percent",
                    "y_clipmax": 100,
                    "ymin": 0,
                    "tt_label": "{point.x}: {point.y:.2f}%",
                    "smooth_points": 300,
                    "smooth_points_sumcounts": False,
                },
            ),
        )

        # Base quality plot
        self.add_section(
            name="Sequence Quality",
            anchor="fastp-seq-quality",
            description="Average sequencing quality over each base of all reads.",
            plot=self.fastp_read_qual_plot(),
        )

        # GC content plot
        self.add_section(
            name="GC Content",
            anchor="fastp-seq-content-gc",
            description="Average GC content over each base of all reads.",
            plot=self.fastp_read_gc_plot(),
        )

        # N content plot
        self.add_section(
            name="N content",
            anchor="fastp-seq-content-n",
            description="Average N content over each base of all reads.",
            plot=self.fastp_read_n_plot(),
        )

        # Overrepresented sequences plot
        self.add_section(
            name="Overrepresented Sequences",
            anchor="fastp-overrepresented-sequences",
            description="Overrepresented sequences in the reads.",
            plot=self.fastp_overrepresented_sequences_plot(),
        )

    def parse_fastp_log(self, f) -> Tuple[Optional[str], Dict]:
        """Parse the JSON output from fastp and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse fastp JSON: '{f['fn']}': {e}, skipping sample")
            return None, {}
        if not isinstance(parsed_json, dict) or "command" not in parsed_json:
            log.warning(f"Could not find 'command' field in JSON: '{f['fn']}', skipping sample")
            return None, {}

        s_name = None
        # Check if we should use filename instead of parsed sample name
        should_use_filename = False
        if isinstance(config.use_filename_as_sample_name, list):
            # Check for module anchor
            if self.anchor in config.use_filename_as_sample_name:
                should_use_filename = True
        elif config.use_filename_as_sample_name is True:
            should_use_filename = True
        elif getattr(config, "fastp", {}).get("s_name_filenames", False):
            # Deprecated option - warn user
            log.warning(
                "The 'fastp.s_name_filenames' config option is deprecated. Use the global 'use_filename_as_sample_name' option instead."
            )
            should_use_filename = True

        if should_use_filename:
            s_name = f["s_name"]

        if s_name is None:
            # Parse the "command" line usually found in the JSON, and use the report title
            # if present, otherwise fall back to the input FastQ file names.
            cmd = parsed_json["command"].strip()

            # First, try to extract --report_title / -R if present
            # The value can contain spaces and should be everything until the next option or end
            m = re.search(r"(-R|--report_title)\s(.+?)(?:\s-|$)", cmd)
            if m:
                s_name = self.clean_s_name(m.group(2), f)
            else:
                # Fall back to input file names
                # On caveat is that the command won't have file names escaped properly,
                # so we need some special logic to account for names with spaces:
                # "fastp -c -g -y -i Sample 1 1.fastq.gz -o ..."
                # "fastp -c -g -y --in1 Sample 1 1.fastq.gz --out1 ..."
                # "fastp -c -g -y --in1 Sample 1 1.fastq.gz --in2 Sample 1_2.fastq.gz --out1 ..."
                #
                # Using a regex that extracts everything between "-i " or "--in1 " and " -".
                # It still won't work exactly right for file names with dashes following a
                # space, but that's a pretty rare case, and will still extract something
                # meaningful.
                s_names = []
                m = re.search(r"(-i|--in1)\s(.+?)(?:\s-|$)", cmd)
                if m:
                    s_names.append(m.group(2))
                    # Second input for paired end?
                    m = re.search(r"(-I|--in2)\s(.+?)(?:\s-|$)", cmd)
                    if m:
                        s_names.append(m.group(2))
                    s_name = self.clean_s_name(s_names, f)
                else:
                    s_name = f["s_name"]
                    log.warning(
                        f"Could not parse sample name from the fastp command:\n{cmd}\n"
                        f"Falling back to extracting it from the file name: "
                        f'"{f["fn"]}" -> "{s_name}"'
                    )

        self.add_data_source(f, s_name)
        return s_name, parsed_json

    def process_parsed_data(self, parsed_json: Dict, s_name: str):
        """Process the JSON extracted from logs"""

        self.fastp_data[s_name] = {}
        self.fastp_duplication_plotdata[s_name] = {}
        self.fastp_insert_size_data[s_name] = {}
        self.fastp_all_data[s_name] = parsed_json
        for k in [
            "read1_before_filtering",
            "read2_before_filtering",
            "read1_after_filtering",
            "read2_after_filtering",
            "merged_and_filtered",
        ]:
            self.fastp_qual_plotdata[k][s_name] = {}
            self.fastp_gc_content_data[k][s_name] = {}
            self.fastp_n_content_data[k][s_name] = {}

        # Parse filtering_result
        try:
            for k in parsed_json["filtering_result"]:
                self.fastp_data[s_name][f"filtering_result_{k}"] = float(parsed_json["filtering_result"][k])
        except KeyError:
            log.debug(f"fastp JSON did not have 'filtering_result' key: '{s_name}'")

        # Parse duplication
        try:
            self.fastp_data[s_name]["pct_duplication"] = float(parsed_json["duplication"]["rate"] * 100.0)
        except KeyError:
            log.debug(f"fastp JSON did not have a 'duplication' key: '{s_name}'")

        # Parse before_filtering
        try:
            for k in parsed_json["summary"]["before_filtering"]:
                self.fastp_data[s_name][f"before_filtering_{k}"] = float(parsed_json["summary"]["before_filtering"][k])
        except KeyError:
            log.debug(f"fastp JSON did not have a 'summary'-'before_filtering' keys: '{s_name}'")

        # Parse after_filtering
        try:
            for k in parsed_json["summary"]["after_filtering"]:
                self.fastp_data[s_name][f"after_filtering_{k}"] = float(parsed_json["summary"]["after_filtering"][k])
        except KeyError:
            log.debug(f"fastp JSON did not have a 'summary'-'after_filtering' keys: '{s_name}'")

        # Parse data required to calculate Pct reads surviving
        try:
            self.fastp_data[s_name]["before_filtering_total_reads"] = float(
                parsed_json["summary"]["before_filtering"]["total_reads"]
            )
        except KeyError:
            log.debug(f"Could not find pre-filtering # reads: '{s_name}'")

        try:
            self.fastp_data[s_name]["pct_surviving"] = (
                self.fastp_data[s_name]["filtering_result_passed_filter_reads"]
                / self.fastp_data[s_name]["before_filtering_total_reads"]
            ) * 100.0
        except (KeyError, ZeroDivisionError) as e:
            log.debug(f"Could not calculate 'pct_surviving' ({e.__class__.__name__}): {s_name}")

        # Parse adapter_cutting
        try:
            for k in parsed_json["adapter_cutting"]:
                try:
                    self.fastp_data[s_name][f"adapter_cutting_{k}"] = float(parsed_json["adapter_cutting"][k])
                except (ValueError, TypeError):
                    pass
        except KeyError:
            log.debug(f"fastp JSON did not have a 'adapter_cutting' key, skipping: '{s_name}'")

        try:
            self.fastp_data[s_name]["pct_adapter"] = (
                self.fastp_data[s_name]["adapter_cutting_adapter_trimmed_reads"]
                / self.fastp_data[s_name]["before_filtering_total_reads"]
            ) * 100.0
        except (KeyError, ZeroDivisionError) as e:
            log.debug(f"Could not calculate 'pct_adapter' ({e.__class__.__name__}): {s_name}")

        # Duplication rate plot data
        try:
            # First count the total read count in the dup analysis
            total_reads = 0
            for v in parsed_json["duplication"]["histogram"]:
                total_reads += v
            if total_reads == 0:
                raise KeyError
            # Calculate percentages
            for i, v in enumerate(parsed_json["duplication"]["histogram"]):
                self.fastp_duplication_plotdata[s_name][i + 1] = (float(v) / float(total_reads)) * 100.0
        except KeyError:
            log.debug(f"No duplication rate plot data: {s_name}")

        # Insert size plot data
        try:
            # First count the total read count in the insert size analysis
            total_reads = 0
            max_i = 0
            for i, v in enumerate(parsed_json["insert_size"]["histogram"]):
                total_reads += v
                if float(v) > 0:
                    max_i = i
            if total_reads == 0:
                raise KeyError
            # Calculate percentages
            for i, v in enumerate(parsed_json["insert_size"]["histogram"]):
                if i <= max_i:
                    self.fastp_insert_size_data[s_name][i + 1] = (float(v) / float(total_reads)) * 100.0
        except KeyError:
            log.debug(f"No insert size plot data: {s_name}")

        for k in [
            "read1_before_filtering",
            "read2_before_filtering",
            "read1_after_filtering",
            "read2_after_filtering",
            "merged_and_filtered",
        ]:
            # Read quality data
            try:
                for i, v in enumerate(parsed_json[k]["quality_curves"]["mean"]):
                    self.fastp_qual_plotdata[k][s_name][i + 1] = float(v)
            except KeyError:
                log.debug(f"Read quality {k} not found: {s_name}")

            # GC and N content plots
            try:
                for i, v in enumerate(parsed_json[k]["content_curves"]["GC"]):
                    self.fastp_gc_content_data[k][s_name][i + 1] = float(v) * 100.0
                for i, v in enumerate(parsed_json[k]["content_curves"]["N"]):
                    self.fastp_n_content_data[k][s_name][i + 1] = float(v) * 100.0
            except KeyError:
                log.debug(f"Content curve data {k} not found: {s_name}")

            # Overrepresented sequences
            try:
                self.fastp_overrepresented_sequences[k][s_name] = parsed_json[k]["overrepresented_sequences"]
            except KeyError:
                log.debug(f"Overrepresented sequences data {k} not found: {s_name}")

        # Remove empty dicts
        if len(self.fastp_data[s_name]) == 0:
            del self.fastp_data[s_name]
        if len(self.fastp_duplication_plotdata[s_name]) == 0:
            del self.fastp_duplication_plotdata[s_name]
        if len(self.fastp_insert_size_data[s_name]) == 0:
            del self.fastp_insert_size_data[s_name]
        if len(self.fastp_all_data[s_name]) == 0:
            del self.fastp_all_data[s_name]
        # Don't delete dicts with subkeys, messes up multi-panel plots

        # Add software version if available
        # Note: this was added to fastp JSON output in v0.22, so it won't be available in older versions
        if "fastp_version" in parsed_json["summary"]:
            self.add_software_version(parsed_json["summary"]["fastp_version"], s_name)

    def fastp_general_stats_table(self, data_by_sample):
        """Take the parsed stats from the fastp report and add it to the
        General Statistics table at the top of the report"""

        self.general_stats_addcols(
            data_by_sample,
            headers={
                "pct_duplication": {
                    "title": "% Duplication",
                    "description": "Duplication rate before filtering",
                    "suffix": "%",
                    "scale": "RdYlGn-rev",
                },
                "after_filtering_q30_rate": {
                    "title": "% > Q30",
                    "description": "Percentage of reads > Q30 after filtering",
                    "modify": lambda x: x * 100.0,
                    "scale": "GnBu",
                    "suffix": "%",
                    "hidden": True,
                },
                "after_filtering_q30_bases": {
                    "title": f"{config.base_count_prefix} Q30 bases",
                    "description": f"Bases > Q30 after filtering ({config.base_count_desc})",
                    "scale": "GnBu",
                    "shared_key": "base_count",
                    "hidden": True,
                },
                "filtering_result_passed_filter_reads": {
                    "title": "Reads After Filtering",
                    "description": f"Total reads after filtering ({config.read_count_desc})",
                    "scale": "Blues",
                    "shared_key": "read_count",
                },
                "after_filtering_gc_content": {
                    "title": "GC content",
                    "description": "GC content after filtering",
                    "suffix": "%",
                    "scale": "Blues",
                    "modify": lambda x: x * 100.0,
                },
                "pct_surviving": {
                    "title": "% PF",
                    "description": "Percent reads passing filter",
                    "suffix": "%",
                    "scale": "RdYlGn",
                    "min": 0,
                    "max": 100,
                },
                "pct_adapter": {
                    "title": "% Adapter",
                    "description": "Percentage adapter-trimmed reads",
                    "suffix": "%",
                    "scale": "RdYlGn-rev",
                    "min": 0,
                    "max": 100,
                },
                "before_filtering_read1_mean_length": {
                    "title": "Mean R1 Length",
                    "description": "Mean read length of R1",
                    "suffix": "bp",
                    "hidden": True,
                },
                "before_filtering_read2_mean_length": {
                    "title": "Mean R2 Length",
                    "description": "Mean read length of R2",
                    "suffix": "bp",
                    "hidden": True,
                },
            },
        )

    def fastp_filtered_reads_chart(self):
        """Function to generate the fastp filtered reads bar plot"""
        # Specify the order of the different possible categories
        keys = {
            "filtering_result_passed_filter_reads": {"name": "Passed Filter"},
            "filtering_result_low_quality_reads": {"name": "Low Quality"},
            "filtering_result_too_many_N_reads": {"name": "Too Many N"},
            "filtering_result_low_complexity_reads": {"name": "Low Complexity"},
            "filtering_result_too_short_reads": {"name": "Too Short"},
            "filtering_result_too_long_reads": {"name": "Too Long"},
        }

        # Config for the plot
        pconfig = {
            "id": "fastp_filtered_reads_plot",
            "title": "Fastp: Filtered Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_zero_cats": False,
            "tt_decimals": 0,
        }
        return bargraph.plot(self.fastp_data, keys, pconfig)

    def fastp_read_qual_plot(self):
        """Make the read quality plot for Fastp"""
        data_labels, pdata = self.filter_pconfig_pdata_subplots(self.fastp_qual_plotdata, "Sequence Quality")
        pconfig = {
            "id": "fastp-seq-quality-plot",
            "title": "Fastp: Sequence Quality",
            "xlab": "Read Position",
            "ylab": "R1 Before filtering: Sequence Quality",
            "ymin": 0,
            "data_labels": data_labels,
        }
        return linegraph.plot(pdata, pconfig)

    def fastp_read_gc_plot(self):
        """Make the read GC plot for Fastp"""
        data_labels, pdata = self.filter_pconfig_pdata_subplots(self.fastp_gc_content_data, "Base Content Percent")
        pconfig = {
            "id": "fastp-seq-content-gc-plot",
            "title": "Fastp: Read GC Content",
            "xlab": "Read Position",
            "ylab": "R1 Before filtering: Base Content Percent",
            "ymax": 100,
            "ymin": 0,
            "tt_label": "{point.x}: {point.y:.2f}%",
            "data_labels": data_labels,
        }
        return linegraph.plot(pdata, pconfig)

    def fastp_read_n_plot(self):
        """Make the read N content plot for Fastp"""
        data_labels, pdata = self.filter_pconfig_pdata_subplots(self.fastp_n_content_data, "Base Content Percent")
        pconfig = {
            "id": "fastp-seq-content-n-plot",
            "title": "Fastp: Read N Content",
            "xlab": "Read Position",
            "ylab": "R1 Before filtering: Base Content Percent",
            "y_clipmax": 100,
            "y_minrange": 5,
            "ymin": 0,
            "tt_label": "{point.x}: {point.y:.2f}%",
            "data_labels": data_labels,
        }
        return linegraph.plot(pdata, pconfig)

    def fastp_overrepresented_sequences_plot(self):
        cnt_by_sample = defaultdict(Counter)
        cnt_by_seq = Counter()
        pct_by_seq = Counter()
        samples_by_seq = defaultdict(set)

        for read_name, by_sample in self.fastp_overrepresented_sequences.items():
            for s_name, by_seq in by_sample.items():
                for seq, count in by_seq.items():
                    cnt_by_seq[seq] += count
                    pct_by_seq[seq] += count / self.fastp_data[s_name]["before_filtering_total_reads"]
                    cnt_by_sample[read_name][s_name] += count
                    samples_by_seq[seq].add(s_name)

        data_labels, cnt_by_sample_pdata = self.filter_pconfig_pdata_subplots(
            cnt_by_sample, "Overrepresented Sequences"
        )
        self.add_section(
            name="Overrepresented sequences by sample",
            anchor="fastp_overrepresented_sequences",
            description="The total amount of overrepresented sequences found in each library.",
            plot=bargraph.plot(
                [
                    {sn: {"overrepresented_sequences": cnt} for sn, cnt in cnt_by_sample.items()}
                    for cnt_by_sample in cnt_by_sample_pdata
                ],
                {
                    "overrepresented_sequences": {"name": "Overrepresented sequences"},
                },
                {
                    "id": "fastp_overrepresented_sequences_plot",
                    "title": "Fastp: Overrepresented sequences",
                    "cpswitch": False,
                    "ylab": "Number of overrepresented sequences",
                    "data_labels": data_labels,
                    "tt_decimals": 0,
                },
            ),
        )

        # Top overrepresented sequences across all samples
        top_n = getattr(config, "fastp_config", {}).get("top_overrepresented_sequences", 20)
        top_seqs = cnt_by_seq.most_common(top_n)
        table_data: Dict[str, Dict[str, Any]] = {
            seq: {
                "sequence": seq,
                "total_percent": pct_by_seq[seq],
                "total_count": cnt_by_seq[seq],
                "samples": len(samples_by_seq[seq]),
            }
            for seq, _ in top_seqs
        }

        self.add_section(
            name="Top overrepresented sequences",
            anchor="fastp_top_overrepresented_sequences",
            description=f"""
            Top overrepresented sequences across all samples. The table shows {top_n}
            most overrepresented sequences across all samples, ranked by the number of occurrences across all samples.
            """,
            plot=table.plot(
                table_data,
                headers={
                    ColumnKey("samples"): {
                        "title": "Reports",
                        "description": "Number of fastp reports where this sequence is found as overrepresented",
                        "scale": "Greens",
                        "min": 0,
                        "format": "{:,d}",
                    },
                    ColumnKey("total_count"): {
                        "title": "Occurrences",
                        "description": "Total number of occurrences of the sequence (among the samples where the sequence is overrepresented)",
                        "scale": "Blues",
                        "min": 0,
                        "format": "{:,d}",
                    },
                    ColumnKey("total_percent"): {
                        "title": "% of all reads",
                        "description": "Total number of occurrences as the percentage of all reads (among samples where the sequence is overrepresented)",
                        "scale": "Blues",
                        "min": 0,
                        "max": 100,
                        "suffix": "%",
                        "format": "{:,.4f}",
                    },
                },
                pconfig={
                    "namespace": self.name,
                    "id": "fastp_top_overrepresented_sequences_table",
                    "title": "fastp: Top overrepresented sequences",
                    "col1_header": "Overrepresented sequence",
                    "sort_rows": False,
                    "rows_are_samples": False,
                },
            ),
        )

    @staticmethod
    def filter_pconfig_pdata_subplots(data, label):
        data_labels = []
        pdata = []
        for k, dl in {
            "read1_before_filtering": {
                "name": "Read 1: Before filtering",
                "ylab": f"R1 Before filtering: {label}",
            },
            "read1_after_filtering": {
                "name": "Read 1: After filtering",
                "ylab": f"R1 After filtering: {label}",
            },
            "read2_before_filtering": {
                "name": "Read 2: Before filtering",
                "ylab": f"R2 Before filtering: {label}",
            },
            "read2_after_filtering": {
                "name": "Read 2: After filtering",
                "ylab": f"R2 After filtering: {label}",
            },
            "merged_and_filtered": {
                "name": "Merged and filtered",
                "ylab": f"{label}",
            },
        }.items():
            if any(data[k][x] for x in data[k]):
                data_labels.append(dl)
                pdata.append(data[k])

        return data_labels, pdata
