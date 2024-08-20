############################################################
######  LOOKING FOR AN EXAMPLE OF HOW MULTIQC WORKS?  ######
############################################################
#### Stop! This is one of the most complicated modules. ####
#### Have a look at Kallisto for a simpler example.     ####
############################################################

import io
import json
import logging
from typing import Dict, Set, List, Union, Literal, Optional

import math
import os
import re
import zipfile
from collections import Counter

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap, linegraph, table
from multiqc import report
from multiqc.plots.plotly.line import Series, LinePlotConfig

log = logging.getLogger(__name__)

VERSION_REGEX = r"FastQC\t([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    """
    FastQC generates an HTML report which is what most people use when
    they run the program. However, it also helpfully generates a file
    called `fastqc_data.txt` which is relatively easy to parse.

    A typical run will produce the following files:

    ```
    mysample_fastqc.html
    mysample_fastqc/
      Icons/
      Images/
      fastqc.fo
      fastqc_data.txt
      fastqc_report.html
      summary.txt
    ```

    Sometimes the directory is zipped, with just `mysample_fastqc.zip`.

    The FastQC MultiQC module looks for files called `fastqc_data.txt`
    or ending in `_fastqc.zip`. If the zip files are found, they are
    read in memory and `fastqc_data.txt` parsed.

    :::note
    The directory and zip file are often both present. To speed
    up MultiQC execution, zip files will be skipped if the file name suggests
    that they will share a sample name with data that has already been parsed.
    :::

    You can customise the patterns used for finding these files in your
    MultiQC config (see [Module search patterns](#module-search-patterns)).
    The below code shows the default file patterns:

    ```yaml
    sp:
      fastqc/data:
        fn: "fastqc_data.txt"
      fastqc/zip:
        fn: "*_fastqc.zip"
    ```

    :::note
    Sample names are discovered by parsing the line beginning
    `Filename` in `fastqc_data.txt`, _not_ based on the FastQC report names.
    :::

    #### Theoretical GC Content

    It is possible to plot a dashed line showing the theoretical GC content for a
    reference genome. MultiQC comes with genome and transcriptome guides for Human
    and Mouse. You can use these in your reports by adding the following MultiQC
    config keys (see [Configuring MultiQC](http://multiqc.info/docs/#configuring-multiqc)):

    ```yaml
    fastqc_config:
      fastqc_theoretical_gc: "hg38_genome"
    ```

    Only one theoretical distribution can be plotted.
    The following guides are available: _(txome = transcriptome)_

    - `hg38_genome`
    - `hg38_txome`
    - `mm10_genome`
    - `mm10_txome`

    Alternatively, a custom theoretical guide can be used in reports. To do this,
    create a file with `fastqc_theoretical_gc` in the filename and place it with your
    analysis files. It should be tab delimited with the following format (column 1 = %GC,
    column 2 = % of genome):

    ```bash
    # FastQC theoretical GC content curve: YOUR REFERENCE NAME
    0	0.005311768
    1	0.004108502
    2	0.004060371
    3	0.005066476
    [...]
    ```

    You can generate these files using an R package called
    [fastqcTheoreticalGC](https://github.com/mikelove/fastqcTheoreticalGC)
    written by [Mike Love](https://github.com/mikelove).
    Please see the [package readme](https://github.com/mikelove/fastqcTheoreticalGC)
    for more details.

    Result files from this package are searched for with the following search pattern
    (can be customised as described above):

    ```yaml
    sp:
      fastqc/theoretical_gc:
        fn: "*fastqc_theoretical_gc*"
    ```

    If you want to always use a specific custom file for MultiQC reports without having to
    add it to the analysis directory, add the full file path to the same MultiQC config
    variable described above:

    ```yaml
    fastqc_config:
      fastqc_theoretical_gc: "/path/to/your/custom_fastqc_theoretical_gc.txt"
    ```

    #### Overrepresented sequences

    The overrepresented sequences table shows the most common sequences found,
    measured by the number of samples they occur as overrepresented. By default, the
    table shows top 20 sequences. This can be customised in the config:

    ```yaml
    fastqc_config:
      top_overrepresented_sequences: 50
    ```

    You can also choose to rank the top sequences by the total number of reads
    rather than by number of samples:

    ```yaml
    fastqc_config:
      top_overrepresented_sequences_by: "total"
    ```

    #### Changing the order of sections

    Remember that it is possible to customise the order in which the different module sections appear
    in the report if you wish.
    See [the docs](https://multiqc.info/docs/#order-of-module-and-module-subsection-output) for more information.

    For example, to show the _Status Checks_ section at the top, use the following config:

    ```yaml
    report_section_order:
      fastqc_status_checks:
        order: -1000
    ```

    #### Showing FastQC status checks

    FastQC uses thresholds to mark samples as "pass", "warn" or "fail" for various checks.
    If you prefer the MultiQC module to ignore those thresholds, and use standard MultiQC
    colors for samples instead, use the following config:

    ```yaml
    fastqc_config:
      status_checks: false
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="FastQC",
            anchor="fastqc",
            href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
            info="Quality control tool for high throughput sequencing data",
            # No publication / DOI // doi=
        )

        self.fastqc_data: Dict = dict()
        self.order_of_duplication_levels: List[float] = []

        # Find and parse unzipped FastQC reports
        for f in self.find_log_files("fastqc/data"):
            s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=os.path.dirname(f["root"]))
            self.parse_fastqc_report(f["f"], s_name, f)

        # Find and parse zipped FastQC reports
        for f in self.find_log_files("fastqc/zip", filecontents=False):
            s_name = f["fn"]
            if s_name.endswith("_fastqc.zip"):
                s_name = s_name[:-11]
            # Skip if we already have this report - parsing zip files is slow
            if s_name in self.fastqc_data.keys():
                log.debug(f"Skipping '{f['fn']}' as already parsed '{s_name}'")
                continue
            try:
                fqc_zip = zipfile.ZipFile(os.path.join(f["root"], f["fn"]))
            except Exception as e:
                log.warning(f"Couldn't read '{f['fn']}' - Bad zip file")
                log.debug(f"Bad zip file error: {e}")
                continue
            # FastQC zip files should have just one directory inside, containing report
            d_name = fqc_zip.namelist()[0]
            try:
                path = os.path.join(d_name, "fastqc_data.txt")
                with fqc_zip.open(path) as fh:
                    r_data_binary = fh.read()
                    try:
                        r_data = r_data_binary.decode("utf8")
                    except UnicodeDecodeError as e:
                        log.debug(f"Could not parse {path} as Unicode: {e}, attempting the latin-1 encoding")
                        try:
                            r_data = r_data_binary.decode("latin-1")
                        except Exception as e:
                            log.warning(f"Error reading FastQC data file {path}: {e}. Skipping sample {s_name}.")
                            continue
                    self.parse_fastqc_report(r_data, s_name, f)
            except KeyError:
                log.warning(f"Error - can't find fastqc_raw_data.txt in {f}")

        # Filter to strip out ignored sample names
        self.fastqc_data = self.ignore_samples(self.fastqc_data)
        if len(self.fastqc_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.fastqc_data)} reports")

        # Write the summary stats to a file
        data = dict()
        for s_name in self.fastqc_data:
            data[s_name] = self.fastqc_data[s_name]["basic_statistics"]
            data[s_name].update(self.fastqc_data[s_name]["statuses"])
        self.write_data_file(data, "multiqc_fastqc")

        # Add to self.css and self.js to be included in template
        self.css = {
            "assets/css/multiqc_fastqc.css": os.path.join(
                os.path.dirname(__file__), "assets", "css", "multiqc_fastqc.css"
            )
        }
        self.js = {
            "assets/js/multiqc_fastqc.js": os.path.join(os.path.dirname(__file__), "assets", "js", "multiqc_fastqc.js")
        }

        # Colours to be used for plotting lines
        self.status_colours = {"pass": "#5cb85c", "warn": "#f0ad4e", "fail": "#d9534f", "default": "#999"}

        # Add to the general statistics table
        self.fastqc_general_stats()

        status_checks = getattr(config, "fastqc_config", {}).get("status_checks", True)

        # Add the statuses to the intro for multiqc_fastqc.js JavaScript to pick up
        statuses: Dict = dict()
        if status_checks:
            for s_name in self.fastqc_data:
                for section, status in self.fastqc_data[s_name]["statuses"].items():
                    try:
                        statuses[section][s_name] = status
                    except KeyError:
                        statuses[section] = {s_name: status}

        self.intro += '<script type="application/json" class="fastqc_passfails">{}</script>'.format(
            json.dumps([self.anchor.replace("-", "_"), statuses])
        )
        if status_checks:
            self.intro += '<script type="text/javascript">load_fastqc_passfails();</script>'

        # Now add each section in order
        self.read_count_plot()
        self.sequence_quality_plot(status_checks)
        self.per_seq_quality_plot(status_checks)
        self.sequence_content_plot()
        self.gc_content_plot(status_checks)
        self.n_content_plot(status_checks)
        self.seq_length_dist_plot(status_checks)
        self.seq_dup_levels_plot(status_checks)
        self.overrepresented_sequences()
        self.adapter_content_plot(status_checks)
        if status_checks:
            self.status_heatmap()
        del self.fastqc_data

    def parse_fastqc_report(self, file_contents, s_name=None, f=None):
        """Takes contents from a fastq_data.txt file and parses out required
        statistics and data. Returns a dict with keys 'stats' and 'data'.
        Data is for plotting graphs, stats are for top table."""

        # Make the sample name from the input filename if we find it
        fn_search = re.search(r"Filename\s+(.+)", file_contents)
        if fn_search:
            s_name = self.clean_s_name(fn_search.group(1), f)

        if s_name in self.fastqc_data.keys():
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.add_data_source(f, s_name)
        self.fastqc_data[s_name] = {"statuses": dict()}

        # Parse the report
        section = None
        s_headers = None
        for line in file_contents.splitlines():
            if line.startswith("##FastQC"):
                version_match = re.search(VERSION_REGEX, line)
                if version_match:
                    self.add_software_version(version_match.group(1), s_name)
            if line == ">>END_MODULE":
                section = None
                s_headers = None
            elif line.startswith(">>"):
                (section, status) = line[2:].split("\t", 1)
                section = section.lower().replace(" ", "_")
                self.fastqc_data[s_name]["statuses"][section] = status
            elif section == "per_tile_sequence_quality":
                # Lots and lots of data. None of it used in MultiQC.
                continue
            elif section is not None:
                if line.startswith("#"):
                    s_headers = line[1:].split("\t")
                    # Special case: Total Deduplicated Percentage header line
                    if s_headers[0] == "Total Deduplicated Percentage":
                        self.fastqc_data[s_name]["basic_statistics"].append(
                            {"measure": "total_deduplicated_percentage", "value": float(s_headers[1])}
                        )
                    else:
                        # Special case: Rename dedup header in old versions of FastQC (v10)
                        if s_headers[1] == "Relative count":
                            s_headers[1] = "Percentage of total"
                        s_headers = [s.lower().replace(" ", "_") for s in s_headers]
                        self.fastqc_data[s_name][section] = list()

                elif s_headers is not None:
                    s = line.split("\t")
                    row = dict()
                    for i, v in enumerate(s):
                        v.replace("NaN", "0")
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                        row[s_headers[i]] = v
                    self.fastqc_data[s_name][section].append(row)
                    # Special case - need to remember order of duplication keys
                    if section == "sequence_duplication_levels":
                        try:
                            level = float(s[0])
                        except ValueError:
                            level = s[0]
                        if level not in self.order_of_duplication_levels:
                            self.order_of_duplication_levels.append(level)

        # Tidy up the Basic Stats
        self.fastqc_data[s_name]["basic_statistics"] = {
            d["measure"]: d["value"] for d in self.fastqc_data[s_name]["basic_statistics"]
        }

        # We sort by the mean of the range, which is effectively sorting ranges in asc order assuming no overlap
        sequence_length_distributions = self.fastqc_data[s_name].get("sequence_length_distribution", [])
        sequence_length_distributions.sort(key=lambda d: _range_bp_to_num(d["length"], method="mean"))

        # Calculate the average sequence length (Basic Statistics gives a range)
        length_reads = 0
        length_bp = 0
        total_count = sum(d["count"] for d in sequence_length_distributions)
        median: Optional[int] = None

        for d in sequence_length_distributions:
            length_reads += d["count"]
            length_bp += d["count"] * _range_bp_to_num(d["length"], method="mean")

            if median is None and length_reads >= total_count / 2:
                # if the distribution-entry is a range, we use the average of the range.
                # this isn't technically correct, because we can't know what the distribution
                # is within that range. Probably good enough though.
                median = int(_range_bp_to_num(d["length"], method="median"))
        if total_count > 0:
            self.fastqc_data[s_name]["basic_statistics"]["avg_sequence_length"] = length_bp / total_count
        if median is not None:
            self.fastqc_data[s_name]["basic_statistics"]["median_sequence_length"] = median

    def fastqc_general_stats(self):
        """Add some single-number stats to the basic statistics
        table at the top of the report"""

        # Prep the data
        data_by_sample: Dict[str, Dict[str, Union[int, float, str]]] = dict()
        for s_name, sd in self.fastqc_data.items():
            data_by_sample[s_name] = dict()
            bs = sd["basic_statistics"]
            # Samples with 0 reads and reports with some skipped sections might be missing things here
            data_by_sample[s_name]["percent_gc"] = bs.get("%GC", 0)
            data_by_sample[s_name]["avg_sequence_length"] = bs.get("avg_sequence_length", 0)
            data_by_sample[s_name]["median_sequence_length"] = bs.get("median_sequence_length", 0)
            data_by_sample[s_name]["total_sequences"] = bs.get("Total Sequences", 0)

            # Log warning about zero-read samples as a courtesy
            if data_by_sample[s_name]["total_sequences"] == 0:
                log.warning(f"Sample had zero reads: '{s_name}'")

            if "total_deduplicated_percentage" in bs:
                # Older versions of FastQC don't have this
                data_by_sample[s_name]["percent_duplicates"] = 100 - bs["total_deduplicated_percentage"]

            # Add count of fail statuses
            num_statuses = 0
            num_fails = 0
            for s in sd["statuses"].values():
                num_statuses += 1
                if s == "fail":
                    num_fails += 1
            if num_statuses > 0:
                # Make sure we have reads, otherwise there are no sample in data
                data_by_sample[s_name]["percent_fails"] = (float(num_fails) / float(num_statuses)) * 100.0

        # Are sequence lengths interesting?
        median_seq_lengths = [x["median_sequence_length"] for x in data_by_sample.values()]
        try:
            hide_seq_length = max(median_seq_lengths) - min(median_seq_lengths) <= 10
        except ValueError:
            # Zero reads
            hide_seq_length = True

        headers = {
            "percent_duplicates": {
                "title": "% Dups",
                "description": "% Duplicate Reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn-rev",
            },
            "percent_gc": {
                "title": "% GC",
                "description": "Average % GC Content",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "PuRd",
                "format": "{:,.0f}",
            },
            "avg_sequence_length": {
                "title": "Average Read Length",
                "description": "Average Read Length (bp)",
                "min": 0,
                "suffix": " bp",
                "scale": "RdYlGn",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "median_sequence_length": {
                "title": "Median Read Length",
                "description": "Median Read Length (bp)",
                "min": 0,
                "suffix": " bp",
                "scale": "RdYlGn",
                "format": "{:,.0f}",
                "hidden": hide_seq_length,
            },
            "percent_fails": {
                "title": "% Failed",
                "description": "Percentage of modules failed in FastQC report (includes those not plotted here)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "Reds",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "total_sequences": {
                "title": f"{config.read_count_prefix} Seqs",
                "description": f"Total Sequences ({config.read_count_desc})",
                "min": 0,
                "scale": "Blues",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

    def read_count_plot(self):
        """Stacked bar plot showing counts of reads"""

        # Main plot config
        pconfig = {
            "id": "fastqc_sequence_counts_plot",
            "title": "FastQC: Sequence Counts",
            "ylab": "Number of reads",
            "cpswitch_counts_label": "Number of reads",
            "hide_empty": False,
        }

        # Calculate the number of unique and duplicate reads if we can
        data_by_sample: Dict[str, Dict[str, int]] = dict()
        has_dups = False
        has_total = False
        for s_name, sd in self.fastqc_data.items():
            pd = sd["basic_statistics"]
            data_by_sample[s_name] = dict()
            try:
                data_by_sample[s_name]["Duplicate Reads"] = int(
                    ((100.0 - float(pd["total_deduplicated_percentage"])) / 100.0) * pd["Total Sequences"]
                )
                data_by_sample[s_name]["Unique Reads"] = (
                    pd["Total Sequences"] - data_by_sample[s_name]["Duplicate Reads"]
                )
                has_dups = True
            # Older versions of FastQC don't have duplicate reads
            # Very sparse data can report -nan: #Total Deduplicated Percentage  -nan
            except (KeyError, ValueError):
                data_by_sample[s_name] = {"Total Sequences": pd["Total Sequences"]}
                has_total = True

        # Configure the cats and config according to what we found
        pcats = list()
        duptext = ""
        if has_total:
            pcats.append("Total Sequences")
        if has_dups:
            pcats.extend(["Unique Reads", "Duplicate Reads"])
            duptext = " Duplicate read counts are an estimate only."
        if has_total and not has_dups:
            pconfig["use_legend"] = False
            pconfig["cpswitch"] = False

        # Add the report section
        self.add_section(
            name="Sequence Counts",
            anchor="fastqc_sequence_counts",
            description="Sequence counts for each sample." + duptext,
            helptext="""
            This plot show the total number of reads, broken down into unique and duplicate
            if possible (only more recent versions of FastQC give duplicate info).

            You can read more about duplicate calculation in the
            [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html).
            A small part has been copied here for convenience:

            _Only sequences which first appear in the first 100,000 sequences
            in each file are analysed. This should be enough to get a good impression
            for the duplication levels in the whole file. Each sequence is tracked to
            the end of the file to give a representative count of the overall duplication level._

            _The duplication detection requires an exact sequence match over the whole length of
            the sequence. Any reads over 75bp in length are truncated to 50bp for this analysis._
            """,
            plot=bargraph.plot(data_by_sample, pcats, pconfig),
        )

    def sequence_quality_plot(self, status_checks=True):
        """Create the HTML for the phred quality score plot"""

        data_by_sample: Dict[str, Dict[int, float]] = dict()
        for s_name, sd in self.fastqc_data.items():
            if "per_base_sequence_quality" not in sd:
                continue
            data_by_sample[s_name] = {
                int(_range_bp_to_num(d["base"], method="start")): d["mean"] for d in sd["per_base_sequence_quality"]
            }
        if len(data_by_sample) == 0:
            log.debug("sequence_quality not found in FastQC reports")
            return None

        pconfig = {
            "id": "fastqc_per_base_sequence_quality_plot",
            "title": "FastQC: Mean Quality Scores",
            "ylab": "Phred Score",
            "xlab": "Position (bp)",
            "ymin": 0,
            "xmin": 0,
            "x_decimals": False,
            "tt_label": "<b>Base {point.x}</b>: {point.y:.2f}",
            "showlegend": False if status_checks else True,
        }
        if status_checks:
            pconfig.update(
                {
                    "colors": self.get_status_cols("per_base_sequence_quality"),
                    "y_bands": [
                        {"from": 28, "to": 100, "color": "#c3e6c3"},
                        {"from": 20, "to": 28, "color": "#e6dcc3"},
                        {"from": 0, "to": 20, "color": "#e6c3c3"},
                    ],
                }
            )

        self.add_section(
            name="Sequence Quality Histograms",
            anchor="fastqc_per_base_sequence_quality",
            description="The mean quality value across each base position in the read.",
            helptext="""
            To enable multiple samples to be plotted on the same graph, only the mean quality
            scores are plotted (unlike the box plots seen in FastQC reports).

            Taken from the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html):

            _The y-axis on the graph shows the quality scores. The higher the score, the better
            the base call. The background of the graph divides the y axis into very good quality
            calls (green), calls of reasonable quality (orange), and calls of poor quality (red).
            The quality of calls on most platforms will degrade as the run progresses, so it is
            common to see base calls falling into the orange area towards the end of a read._
            """,
            plot=linegraph.plot(data_by_sample, pconfig),
        )

    def per_seq_quality_plot(self, status_checks=True):
        """Create the HTML for the per sequence quality score plot"""

        data_by_sample: Dict[str, Dict[int, float]] = dict()
        for s_name, sd in self.fastqc_data.items():
            if "per_sequence_quality_scores" not in sd:
                continue
            data_by_sample[s_name] = {d["quality"]: d["count"] for d in sd["per_sequence_quality_scores"]}
        if len(data_by_sample) == 0:
            log.debug("per_seq_quality not found in FastQC reports")
            return None

        pconfig = {
            "id": "fastqc_per_sequence_quality_scores_plot",
            "title": "FastQC: Per Sequence Quality Scores",
            "ylab": "Count",
            "xlab": "Mean Sequence Quality (Phred Score)",
            "ymin": 0,
            "xmin": 0,
            "x_decimals": False,
            "tt_label": "<b>Phred {point.x}</b>: {point.y} reads",
            "showlegend": False if status_checks else True,
        }
        if status_checks:
            pconfig.update(
                {
                    "colors": self.get_status_cols("per_sequence_quality_scores"),
                    "x_bands": [
                        {"from": 28, "to": 100, "color": "#c3e6c3"},
                        {"from": 20, "to": 28, "color": "#e6dcc3"},
                        {"from": 0, "to": 20, "color": "#e6c3c3"},
                    ],
                }
            )
        self.add_section(
            name="Per Sequence Quality Scores",
            anchor="fastqc_per_sequence_quality_scores",
            description="The number of reads with average quality scores. Shows if a subset of reads has poor quality.",
            helptext="""
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html):

            _The per sequence quality score report allows you to see if a subset of your
            sequences have universally low quality values. It is often the case that a
            subset of sequences will have universally poor quality, however these should
            represent only a small percentage of the total sequences._
            """,
            plot=linegraph.plot(data_by_sample, pconfig),
        )

    def sequence_content_plot(self):
        """Create the epic HTML for the FastQC sequence content heatmap"""

        data_by_sample: Dict[str, Dict[int, Dict[str, int]]] = dict()
        for s_name in sorted(self.fastqc_data.keys()):
            if "per_base_sequence_content" not in self.fastqc_data[s_name]:
                continue

            data_by_sample[s_name] = {
                int(_range_bp_to_num(d["base"], method="start")): d
                for d in self.fastqc_data[s_name]["per_base_sequence_content"]
            }

            # Replace NaN with 0
            for b in data_by_sample[s_name]:
                for base in ["a", "c", "t", "g"]:
                    if math.isnan(float(data_by_sample[s_name][b][base])):
                        data_by_sample[s_name][b][base] = 0

        if len(data_by_sample) == 0:
            log.debug("sequence_content not found in FastQC reports")
            return None

        html = """<div id="fastqc_per_base_sequence_content_plot_div">
            <div class="alert alert-info">
               <span class="glyphicon glyphicon-hand-up"></span>
               Click a sample row to see a line plot for that dataset.
            </div>
            <h5><span class="s_name text-primary"><span class="glyphicon glyphicon-info-sign"></span> Rollover for sample name</span></h5>
            <div class="fastqc_seq_heatmap_key">
                Position: <span id="fastqc_seq_heatmap_key_pos">-</span>
                <div><span id="fastqc_seq_heatmap_key_t"> %T: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_c"> %C: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_a"> %A: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_g"> %G: <span>-</span></span></div>
            </div>
            <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">
                <div id="{id}" class="fastqc_per_base_sequence_content_plot hc-plot has-custom-export">
                    <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas>
                </div>
            </div>
            <div class="clearfix"></div>
        </div>
        <script type="application/json" class="fastqc_seq_content">{d}</script>
        """.format(
            # Generate unique plot ID, needed in mqc_export_selectplots
            id=report.save_htmlid("fastqc_per_base_sequence_content_plot"),
            d=json.dumps([self.anchor.replace("-", "_"), data_by_sample]),
        )

        self.add_section(
            name="Per Base Sequence Content",
            anchor="fastqc_per_base_sequence_content",
            description="The proportion of each base position for which each of the four normal DNA bases has been called.",
            helptext="""
            To enable multiple samples to be shown in a single plot, the base composition data
            is shown as a heatmap. The colours represent the balance between the four bases:
            an even distribution should give an even muddy brown colour. Hover over the plot
            to see the percentage of the four bases under the cursor.

            **To see the data as a line plot, as in the original FastQC graph, click on a sample track.**

            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html):

            _Per Base Sequence Content plots out the proportion of each base position in a
            file for which each of the four normal DNA bases has been called._

            _In a random library you would expect that there would be little to no difference
            between the different bases of a sequence run, so the lines in this plot should
            run parallel with each other. The relative amount of each base should reflect
            the overall amount of these bases in your genome, but in any case they should
            not be hugely imbalanced from each other._

            _It's worth noting that some types of library will always produce biased sequence
            composition, normally at the start of the read. Libraries produced by priming
            using random hexamers (including nearly all RNA-Seq libraries) and those which
            were fragmented using transposases inherit an intrinsic bias in the positions
            at which reads start. This bias does not concern an absolute sequence, but instead
            provides enrichement of a number of different K-mers at the 5' end of the reads.
            Whilst this is a true technical bias, it isn't something which can be corrected
            by trimming and in most cases doesn't seem to adversely affect the downstream
            analysis._
            """,
            content=html,
        )

    def gc_content_plot(self, status_checks=True):
        """Create the HTML for the FastQC GC content plot"""

        data_by_sample: Dict[str, Dict[int, float]] = dict()
        data_norm_by_sample: Dict[str, Dict[int, float]] = dict()
        for s_name, sd in self.fastqc_data.items():
            if "per_sequence_gc_content" not in sd:
                continue

            data_by_sample[s_name] = {d["gc_content"]: d["count"] for d in sd["per_sequence_gc_content"]}
            data_norm_by_sample[s_name] = dict()
            total = sum([c for c in data_by_sample[s_name].values()])
            for gc, count in data_by_sample[s_name].items():
                if total == 0:
                    data_norm_by_sample[s_name][gc] = 0
                else:
                    data_norm_by_sample[s_name][gc] = (count / total) * 100
        if len(data_by_sample) == 0:
            log.debug("per_sequence_gc_content not found in FastQC reports")
            return None

        pconfig = {
            "id": "fastqc_per_sequence_gc_content_plot",
            "title": "FastQC: Per Sequence GC Content",
            "xlab": "% GC",
            "ylab": "Percentage",
            "ymin": 0,
            "xmax": 100,
            "xmin": 0,
            "tt_label": "<b>{point.x}% GC</b>: {point.y}",
            "data_labels": [
                {"name": "Percentages", "ylab": "Percentage", "tt_suffix": "%"},
                {"name": "Counts", "ylab": "Count", "tt_suffix": ""},
            ],
            "showlegend": False if status_checks else True,
        }
        if status_checks:
            pconfig.update(
                {
                    "colors": self.get_status_cols("per_sequence_gc_content"),
                }
            )

        # Try to find and plot a theoretical GC line
        theoretical_gc = None
        theoretical_gc_raw = None
        theoretical_gc_name = None
        for f in self.find_log_files("fastqc/theoretical_gc"):
            if theoretical_gc_raw is not None:
                log.warning(f"Multiple FastQC Theoretical GC Content files found, now using {f['fn']}")
            theoretical_gc_raw = f["f"]
            theoretical_gc_name = f["fn"]
        if theoretical_gc_raw is None:
            tgc = getattr(config, "fastqc_config", {}).get("fastqc_theoretical_gc")
            if tgc is not None and isinstance(tgc, str):
                theoretical_gc_name = os.path.basename(tgc)
                tgc_fn = f"fastqc_theoretical_gc_{tgc}.txt"
                tgc_path = os.path.join(os.path.dirname(__file__), "fastqc_theoretical_gc", tgc_fn)
                if not os.path.isfile(tgc_path):
                    tgc_path = tgc
                try:
                    with io.open(tgc_path, "r", encoding="utf-8") as f:
                        theoretical_gc_raw = f.read()
                except IOError:
                    log.warning(f"Couldn't open FastQC Theoretical GC Content file {tgc_path}")
                    theoretical_gc_raw = None
        if theoretical_gc_raw is not None:
            theoretical_gc = list()
            for line in theoretical_gc_raw.splitlines():
                if "# FastQC theoretical GC content curve:" in line:
                    theoretical_gc_name = line[39:]
                elif not line.startswith("#"):
                    s = line.split()
                    try:
                        theoretical_gc.append((float(s[0]), float(s[1])))
                    except (TypeError, IndexError):
                        pass

        desc = """The average GC content of reads. Normal random library typically have a
        roughly normal distribution of GC content."""
        if theoretical_gc is not None:
            # Calculate the count version of the theoretical data based on the largest data store
            max_total = max([sum(d.values()) for d in data_by_sample.values()])
            extra_series_config = {
                "name": "Theoretical GC Content",
                "dash": "dash",
                "width": 2,
                "color": "black",
                "showlegend": False if status_checks else True,
            }
            s1 = Series(pairs=theoretical_gc, **extra_series_config)
            s2 = Series(pairs=[(d[0], (d[1] / 100.0) * max_total) for d in theoretical_gc], **extra_series_config)
            pconfig["extra_series"] = [[s1], [s2]]
            desc = f" **The dashed black line shows theoretical GC content:** `{theoretical_gc_name}`"

        self.add_section(
            name="Per Sequence GC Content",
            anchor="fastqc_per_sequence_gc_content",
            description=desc,
            helptext="""
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html):

            _This module measures the GC content across the whole length of each sequence
            in a file and compares it to a modelled normal distribution of GC content._

            _In a normal random library you would expect to see a roughly normal distribution
            of GC content where the central peak corresponds to the overall GC content of
            the underlying genome. Since we don't know the GC content of the genome the
            modal GC content is calculated from the observed data and used to build a
            reference distribution._

            _An unusually shaped distribution could indicate a contaminated library or
            some other kinds of biased subset. A normal distribution which is shifted
            indicates some systematic bias which is independent of base position. If there
            is a systematic bias which creates a shifted normal distribution then this won't
            be flagged as an error by the module since it doesn't know what your genome's
            GC content should be._
            """,
            plot=linegraph.plot([data_norm_by_sample, data_by_sample], pconfig),
        )

    def n_content_plot(self, status_checks=True):
        """Create the HTML for the per base N content plot"""

        data_by_sample: Dict[str, Dict[int, int]] = dict()
        for s_name, sd in self.fastqc_data.items():
            if "per_base_n_content" not in sd:
                continue
            data_by_sample[s_name] = {
                int(_range_bp_to_num(d["base"], method="start")): d["n-count"] for d in sd["per_base_n_content"]
            }
        if len(data_by_sample) == 0:
            log.debug("per_base_n_content not found in FastQC reports")
            return None

        pconfig = {
            "id": "fastqc_per_base_n_content_plot",
            "title": "FastQC: Per Base N Content",
            "ylab": "Percentage N-Count",
            "xlab": "Position in Read (bp)",
            "y_clipmax": 100,
            "y_minrange": 5,
            "ymin": 0,
            "xmin": 0,
            "tt_label": "<b>Base {point.x}</b>: {point.y:.2f}%",
            "showlegend": False if status_checks else True,
        }
        if status_checks:
            pconfig.update(
                {
                    "colors": self.get_status_cols("per_base_n_content"),
                    "y_bands": [
                        {"from": 20, "to": 100, "color": "#e6c3c3"},
                        {"from": 5, "to": 20, "color": "#e6dcc3"},
                        {"from": 0, "to": 5, "color": "#c3e6c3"},
                    ],
                }
            )

        self.add_section(
            name="Per Base N Content",
            anchor="fastqc_per_base_n_content",
            description="The percentage of base calls at each position for which an `N` was called.",
            helptext="""
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html):

            _If a sequencer is unable to make a base call with sufficient confidence then it will
            normally substitute an `N` rather than a conventional base call. This graph shows the
            percentage of base calls at each position for which an `N` was called._

            _It's not unusual to see a very low proportion of Ns appearing in a sequence, especially
            nearer the end of a sequence. However, if this proportion rises above a few percent
            it suggests that the analysis pipeline was unable to interpret the data well enough to
            make valid base calls._
            """,
            plot=linegraph.plot(data_by_sample, pconfig),
        )

    def seq_length_dist_plot(self, status_checks):
        """Create the HTML for the Sequence Length Distribution plot"""

        cnt_by_range_by_sample: Dict[str, Dict[int, int]] = dict()
        all_ranges_across_samples: Set[int] = set()
        only_single_length: bool = True
        for s_name, sd in self.fastqc_data.items():
            if "sequence_length_distribution" not in sd:
                continue
            cnt_by_range_by_sample[s_name] = {
                int(_range_bp_to_num(d["length"], method="start")): d["count"]
                for d in sd["sequence_length_distribution"]
            }
            sample_ranges_set = set(cnt_by_range_by_sample[s_name].keys())
            if len(sample_ranges_set) > 1:
                only_single_length = False
            all_ranges_across_samples.update(sample_ranges_set)
        if len(cnt_by_range_by_sample) == 0:
            log.debug("sequence_length_distribution not found in FastQC reports")
            return None

        if only_single_length:
            lengths_line = ", ".join([f"{length:,.0f}bp" for length in list(all_ranges_across_samples)])
            desc = f"All samples have sequences of a single length ({lengths_line})"
            if len(all_ranges_across_samples) > 1:
                desc += ' See the <a href="#general_stats">General Statistics Table</a>.'
            self.add_section(
                name="Sequence Length Distribution",
                anchor="fastqc_sequence_length_distribution",
                description=f'<div class="alert alert-info">{desc}</div>',
            )
        else:
            pconfig = LinePlotConfig(
                id="fastqc_sequence_length_distribution_plot",
                title="FastQC: Sequence Length Distribution",
                ylab="Read Count",
                xlab="Sequence Length (bp)",
                ymin=0,
                tt_label="<b>{point.x} bp</b>: {point.y}",
                showlegend=False if status_checks else True,
            )
            if status_checks:
                pconfig.colors = self.get_status_cols("sequence_length_distribution")
            self.add_section(
                name="Sequence Length Distribution",
                anchor="fastqc_sequence_length_distribution",
                description="The distribution of fragment sizes (read lengths) found. See the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)",
                plot=linegraph.plot(cnt_by_range_by_sample, pconfig),
            )

    def seq_dup_levels_plot(self, status_checks):
        """Create the HTML for the Sequence Duplication Levels plot"""

        data: Dict = dict()
        max_dup_val = 0
        for s_name in self.fastqc_data:
            try:
                thisdata = {}
                for d in self.fastqc_data[s_name]["sequence_duplication_levels"]:
                    thisdata[d["duplication_level"]] = d["percentage_of_total"]
                    max_dup_val = max(max_dup_val, d["percentage_of_total"])
                data[s_name] = {}
                for k in self.order_of_duplication_levels:
                    try:
                        data[s_name][k] = thisdata[k]
                    except KeyError:
                        pass
            except KeyError:
                pass
        if len(data) == 0:
            log.debug("sequence_length_distribution not found in FastQC reports")
            return None
        pconfig = {
            "id": "fastqc_sequence_duplication_levels_plot",
            "title": "FastQC: Sequence Duplication Levels",
            "categories": True,
            "ylab": "% of Library",
            "xlab": "Sequence Duplication Level",
            "ymax": 100 if max_dup_val <= 100.0 else None,
            "ymin": 0,
            "tt_decimals": 2,
            "tt_suffix": "%",
            "showlegend": False if status_checks else True,
        }
        if status_checks:
            pconfig.update(
                {
                    "colors": self.get_status_cols("sequence_duplication_levels"),
                }
            )

        self.add_section(
            name="Sequence Duplication Levels",
            anchor="fastqc_sequence_duplication_levels",
            description="The relative level of duplication found for every sequence.",
            helptext="""
            From the [FastQC Help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html):

            _In a diverse library most sequences will occur only once in the final set.
            A low level of duplication may indicate a very high level of coverage of the
            target sequence, but a high level of duplication is more likely to indicate
            some kind of enrichment bias (e.g. PCR over amplification). This graph shows
            the degree of duplication for every sequence in a library: the relative
            number of sequences with different degrees of duplication._

            _Only sequences which first appear in the first 100,000 sequences
            in each file are analysed. This should be enough to get a good impression
            for the duplication levels in the whole file. Each sequence is tracked to
            the end of the file to give a representative count of the overall duplication level._

            _The duplication detection requires an exact sequence match over the whole length of
            the sequence. Any reads over 75bp in length are truncated to 50bp for this analysis._

            _In a properly diverse library most sequences should fall into the far left of the
            plot in both the red and blue lines. A general level of enrichment, indicating broad
            oversequencing in the library will tend to flatten the lines, lowering the low end
            and generally raising other categories. More specific enrichments of subsets, or
            the presence of low complexity contaminants will tend to produce spikes towards the
            right of the plot._
            """,
            plot=linegraph.plot(data, pconfig),
        )

    def overrepresented_sequences(self):
        """Sum the percentages of overrepresented sequences and display them in a bar plot"""

        data: Dict = dict()
        # Count the number of samples where a sequence is overrepresented
        overrep_by_sample: Counter = Counter()
        overrep_total_cnt: Counter = Counter()
        for s_name in self.fastqc_data:
            data[s_name] = dict()
            try:
                max_pcnt = max([float(d["percentage"]) for d in self.fastqc_data[s_name]["overrepresented_sequences"]])
                total_pcnt = sum(
                    [float(d["percentage"]) for d in self.fastqc_data[s_name]["overrepresented_sequences"]]
                )
                data[s_name]["total_overrepresented"] = total_pcnt
                data[s_name]["top_overrepresented"] = max_pcnt
                data[s_name]["remaining_overrepresented"] = total_pcnt - max_pcnt
                for d in self.fastqc_data[s_name]["overrepresented_sequences"]:
                    overrep_by_sample[d["sequence"]] += 1
                    overrep_total_cnt[d["sequence"]] += int(d["count"])
            except KeyError:
                if self.fastqc_data[s_name]["statuses"].get("overrepresented_sequences") == "pass":
                    data[s_name]["total_overrepresented"] = 0
                    data[s_name]["top_overrepresented"] = 0
                    data[s_name]["remaining_overrepresented"] = 0
                    data[s_name]["overrepresented_sequences"] = []
                else:
                    del data[s_name]
                    log.debug(f"Couldn't find data for {s_name}, invalid Key")

        if all(len(data.get(s_name, {})) == 0 for s_name in self.fastqc_data):
            log.debug("overrepresented_sequences not found in FastQC reports")
            return None

        cats = {
            "top_overrepresented": {"name": "Top overrepresented sequence"},
            "remaining_overrepresented": {"name": "Sum of remaining overrepresented sequences"},
        }

        # Config for the plot
        pconfig = {
            "id": "fastqc_overrepresented_sequences_plot",
            "title": "FastQC: Overrepresented sequences sample summary",
            "ymin": 0,
            "y_clipmax": 100,
            "tt_decimals": 2,
            "tt_suffix": "%",
            "cpswitch": False,
            "ylab": "Percentage of Total Sequences",
        }

        plot = None
        content = None
        # Check if any samples have more than 1% overrepresented sequences, else don't make plot.
        if max([x["total_overrepresented"] for x in data.values()]) < 1:
            content = f'<div class="alert alert-info">{len(data)} samples had less than 1% of reads made up of overrepresented sequences</div>'
        else:
            plot = bargraph.plot(data, cats, pconfig)

        self.add_section(
            name="Overrepresented sequences by sample",
            anchor="fastqc_overrepresented_sequences",
            description="The total amount of overrepresented sequences found in each library.",
            helptext="""
            FastQC calculates and lists overrepresented sequences in FastQ files. It would not be
            possible to show this for all samples in a MultiQC report, so instead this plot shows
            the _number of sequences_ categorized as overrepresented.

            Sometimes, a single sequence  may account for a large number of reads in a dataset.
            To show this, the bars are split into two: the first shows the overrepresented reads
            that come from the single most common sequence. The second shows the total count
            from all remaining overrepresented sequences.

            From the [FastQC Help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html):

            _A normal high-throughput library will contain a diverse set of sequences, with no
            individual sequence making up a tiny fraction of the whole. Finding that a single
            sequence is very overrepresented in the set either means that it is highly biologically
            significant, or indicates that the library is contaminated, or not as diverse as you expected._

            _FastQC lists all the sequences which make up more than 0.1% of the total.
            To conserve memory only sequences which appear in the first 100,000 sequences are tracked
            to the end of the file. It is therefore possible that a sequence which is overrepresented
            but doesn't appear at the start of the file for some reason could be missed by this module._
            """,
            plot=plot,
            content=content,
        )

        # Add a table of the top overrepresented sequences
        # Recalculate counts to percentages for readability:
        total_read_count = sum([int(d["basic_statistics"]["Total Sequences"]) for d in self.fastqc_data.values()])
        overrep_total_pct = {seq: (cnt / total_read_count) * 100 for seq, cnt in overrep_total_cnt.items()}

        # Top overrepresented sequences across all samples
        top_n = getattr(config, "fastqc_config", {}).get("top_overrepresented_sequences", 20)
        by = getattr(config, "fastqc_config", {}).get("top_overrepresented_sequences_by", "samples")
        if by == "samples":
            top_seqs = overrep_by_sample.most_common(top_n)
        else:
            top_seqs = overrep_total_cnt.most_common(top_n)
        headers = {
            "samples": {
                "title": "Samples",
                "description": "Number of samples where this sequence is overrepresented",
                "scale": "Greens",
                "min": 0,
                "format": "{:,d}",
            },
            "total_count": {
                "title": "Occurrences",
                "description": "Total number of occurrences of the sequence (among the samples where the sequence is overrepresented)",
                "scale": "Blues",
                "min": 0,
                "format": "{:,d}",
            },
            "total_percent": {
                "title": "% of all reads",
                "description": "Total number of occurrences as the percentage of all reads (among samples where the sequence is overrepresented)",
                "scale": "Blues",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "format": "{:,.4f}",
            },
        }
        data = {
            seq: {
                "sequence": seq,
                "total_percent": overrep_total_pct[seq],
                "total_count": overrep_total_cnt[seq],
                "samples": overrep_by_sample[seq],
            }
            for seq, _ in top_seqs
        }

        ranked_by = (
            "the number of samples they occur in" if by == "samples" else "the number of occurrences across all samples"
        )
        self.add_section(
            name="Top overrepresented sequences",
            anchor="fastqc_top_overrepresented_sequences",
            description=f"""
            Top overrepresented sequences across all samples. The table shows {top_n} 
            most overrepresented sequences across all samples, ranked by {ranked_by}.
            """,
            plot=table.plot(
                data,
                headers,
                {
                    "namespace": self.name,
                    "id": "fastqc_top_overrepresented_sequences_table",
                    "title": "FastQC: Top overrepresented sequences",
                    "col1_header": "Overrepresented sequence",
                    "sort_rows": False,
                },
            ),
        )

    def adapter_content_plot(self, status_checks=True):
        """Create the HTML for the FastQC adapter plot"""

        pct_by_pos_by_sample: Dict[str, Dict[int, int]] = dict()
        for s_name, data_by_sample in self.fastqc_data.items():
            if "adapter_content" not in data_by_sample:
                continue
            for adapters in data_by_sample["adapter_content"]:
                pos = int(
                    _range_bp_to_num(adapters["position"], method="start")
                )  # split ranges like "10-15", take start
                for adapter_name, percent in adapters.items():
                    k = f"{s_name} - {adapter_name}"
                    if adapter_name != "position":
                        pct_by_pos_by_sample.setdefault(k, {})[pos] = percent
        if len(pct_by_pos_by_sample) == 0:
            log.debug("adapter_content not found in FastQC reports")
            return None

        # Lots of these datasets will be all zeros.
        # Only take datasets with > 0.1% adapter contamination
        pct_by_pos_by_sample = {
            k: d for k, d in pct_by_pos_by_sample.items() if max(pct_by_pos_by_sample[k].values()) >= 0.1
        }

        pconfig = {
            "id": "fastqc_adapter_content_plot",
            "title": "FastQC: Adapter Content",
            "ylab": "% of Sequences",
            "xlab": "Position (bp)",
            "y_clipmax": 100,
            "y_minrange": 5,
            "ymin": 0,
            "tt_label": "<b>Base {point.x}</b>: {point.y:.2f}%",
            "hide_empty": True,
        }
        if status_checks:
            pconfig.update(
                {
                    "y_bands": [
                        {"from": 20, "to": 100, "color": "#e6c3c3"},
                        {"from": 5, "to": 20, "color": "#e6dcc3"},
                        {"from": 0, "to": 5, "color": "#c3e6c3"},
                    ],
                }
            )

        plot = None
        content = None
        if len(pct_by_pos_by_sample) > 0:
            plot = linegraph.plot(pct_by_pos_by_sample, pconfig)
        else:
            content = '<div class="alert alert-info">No samples found with any adapter contamination > 0.1%</div>'

        # Note - colours are messy as we've added adapter names here. Not
        # possible to break down pass / warn / fail for each adapter, which
        # could lead to misleading labelling (fails on adapter types with
        # little or no contamination)

        self.add_section(
            name="Adapter Content",
            anchor="fastqc_adapter_content",
            description="""The cumulative percentage count of the proportion of your
            library which has seen each of the adapter sequences at each position.""",
            helptext="""
            Note that only samples with ≥ 0.1% adapter contamination are shown.

            There may be several lines per sample, as one is shown for each adapter
            detected in the file.

            From the [FastQC Help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html):

            _The plot shows a cumulative percentage count of the proportion
            of your library which has seen each of the adapter sequences at each position.
            Once a sequence has been seen in a read it is counted as being present
            right through to the end of the read so the percentages you see will only
            increase as the read length goes on._
            """,
            plot=plot,
            content=content,
        )

    def status_heatmap(self):
        """Heatmap showing all statuses for every sample"""
        status_numbers = {"pass": 1, "warn": 0.5, "fail": 0.25}
        data = []
        s_names = []
        status_cats = {}
        for s_name in sorted(self.fastqc_data.keys()):
            s_names.append(s_name)
            for status_cat, status in self.fastqc_data[s_name]["statuses"].items():
                if status_cat not in status_cats:
                    status_cats[status_cat] = status_cat.replace("_", " ").title().replace("Gc", "GC")
        for s_name in s_names:
            row = []
            for status_cat in status_cats:
                try:
                    row.append(status_numbers[self.fastqc_data[s_name]["statuses"][status_cat]])
                except KeyError:
                    row.append(0)
            data.append(row)

        pconfig = {
            "id": "fastqc-status-check-heatmap",
            "title": "FastQC: Status Checks",
            "xlab": "Section Name",
            "ylab": "Sample",
            "min": 0,
            "max": 1,
            "square": False,
            "colstops": [
                [0, "#ffffff"],
                [0.25, "#d9534f"],
                [0.5, "#fee391"],
                [1, "#5cb85c"],
            ],
            "tt_decimals": 1,
            "legend": False,
            "display_values": False,
            "xcats_samples": False,
            "angled_xticks": False,
        }

        self.add_section(
            name="Status Checks",
            anchor="fastqc_status_checks",
            description="""
                Status for each FastQC section showing whether results seem entirely normal (green),
                slightly abnormal (orange) or very unusual (red).
            """,
            helptext="""
                FastQC assigns a status for each section of the report.
                These give a quick evaluation of whether the results of the analysis seem
                entirely normal (green), slightly abnormal (orange) or very unusual (red).

                It is important to stress that although the analysis results appear to give a pass/fail result,
                these evaluations must be taken in the context of what you expect from your library.
                A 'normal' sample as far as FastQC is concerned is random and diverse.
                Some experiments may be expected to produce libraries which are biased in particular ways.
                You should treat the summary evaluations therefore as pointers to where you should concentrate
                your attention and understand why your library may not look random and diverse.

                Specific guidance on how to interpret the output of each module can be found in the relevant
                report section, or in the [FastQC help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

                In this heatmap, we summarise all of these into a single heatmap for a quick overview.
                Note that not all FastQC sections have plots in MultiQC reports, but all status checks
                are shown in this heatmap.
            """,
            plot=heatmap.plot(data, list(status_cats.values()), s_names, pconfig),
        )

    def get_status_cols(self, section):
        """Helper function - returns a list of colours according to the FastQC
        status of this module for each sample."""
        colours = dict()
        for s_name in self.fastqc_data:
            status = self.fastqc_data[s_name]["statuses"].get(section, "default")
            colours[s_name] = self.status_colours[status]
        return colours


def _range_bp_to_num(bp: Union[str, int], method: Literal["start", "median", "mean"]) -> Union[int, float]:
    """
    Helper function - FastQC often gives base pair ranges (e.g. 10-15) which are not helpful when plotting.
    This function returns either the median or the start of the interval. If not a range, just returns the int.
    """

    if isinstance(bp, str) and "-" in bp:
        try:
            maxlen = int(bp.split("-", 1)[1])
            minlen = int(bp.split("-", 1)[0])
            if method == "start":
                return minlen
            elif method == "median":
                return ((maxlen - minlen) // 2) + minlen
            else:
                return ((maxlen + minlen) / 2) + minlen
        except TypeError:
            pass
    return int(bp)
