'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers vc_hethom_ratio_metrics data and prepares it for the output report.

It relies on the following official source:
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/QCMetricsCoverageReports.htm?Highlight=Het%2FHom
'''

import logging
import re
from collections import OrderedDict, defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph
from multiqc.utils.util_functions import write_data_file

from .utils import check_duplicate_samples, make_log_report

# Initialise the logger.
log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN VC Het/Hom ratio"


""" These 2 variables are to control html section's appearance. """
BARGRAPHS_ENABLED = False  # True to include bargraphs.

# True to show only 1 plot. Others can be chosen from drop-down list.
# Common plots' bugs associated with resizing/reloading are present if True.
# More about bugs in make_html_for_own_section.
SQUASH_SECTIONS = True


class DragenVCHetHomRatioMetrics(BaseMultiqcModule):
    """Public members of the DragenVCHetHomRatioMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined."""

    def add_vc_hethom_ratio_metrics(self):
        """The main function of the dragen hethom metrics module.
        Public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns keys with sample names
        or an empty set if no samples were found or all samples were ignored."""

        hethom_data = {}

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones.
        all_samples = defaultdict(list)

        # Metric IDs with original consistent metric strings.
        all_metrics = {}

        for file in self.find_log_files("dragen/vc_hethom_ratio_metrics"):
            out = hethom_parser(file)
            if out["success"]:
                self.add_data_source(file, section="stats")
                cleaned_sample = self.clean_s_name(out["sample_name"], file)
                all_samples[cleaned_sample].append(file)
                hethom_data[cleaned_sample] = out["data"]  # Add/overwrite the sample.
                all_metrics.update(out["metric_IDs"])

        # Filter to strip out ignored sample names.
        hethom_data = self.ignore_samples(hethom_data)
        if not hethom_data:
            return set()

        check_duplicate_samples(all_samples, log, "dragen/vc_hethom_ratio_metrics")

        # Write data to file.
        txt_data = make_txt_data(hethom_data, all_metrics)
        for section in txt_data:
            self.write_data_file(txt_data[section], section)

        # Currently everything is placed in a single html block.
        html = make_html_for_own_section(hethom_data)
        self.add_section(
            name="VC Het/Hom Ratio Metrics",
            anchor="dragen-vc-hethom-ratio-metrics-metrics",
            description="When the germline small variant caller is executed, "
            "DRAGEN calculates a per het/hom ratio per contig. DRAGEN reports the "
            "ratios for both the raw (PREFILTER) and hard-filtered (POSTFILTER) VCF. "
            "Press the `Help` button for details.",
            helptext="The metrics are output to the &#60;output-file-prefix&#62;.vc_hethom_ratio_metrics.csv file. "
            "The file contains the following values for each primary contig processed:\n\n"
            "* Contig\n\n"
            "* Number of heterozygous variants\n\n"
            "* Number of homozygous variants\n\n"
            "* Het/Hom ratio\n\n"
            "Notes:\n\n"
            "* If there is no point in the linegraph for some sample, then it is probably "
            "because its value is inf. It shall happen only for the 'Y' chromosome.\n\n"
            "* Value of -1 in the linegraph means that metric could not be extracted.\n\n"
            "* If SQUASH_SECTIONS is True, then resizing/reloading can break tables.",
            content=html,
        )
        # Report found info/warnings/errors. You can disable it anytime, if it is not wanted.
        make_log_report(log_data, log, "vc_hethom_ratio_metrics")

        return hethom_data.keys()


def make_txt_data(data, metric_IDs):
    """Reorganize data for storing in .txt output."""
    out_data = defaultdict(lambda: defaultdict(dict))
    for sample in data:
        for section in data[sample]:
            for contig in data[sample][section]:
                for metric in data[sample][section][contig]:
                    # Make monolithic sample and section names by replacing any
                    # sequence of spaces/hyphens/dots/underscores by single underscore.
                    new_sample = re.sub("(\s|-|\.|_)+", "_", sample + "_" + contig)
                    new_section = "dragen_vc_hethom_ratio_metrics_" + re.sub("(\s|-|\.|_)+", "_", section)
                    metric_id = metric_IDs[metric] if metric in metric_IDs else metric
                    out_data[new_section][new_sample][metric_id] = data[sample][section][contig][metric]
    return out_data


def make_html_for_own_section(DATA):
    """Create HTML block ready to be inserted into the output report."""

    # The CHROMOSOMES is a list with 1, 2,..., 22, x, y.
    # Must be in lower case, because of how make_consistent_metric works.
    CHROMOSOMES = [str(n) for n in range(1, 23)] + ["x", "y"]

    ## Prepare data and config for linegraphs.

    linegraph_data = defaultdict(lambda: defaultdict(dict))

    # The sole purpose of the chromosome_and_metric array is to guarantee the sequence of metrics.
    # Creating the output data by iterating through the input data would result in wrong output if
    # the first catched sample has a sequence of metrics, which is distinct from the sequence of
    # the same metrics in other samples.
    chromosome_and_metric = [(CHR.upper(), CHR + " het/hom ratio") for CHR in CHROMOSOMES]
    for sample in DATA:
        for section in DATA[sample]:
            for contig in DATA[sample][section]:
                metrics = DATA[sample][section][contig]
                new_sample = sample + "_" + contig
                for chrom, metric in chromosome_and_metric:
                    # Check if metric is in the sample and its value is a number.
                    # MultiQC can crash without the latter.
                    if metric in metrics and isinstance(metrics[metric], (float, int)):
                        linegraph_data[section][new_sample][chrom] = metrics[metric]
                    # Additional safety measure to prevent producing broken output
                    # in case the first catched sample does not contain all metrics.
                    else:
                        linegraph_data[section][new_sample][chrom] = -1

    LINEGRAPH_CONFIG = {
        "id": "dragen_vc_hethom_ratio_metrics_linegraph_",
        "title": "Dragen: Het/Hom ratio. Section - {}",
        "ylab": "Het/Hom ratio",
        "xlab": "Chromosome",
        "categories": True,  # Set to True to use x values as categories instead of numbers.
        "tt_label": "{point.category} Het/Hom ratio: {point.y}%",  # Use to customise tooltip label
    }

    ## Prepare data, cats and config for bargraphs.

    if BARGRAPHS_ENABLED:
        bargraph_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        het_metrics = [CHR + " heterozygous" for CHR in CHROMOSOMES]
        hom_metrics = [CHR + " homozygous" for CHR in CHROMOSOMES]
        for sample in DATA:
            for section in DATA[sample]:
                for contig in DATA[sample][section]:
                    metrics = DATA[sample][section][contig]
                    new_sample = sample + "_" + contig
                    for metric in metrics:
                        if metric in het_metrics:
                            bargraph_data[section]["het"][new_sample][metric] = metrics[metric]
                        elif metric in hom_metrics:
                            bargraph_data[section]["hom"][new_sample][metric] = metrics[metric]
                        # Other metrics (potentially may be added in the future) are not supported.
                        else:
                            pass

        CATS = [OrderedDict(), OrderedDict()]
        for CHR in CHROMOSOMES:
            CATS[0][CHR + " heterozygous"] = {"name": CHR.upper() + " Het"}
            CATS[1][CHR + " homozygous"] = {"name": CHR.upper() + " Hom"}

        BARGRAPH_CONFIG = {
            "id": "dragen_vc_hethom_ratio_metrics_bargraph_",
            "title": "Dragen: Het/Hom Metrics. Section - {}",
            "tt_decimals": 0,
            "ylab": "",  # to pass lint test. ylabs are defined in data_labels.
            "data_labels": [
                {
                    "name": "Het",
                    "ylab": "Number of heterozygous variants",
                },
                {
                    "name": "Hom",
                    "ylab": "Number of homozygous variants",
                },
            ],
        }

    ## Now build the html.

    html = ""

    BASE_PREFIX_ID = "DRAGEN_hethom_"
    BASE_SECTION_ID = BASE_PREFIX_ID + "inner_section_"
    SELECT_TAG_ID = BASE_PREFIX_ID + "select_tag"

    # Create a set of sections and convert it to list to guarantee the order.
    all_sections = list({section for sample in DATA for section in DATA[sample]})

    """ Plots' responsiveness bug.
    Resizing/reloading can affect table's attributes (eg width, height).
    The following was noticed:
    - size changes (eg can become very large or just wider).
    - collapsing with other blocks (eg own scroll, left navigation bar, other sections)
    Changing the 'display' attribute can cause such problems. The root of the issue is unknown.
    It might be some strange combination of some working/logical aspects of Bootstrap, Highcharts and
    maybe browser engines. It seems to be the common problem for plots (eg dragen_fastqc has the same issues).
    """
    if SQUASH_SECTIONS:
        ## First add the drop-down list with sections.
        html += "<label for='{}'>Choose a section:</label>".format(SELECT_TAG_ID)
        # "display" attr is set to block, because select is an inline element.
        html += "<select id='{}' style='display:block'>".format(SELECT_TAG_ID)

        selected_section = all_sections[0]
        html += "<option selected='selected' value='{sec}'>{sec}</option>".format(sec=selected_section)

        # Now add all other sections as options.
        for section in all_sections[1:]:
            html += "<option value='{sec}'>{sec}</option>".format(sec=section)
        html += "</select>"

    ## Append sections as div blocks.
    BASE_DIV_STRUCTURE = "<div id='{}' class='mqc-section-plot'>"
    EMPTY_LINEGRAPH_TAG = (
        "<p style='text-align:center;'>"
        "Linegraph could not be created because 'Chromosome Het/Hom ratio' metrics could not be found"
        "</p>"
    )
    EMPTY_BARGRAPH_TAG = (
        "<p style='text-align:center;'>"
        "Bargraph could not be created because 'Chromosome Heterozygous' and 'Chromosome Homozygous' metrics could not be found"
        "</p>"
    )
    for section in all_sections:
        section_id = re.sub("\s+", "_", section)
        html += BASE_DIV_STRUCTURE.format(BASE_SECTION_ID + section_id)

        if section in linegraph_data:
            line_config = LINEGRAPH_CONFIG.copy()
            line_config["id"] += section_id
            line_config["title"] = line_config["title"].format(section)
            html += linegraph.plot(linegraph_data[section], line_config)

        # Extremely unlikely to happen, but nevertheless shall be checked and handled appropriately.
        else:
            html += EMPTY_LINEGRAPH_TAG

        if BARGRAPHS_ENABLED:
            html += "<br>"  # To separate both plots from each other.

            if section in bargraph_data:
                bar_config = BARGRAPH_CONFIG.copy()
                bar_config["id"] += section_id
                bar_config["title"] = bar_config["title"].format(section)
                # There is a tiny chance that het and/or hom is not present in data. Not checking it would crash the program.
                het_data = bargraph_data[section]["het"] if "het" in bargraph_data[section] else {}
                hom_data = bargraph_data[section]["hom"] if "hom" in bargraph_data[section] else {}
                if het_data or hom_data:
                    html += bargraph.plot([het_data, hom_data], CATS, bar_config)
                else:
                    html += EMPTY_BARGRAPH_TAG
            # Extremely unlikely to happen, but nevertheless shall be checked and handled appropriately.
            else:
                html += EMPTY_BARGRAPH_TAG

        html += "</div>"  # "Close" the current section.

    if SQUASH_SECTIONS:
        ## Now we can add JS.
        # Linefeeds have only 1 purpose: make the following block readable in output html.
        # It may be useful for checking results, otherwise it is just a long line.
        # setTimeout seems to improve the mentioned bugs, but does not solve them.
        # hethom_sections_display_values stores original display values to restore them properly.
        html += (
            "<script>\n"
            + "window.addEventListener('DOMContentLoaded', () => {\n"
            + "  let select_tag_ref = document.getElementById('{}');\n".format(SELECT_TAG_ID)
            + "  let chosen_section_id = select_tag_ref.value.replaceAll(' ','_');\n"
            + "  const BASE_SECTION_ID = '{}';\n".format(BASE_SECTION_ID)
            + "  function change_section(){\n"
            + "    let new_chosen_section_id = select_tag_ref.value.replaceAll(' ','_');\n"
            + "    let current_plot = document.getElementById(BASE_SECTION_ID + chosen_section_id);\n"
            + "    let new_plot = document.getElementById(BASE_SECTION_ID + new_chosen_section_id);\n"
            + "    current_plot.style.display = 'none';\n"
            + "    new_plot.style.display = window.hethom_sections_display_values[BASE_SECTION_ID + new_chosen_section_id];\n"
            + "    chosen_section_id = new_chosen_section_id;\n"
            + "  };\n"
            + "  select_tag_ref.addEventListener('change', change_section);\n"
            + "});\n"
            + "window.addEventListener('load', () => {\n"
            + "  setTimeout(() => {\n"
            + "    const BASE_SECTION_ID = '{}';\n".format(BASE_SECTION_ID)
            + "    window.hethom_sections_display_values = {};\n"
            + "    let section_id;\n"
        )
        ignore_first_section = True
        for section in all_sections:
            section_id = re.sub("\s+", "_", section)
            html += (
                "section_id = BASE_SECTION_ID+'{}';\n".format(section_id)
                + "hethom_sections_display_values[section_id] = document.getElementById(section_id).style.display;\n"
            )
            if ignore_first_section:
                ignore_first_section = False
            else:
                html += "document.getElementById(section_id).style.display='none';\n"

        html += "},1000)});\n"
        html += "</script>"  # "Close" script block.

    return html


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The log_data is the container for found info, warnings and errors.
Collected data is stored in groups to make the output log file more readable.

warnings:
- invalid_file_names, which do not conform to the:
    <output-file-prefix>.vc_hethom_ratio_metrics.csv

debug/info:
- invalid_file_lines, which do not conform to:
  <section>,<contig>,<metric>,<value>
- unknown_metrics are those, which are not present in the WELL_KNOWN_METRICS.
  They are not included in the html report.
- unusual_values are those except for int, float and inf.
'''
# "1,2,3...,22,X,Y Heterozygous / Homozygous / Het/Hom ratio" metrics
WELL_KNOWN_METRICS = [
    str(chrom) + het_hom_rat
    for chrom in list(range(23)) + ["x", "y"]
    for het_hom_rat in [" heterozygous", " homozygous", " het/hom ratio"]
]
log_data = {
    "invalid_file_names": defaultdict(list),
    "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
    "unknown_metrics": set(),
    "unusual_values": defaultdict(lambda: defaultdict(dict)),
}


def create_parser():
    """Isolation for the parsing codeblock. Returns the parse_vc_hethom_ratio_metrics_file."""

    def make_consistent_metric(metric):
        """Tries to guarantee consistency to a certain degree.
        Metric-specific peculiarities can be handled here."""
        # metric will be stored in .txt output.
        metric = re.sub("\s+", " ", metric).strip()
        # metric_id will be used for internal storage.
        metric_id = metric.lower()
        return metric, metric_id

    def make_consistent_section(section):
        """Tries to guarantee consistency for section string."""
        return re.sub("\s+", " ", section).strip().upper()

    # The official/accepted structure of the input file:
    # <output-file-prefix>.vc_hethom_ratio_metrics.csv
    FILE_RGX = re.compile("(.+)\.vc_hethom_ratio_metrics\.csv$")

    # The following regex is used to check input lines.
    LINE_RGX = re.compile("^([^,]+),([^,]+),([^,]+),([^,]+)$")

    def parse_vc_hethom_ratio_metrics_file(file_handler):
        """
        Parser for *.vc_hethom_ratio_metrics.csv files.

        Input:  file_handler with necessary info - file name/content/root.

        Output: {"success": False} if file name test failed. Otherwise:
                {"success": False/True,
                 "sample_name": <output-file-prefix>,
                 "data": {section: {contig: {metric: value}}},
                 "metric_IDs": {metric_ID: original_metric}
                }, with success False if all file's lines are invalid.
        """
        file, root = file_handler["fn"], file_handler["root"]

        file_match = FILE_RGX.search(file)
        if not file_match:
            log_data["invalid_file_names"][root].append(file)
            return {"success": False}

        sample = file_match.group(1)

        success = False
        data = defaultdict(lambda: defaultdict(dict))
        metric_IDs = {}

        for line in file_handler["f"].splitlines():
            # Check the general line structure.
            line_match = LINE_RGX.search(line)

            # If line is fine then extract the necessary fields.
            if line_match:
                success = True
                section, contig, metric, value = line_match.groups()

            # Otherwise check if line is empty. If not then report it. Go to the next line.
            else:
                if not re.search("^\s*$", line):
                    log_data["invalid_file_lines"][root][file].append(line)
                continue

            # Check the extracted value. Shall be int, float or inf.
            # inf is catched by the second try, since float('inf') produces inf
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    log_data["unusual_values"][root][file][metric] = value

            # Modify extracted parts to improve consistency.
            consistent_section = make_consistent_section(section)
            consistent_metric, metric_id = make_consistent_metric(metric)

            if metric_id not in WELL_KNOWN_METRICS:
                log_data["unknown_metrics"].add(metric)

            metric_IDs[metric_id] = consistent_metric
            data[consistent_section][contig][metric_id] = value

        # Return results of parsing the input file.
        return {
            "success": success,
            "sample_name": sample,
            "data": data,
            "metric_IDs": metric_IDs,
        }

    # Return the parser closure.
    return parse_vc_hethom_ratio_metrics_file


hethom_parser = create_parser()
