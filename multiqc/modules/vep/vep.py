#!/usr/bin/env python

""" MultiQC module to parse output from Samblaster """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule
import ast

# Initialise the logger
log = logging.getLogger(__name__)


def extract_vep_html_data(chart_title, html_content):
    """Function for finding and extracting VEP stats that were stored as javascript arrays"""
    found_matches = re.findall("{}.*google.visualization.*;".format(chart_title), html_content)
    if len(found_matches) > 0:
        array_content = re.search(r"\[\[.*]]", found_matches[0])
        if array_content:
            x = array_content.group(0)
            try:
                rv = ast.literal_eval(x)
                return rv
            except ValueError:
                return None
        else:
            return None
    else:
        return None


class MultiqcModule(BaseMultiqcModule):
    """VEP"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="VEP",
            anchor="vep",
            href="https://www.ensembl.org/info/docs/tools/vep/index.html",
            info="Ensembl VEP determines the effect of your variants on genes, transcripts and protein sequences, "
            "as well as regulatory regions.",
        )

        self.vep_data = dict()
        # Scan for VEP stats in plain html format
        for f in self.find_log_files("vep/vep_html", filehandles=True):
            self.parse_vep_html(f)

        # Scan for VEP stats in plain text format
        for f in self.find_log_files("vep/vep_txt", filehandles=True):
            self.parse_vep_txt(f)

        # Filter to strip out ignored sample names
        self.vep_data = self.ignore_samples(self.vep_data)

        if len(self.vep_data) > 0:
            # Add general stats table
            self.add_stats_table()
            # Generate bar graphs if any data was found
            self.generate_bar_graphs()
            log.info("Found {} VEP summaries".format(len(self.vep_data)))
        else:
            raise UserWarning

    def parse_vep_html(self, f):
        """This Function will parse VEP summary files with HTML format"""
        # Initialise vep_data dictionary for the current sample
        if f["s_name"] not in self.vep_data:
            self.vep_data[f["s_name"]] = {}
        vepf = f["f"]
        html_content = vepf.read()
        # The tables with the titles given below have common format inside the javascript section
        titles = [
            "Variant classes",
            "Consequences \(most severe\)",
            "Consequences \(all\)",
            "Coding consequences",
            "SIFT summary",
            "PolyPhen summary",
            "Variants by chromosome",
            "Position in protein",
        ]
        # Grab values of tables with the titles given above
        for title in titles:
            # Remove backslashes which were required for re
            clean_title = title.replace("\\", "")
            extracted_data = extract_vep_html_data(title, html_content)
            if extracted_data:
                if clean_title not in self.vep_data[f["s_name"]]:
                    self.vep_data[f["s_name"]][clean_title] = {}
                plot_data = {}
                for i in range(1, len(extracted_data)):
                    plot_data[extracted_data[i][0]] = extracted_data[i][1]
                self.vep_data[f["s_name"]][clean_title] = plot_data

        # Finally parse general stats table which is an HTML table
        title = "General statistics"
        hits = re.search(r'General statistics</h3><table class="stats_table">.*?</table>', html_content)
        if hits:
            if title not in self.vep_data[f["s_name"]]:
                self.vep_data[f["s_name"]][title] = {}
            hit = hits.group(0)
            rows = re.findall(r"<tr>.*?</tr>", hit)
            for r in rows:
                cells = re.findall(r"<td>.*?</td>", r)
                if len(cells) == 2:
                    key = cells[0][4:-5]
                    value = cells[1][4:-5]
                    if key == "Novel / existing variants":
                        values = value.split("/")
                        novel = values[0].split("(")[0].replace(" ", "")
                        existing = values[1].split("(")[0].replace(" ", "")
                        self.vep_data[f["s_name"]][title]["Novel variants"] = int(novel)
                        self.vep_data[f["s_name"]][title]["Existing variants"] = int(existing)
                    else:
                        self.vep_data[f["s_name"]][title][key] = int(value)

    def parse_vep_txt(self, f):
        """This Function will parse VEP summary files with plain text format"""
        # Initialise vep_data dictionary for the current sample
        if f["s_name"] not in self.vep_data:
            self.vep_data[f["s_name"]] = {}
        vepf = f["f"]
        txt_data = {}
        title = ""
        for line in vepf.readlines():
            line = line.rstrip("\n")
            if len(line) < 1:
                continue
            if line[0] == "[":
                title = line.replace("[", "").replace("]", "")
                txt_data[title] = {}
                continue
            key, value = line.split("\t")
            if key == "Novel / existing variants":
                values = value.split("/")
                novel = values[0].split("(")[0].replace(" ", "")
                existing = values[1].split("(")[0].replace(" ", "")
                txt_data[title]["Novel variants"] = int(novel)
                txt_data[title]["Existing variants"] = int(existing)
            elif title != "VEP run statistics":
                txt_data[title][key] = int(value)
            else:
                txt_data[title][key] = value

        for title in txt_data.keys():
            if title not in self.vep_data[f["s_name"]]:
                self.vep_data[f["s_name"]][title] = {}
            self.vep_data[f["s_name"]][title] = txt_data[title]

    def add_stats_table(self):
        """Add a section with VEP General Statistics"""
        table_config = {
            "id": "vep-general-stats",  # ID used for the table
            "table_title": "VEP General Statistics",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "sortRows": True,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "Sample Name",  # The header used for the first column
            "no_beeswarm": True,  # Force a table to always be plotted (beeswarm by default if many rows)
        }
        title = "General statistics"
        plot_cats = OrderedDict()
        plot_data = OrderedDict()
        color_list = ["Oranges", "Blues", "Reds", "Greens"]
        for sample_name in self.vep_data:
            if self.vep_data[sample_name][title] and len(self.vep_data[sample_name][title]) > 0:
                order = 0
                plot_data[sample_name] = {}
                for key in self.vep_data[sample_name][title]:
                    order += 1
                    if key not in plot_cats:
                        plot_cats[key] = {"name": key, "format": "{:,.0f}", "scale": color_list[order % 4]}
                    plot_data[sample_name][key] = self.vep_data[sample_name][title][key]
        if len(plot_data) > 0:
            self.add_section(
                name="General Statistics",
                anchor="vep-general-statistics",
                description="VEP General statistics",
                helptext="Table showing general statistics of VEP annotaion run",
                plot=table.plot(plot_data, plot_cats, table_config),
            )

    def generate_bar_graphs(self):
        """Add bar graphs from the given VEP stats"""
        titles = [
            "Variant classes",
            "Consequences (most severe)",
            "Consequences (all)",
            "Coding consequences",
            "SIFT summary",
            "PolyPhen summary",
            "Variants by chromosome",
            "Position in protein",
        ]
        for title in titles:
            htmlid = re.sub("\W*", "", title)
            plotid = htmlid + "plot"
            plot_config = {"id": plotid, "title": "VEP: {}".format(title), "ylab": "Counts", "stacking": "normal"}
            plot_cats = OrderedDict()
            plot_data = OrderedDict()
            for sample_name in self.vep_data:
                if self.vep_data[sample_name][title] and len(self.vep_data[sample_name][title]) > 0:
                    plot_data[sample_name] = {}
                    for key in self.vep_data[sample_name][title]:
                        if key not in plot_cats:
                            plot_cats[key] = {"name": key}
                        plot_data[sample_name][key] = self.vep_data[sample_name][title][key]

            if len(plot_data) > 0:
                self.add_section(
                    name=title,
                    anchor=htmlid,
                    description='Bar graph showing VEP "{}" statistics'.format(title),
                    helptext='Stacked bar graph showing VEP "{}" statistics'.format(title),
                    plot=bargraph.plot(plot_data, cats=plot_cats, pconfig=plot_config),
                )
