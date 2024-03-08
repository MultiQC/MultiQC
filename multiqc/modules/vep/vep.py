""" MultiQC module to parse output from VEP """


import ast
import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.utils import mqc_colour

# Initialise the logger
log = logging.getLogger(__name__)


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
            doi="10.1186/s13059-016-0974-4",
        )

        self.vep_data = dict()

        # Scan for VEP stats in plain html format
        for f in self.find_log_files("vep/vep_html", filehandles=True):
            self.parse_vep_html(f)
            self.add_data_source(f)

        # Scan for VEP stats in plain text format
        for f in self.find_log_files("vep/vep_txt", filehandles=True):
            self.parse_vep_txt(f)
            self.add_data_source(f)

        # Add version information
        for sample, data in self.vep_data.items():
            if "VEP run statistics" not in data:
                continue

            vep_version, api_version = data["VEP run statistics"]["VEP version (API)"].strip().split(" ")
            api_version = api_version.replace("(", "").replace(")", "")
            self.add_software_version(vep_version, sample)
            # Only add API version if it's different to VEP version
            if vep_version != api_version:
                self.add_software_version(api_version, sample, "VEP API")
        # Filter to strip out ignored sample names
        self.vep_data = self.ignore_samples(self.vep_data)

        # Stop if we didn't get any samples
        if len(self.vep_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.vep_data)} VEP summaries")

        # Write data to file
        self.write_data_file(self.vep_data, "vep")

        # Add general stats table
        self.add_stats_table()

        # Bar graphs
        self.bar_graph_variant_classes()
        self.bar_graph_consequences()
        self.bar_graph_sift()
        self.bar_graph_polyphen()
        self.bar_graph_variants_by_chromosome()
        self.bar_graph_position_in_protein()

    def extract_vep_html_data(self, chart_title, html_content):
        """Function for finding and extracting VEP stats that were stored as javascript arrays"""
        found_matches = re.findall(f"{chart_title}.*google.visualization.*;", html_content)
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
            r"Consequences \(most severe\)",
            r"Consequences \(all\)",
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
            extracted_data = self.extract_vep_html_data(title, html_content)
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
                    if value == "-":
                        continue
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
            if value == "-":
                continue
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
            "id": "vep-general-stats",
            "namespace": "VEP",
            "table_title": "VEP General Statistics",
        }
        table_data = {s_name: self.vep_data[s_name]["General statistics"] for s_name in self.vep_data}

        # Build the header configs
        cat_names = [
            "Overlapped regulatory features",
            "Overlapped transcripts",
            "Overlapped genes",
            "Existing variants",
            "Novel variants",
            "Variants filtered out",
            "Variants processed",
            "Lines of input read",
        ]
        # Set up the base config for each column
        table_cats = dict()
        color_list = ["Oranges", "Reds", "Blues", "Greens"]
        for order, header in enumerate(cat_names):
            table_cats[header] = {
                "name": header,
                "format": "{:,.0f}",
                "scale": color_list[order % 4],
            }
        # Column-specific customisation
        table_cats["Lines of input read"]["hidden"] = True
        table_cats["Variants processed"]["shared_key"] = "variants"
        table_cats["Variants filtered out"]["shared_key"] = "variants"
        table_cats["Novel variants"]["shared_key"] = "variants"
        table_cats["Existing variants"]["shared_key"] = "variants"

        self.add_section(
            name="General Statistics",
            anchor="vep-general-statistics",
            helptext="Table showing general statistics of VEP annotaion run",
            plot=table.plot(table_data, table_cats, table_config),
        )

    def bar_graph_variant_classes(self):
        title = "Variant classes"
        plot_data, plot_cats, plot_config = self._prep_bar_graph(title)
        htmlid = re.sub(r"\W+", "_", title).lower()
        if len(plot_data) == 0:
            return

        plot_config["ylab"] = "Number of variants"
        self.add_section(
            name=title,
            anchor=htmlid,
            description="Classes of variants found in the data.",
            plot=bargraph.plot(plot_data, plot_cats, plot_config),
        )

    def bar_graph_consequences(self):
        p_data = []
        p_cats = []
        p_config = {"data_labels": []}
        for title in ["Consequences (most severe)", "Coding consequences", "Consequences (all)"]:
            plot_data, plot_cats, plot_config = self._prep_bar_graph(title)
            p_data.append(plot_data)
            p_cats.append(plot_cats)
            p_config.update(plot_config)
            p_config["data_labels"].append(title)
        p_config["title"] = "VEP: Variant Consequences"
        p_config["ylab"] = p_config["data_labels"][0]

        if len(plot_data) == 0 or max([len(d) for d in plot_data]) == 0:
            return

        self.add_section(
            name="Consequences",
            anchor="consequences",
            description="Predicted consequences of variations.",
            plot=bargraph.plot(p_data, p_cats, p_config),
        )

    def bar_graph_sift(self):
        title = "SIFT summary"
        plot_data, plot_cats, plot_config = self._prep_bar_graph(title)
        htmlid = re.sub(r"\W+", "_", title).lower()
        if len(plot_data) == 0:
            return

        # Customise order and colours of categories
        p_cats = {
            "tolerated": {"color": "#59ae61"},
            "tolerated_low_confidence": {"color": "#a6db9f"},
            "deleterious_low_confidence": {"color": "#fec44f"},
            "deleterious": {"color": "#d53e4f"},
        }
        for c, cat in plot_cats.items():
            if c not in p_cats:
                p_cats[c] = cat
            p_cats[c].update(cat)
            p_cats[c]["name"] = p_cats[c]["name"].replace("_", " ").capitalize()

        plot_config["ylab"] = "Number of variants"

        self.add_section(
            name=title,
            anchor=htmlid,
            description="SIFT variant effect prediction.",
            plot=bargraph.plot(plot_data, p_cats, plot_config),
        )

    def bar_graph_polyphen(self):
        title = "PolyPhen summary"
        plot_data, plot_cats, plot_config = self._prep_bar_graph(title)
        htmlid = re.sub(r"\W+", "_", title).lower()
        if len(plot_data) == 0:
            return

        # Customise order and colours of categories
        p_cats = {
            "benign": {"color": "#a6db9f"},
            "possibly_damaging": {"color": "#fec44f"},
            "probably_damaging": {"color": "#d53e4f"},
            "unknown": {"color": "#d9d9d9"},
        }
        for c, cat in plot_cats.items():
            if c not in p_cats:
                p_cats[c] = cat
            p_cats[c].update(cat)
            p_cats[c]["name"] = p_cats[c]["name"].replace("_", " ").capitalize()

        plot_config["ylab"] = "Number of variants"

        self.add_section(
            name=title,
            anchor=htmlid,
            description="PolyPhen variant effect prediction.",
            plot=bargraph.plot(plot_data, p_cats, plot_config),
        )

    def bar_graph_variants_by_chromosome(self):
        title = "Variants by chromosome"
        plot_data, plot_cats, plot_config = self._prep_bar_graph(title)
        htmlid = re.sub(r"\W+", "_", title).lower()
        if len(plot_data) == 0:
            return

        # Sort the chromosomes numerically (almost - feel free to improve)
        chrs = {chr: k["name"].replace("chr", "").split("_")[0].rjust(20, "0") for chr, k in plot_cats.items()}
        p_cats = {}
        for chr in sorted(chrs, key=chrs.get):
            p_cats[chr] = plot_cats[chr]

        plot_config["cpswitch_c_active"] = False

        self.add_section(
            name=title,
            anchor=htmlid,
            description="Number of variants found on each chromosome.",
            plot=bargraph.plot(plot_data, p_cats, plot_config),
        )

    def bar_graph_position_in_protein(self):
        title = "Position in protein"
        plot_data, plot_cats, plot_config = self._prep_bar_graph(title)
        htmlid = re.sub(r"\W+", "_", title).lower()
        if len(plot_data) == 0:
            return

        # Nice graduated colours for categories
        colour_scale = mqc_colour.mqc_colour_scale("YlGnBu", 0, len(plot_cats) - 1)
        for idx, k in enumerate(plot_cats):
            plot_cats[k]["color"] = colour_scale.get_colour(idx, lighten=0.8)

        plot_config["cpswitch_c_active"] = False

        self.add_section(
            name=title,
            anchor=htmlid,
            description="Relative position of affected amino acids in protein.",
            plot=bargraph.plot(plot_data, plot_cats, plot_config),
        )

    def _prep_bar_graph(self, title):
        plot_data = dict()
        for s_name in self.vep_data:
            if title in self.vep_data[s_name]:
                plot_data[s_name] = self.vep_data[s_name][title]
        plot_cats = dict()
        htmlid = re.sub(r"\W+", "_", title).lower()
        plotid = f"{htmlid}_plot"
        plot_config = {
            "id": plotid,
            "title": f"VEP: {title}",
            "ylab": "Number of variants",
        }
        if len(plot_data) == 0:
            return [plot_data, plot_cats, plot_config]

        # Build plot categories
        for d in plot_data.values():
            for k in d.keys():
                if k not in plot_cats:
                    plot_cats[k] = {"name": k}

        return [plot_data, plot_cats, plot_config]
