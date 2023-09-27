""" MultiQC module to parse output from FastQ Screen """


import json
import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph
from multiqc.utils import report

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"Fastq_screen version: ([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FastQ Screen",
            anchor="fastq_screen",
            href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/",
            info="allows you to screen a library of sequences in FastQ format against"
            " a set of sequence databases so you can see if the composition of the"
            " library matches with what you expect.",
            doi="10.12688/f1000research.15931.2",
        )

        # Find and load any FastQ Screen reports
        self.fq_screen_data = dict()
        self.num_orgs = 0
        for f in self.find_log_files("fastq_screen", filehandles=True):
            parsed_data = self.parse_fqscreen(f)
            if parsed_data is not None:
                if f["s_name"] in self.fq_screen_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.add_data_source(f)
                self.fq_screen_data[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.fq_screen_data = self.ignore_samples(self.fq_screen_data)

        if len(self.fq_screen_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.fq_screen_data)))

        # Check whether we have a consistent number of organisms across all samples
        num_orgs = set([len(orgs) for orgs in self.fq_screen_data.values()])

        # Section 1 - Alignment Profiles
        # Posh plot only works for around 20 samples, 8 organisms. If all samples have the same number of organisms.
        if (
            len(num_orgs) == 1
            and len(self.fq_screen_data) * self.num_orgs <= 160
            and not config.plots_force_flat
            and not getattr(config, "fastqscreen_simpleplot", False)
        ):
            self.add_section(name="Mapped Reads", anchor="fastq_screen_mapped_reads", content=self.fqscreen_plot())
        # Use simpler plot that works with many samples
        else:
            self.add_section(name="Mapped Reads", anchor="fastq_screen_mapped_reads", plot=self.fqscreen_simple_plot())

        # Section 2 - Optional bisfulfite plot
        self.fqscreen_bisulfite_plot()

        # Write the total counts and percentages to files
        self.write_data_file(self.parse_csv(), "multiqc_fastq_screen")

    def parse_fqscreen(self, f):
        """Parse the FastQ Screen output into a 3D dict"""
        parsed_data = OrderedDict()
        nohits_pct = None
        headers = None
        bs_headers = None
        for l in f["f"]:
            # Skip comment lines
            if l.startswith("#"):
                version_match = re.search(VERSION_REGEX, l)
                if version_match:
                    self.add_software_version(version_match.group(1), f["s_name"])
                continue
            if l.startswith("%Hit_no_genomes:") or l.startswith("%Hit_no_libraries:"):
                nohits_pct = float(l.split(":", 1)[1])
                parsed_data["No hits"] = {"percentages": {"one_hit_one_library": nohits_pct}}
            else:
                s = l.strip().split("\t")

                # Regular FastQ Screen table section
                if len(s) == 12:
                    if headers is None:
                        headers = s
                    else:
                        # Can't use #Reads in subset as varies. #Reads_processed should be same for all orgs in a sample
                        parsed_data["total_reads"] = int(s[1])
                        # Loop through all columns
                        parsed_data[s[0]] = {"percentages": {}, "counts": {}}
                        for idx, h in enumerate(headers):
                            if idx == 0:
                                continue
                            dtype = "percentages" if h.startswith("%") else "counts"
                            field = (
                                h.replace("%", "")
                                .replace("#", "")
                                .replace("genomes", "libraries")
                                .replace("genome", "library")
                                .lower()
                            )
                            parsed_data[s[0]][dtype][field] = float(s[idx])

                # Optional Bisulfite table section
                elif len(s) == 9:
                    if bs_headers is None:
                        bs_headers = s
                    else:
                        # Loop through all columns
                        parsed_data[s[0]]["bisulfite_percentages"] = {}
                        parsed_data[s[0]]["bisulfite_counts"] = {}
                        for idx, h in enumerate(bs_headers):
                            if idx == 0:
                                continue
                            dtype = "bisulfite_percentages" if h.startswith("%") else "bisulfite_counts"
                            field = h.replace("%", "").replace("#", "").lower()
                            parsed_data[s[0]][dtype][field] = float(s[idx])

        if len(parsed_data) == 0:
            return None

        # Calculate no hits counts
        if parsed_data.get("total_reads") is not None and nohits_pct is not None:
            parsed_data["No hits"]["counts"] = {
                "one_hit_one_library": int((nohits_pct / 100.0) * float(parsed_data["total_reads"]))
            }
        else:
            log.warning("Couldn't find number of reads with no hits for '{}'".format(f["s_name"]))

        self.num_orgs = max(len(parsed_data), self.num_orgs)
        return parsed_data

    def parse_csv(self):
        totals = OrderedDict()
        for s in sorted(self.fq_screen_data.keys()):
            totals[s] = OrderedDict()
            for org in self.fq_screen_data[s]:
                if org == "total_reads":
                    totals[s]["total_reads"] = self.fq_screen_data[s][org]
                    continue
                try:
                    k = "{} counts".format(org)
                    totals[s][k] = self.fq_screen_data[s][org]["counts"]["one_hit_one_library"]
                    totals[s][k] += self.fq_screen_data[s][org]["counts"].get("multiple_hits_one_library", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["counts"].get("one_hit_multiple_libraries", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["counts"].get("multiple_hits_multiple_libraries", 0)
                except KeyError:
                    pass
                try:
                    k = "{} percentage".format(org)
                    totals[s][k] = self.fq_screen_data[s][org]["percentages"]["one_hit_one_library"]
                    totals[s][k] += self.fq_screen_data[s][org]["percentages"].get("multiple_hits_one_library", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["percentages"].get("one_hit_multiple_libraries", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["percentages"].get(
                        "multiple_hits_multiple_libraries", 0
                    )
                except KeyError:
                    pass
        return totals

    def fqscreen_plot(self):
        """Makes a fancy custom plot which replicates the plot seen in the main
        FastQ Screen program. Not useful if lots of samples as gets too wide."""

        categories = list()
        getCats = True
        data = list()
        p_types = OrderedDict()
        p_types["multiple_hits_multiple_libraries"] = {"col": "#7f0000", "name": "Multiple Hits, Multiple Genomes"}
        p_types["one_hit_multiple_libraries"] = {"col": "#ff0000", "name": "One Hit, Multiple Genomes"}
        p_types["multiple_hits_one_library"] = {"col": "#00007f", "name": "Multiple Hits, One Genome"}
        p_types["one_hit_one_library"] = {"col": "#0000ff", "name": "One Hit, One Genome"}
        for k, t in p_types.items():
            first = True
            for s in sorted(self.fq_screen_data.keys()):
                thisdata = list()
                if len(categories) > 0:
                    getCats = False
                for org in sorted(self.fq_screen_data[s]):
                    if org == "total_reads":
                        continue
                    try:
                        thisdata.append(self.fq_screen_data[s][org]["percentages"][k])
                    except KeyError:
                        thisdata.append(None)
                    if getCats:
                        categories.append(org)
                td = {"name": t["name"], "stack": s, "data": thisdata, "color": t["col"]}
                if first:
                    first = False
                else:
                    td["linkedTo"] = ":previous"
                data.append(td)

        plot_id = report.save_htmlid("fq_screen_plot")
        html = """<div id={plot_id} class="fq_screen_plot hc-plot"></div>
        <script type="application/json" class="fq_screen_dict">{dict}</script>
        """.format(
            plot_id=json.dumps(plot_id),
            dict=json.dumps({"plot_id": plot_id, "data": data, "categories": categories}),
        )

        html += """<script type="text/javascript">
            fq_screen_dict = { }; // { <plot_id>: data, categories }
            $('.fq_screen_dict').each(function (i, elem) {
                var dict = JSON.parse(elem.innerHTML);
                fq_screen_dict[dict.plot_id] = dict;
            });

            $(function () {
                // In case of repeated modules: #fq_screen_plot, #fq_screen_plot-1, ..
                $(".fq_screen_plot").each(function () {
                    var plot_id = $(this).attr('id');

                    $(this).highcharts({
                        chart: { type: "column", backgroundColor: null },
                        title: { text: "FastQ Screen Results" },
                        xAxis: { categories: fq_screen_dict[plot_id].categories },
                        yAxis: {
                            max: 100,
                            min: 0,
                            title: { text: "Percentage Aligned" }
                        },
                        tooltip: {
                            formatter: function () {
                                return "<b>" + this.series.stackKey.replace("column","") + " - " + this.x + "</b><br/>" +
                                    this.series.name + ": " + this.y + "%<br/>" +
                                    "Total Alignment: " + this.point.stackTotal + "%";
                            },
                        },
                        plotOptions: {
                            column: {
                                pointPadding: 0,
                                groupPadding: 0.02,
                                stacking: "normal"
                            }
                        },
                        series: fq_screen_dict[plot_id].data
                    });
                });
            });
        </script>"""

        return html

    def fqscreen_simple_plot(self):
        """Makes a simple bar plot with summed alignment counts for
        each species, stacked."""

        # First, sum the different types of alignment counts
        data = {}
        org_counts = {}
        for s_name in sorted(self.fq_screen_data):
            data[s_name] = OrderedDict()
            sum_alignments = 0
            for org in self.fq_screen_data[s_name]:
                if org == "total_reads":
                    continue
                try:
                    data[s_name][org] = self.fq_screen_data[s_name][org]["counts"]["one_hit_one_library"]
                except KeyError:
                    log.error(
                        "No counts found for '{}' ('{}'). Could be malformed or very old FastQ Screen results.".format(
                            org, s_name
                        )
                    )
                    continue
                try:
                    data[s_name][org] += self.fq_screen_data[s_name][org]["counts"]["multiple_hits_one_library"]
                except KeyError:
                    pass
                sum_alignments += data[s_name][org]
                org_counts[org] = org_counts.get(org, 0) + data[s_name][org]

            # Calculate hits in multiple genomes
            if "total_reads" in self.fq_screen_data[s_name]:
                data[s_name]["Multiple Genomes"] = self.fq_screen_data[s_name]["total_reads"] - sum_alignments

        # Sort the categories by the total read counts
        cats = OrderedDict()
        for org in sorted(org_counts, key=org_counts.get, reverse=True):
            if org not in cats and org != "No hits":
                cats[org] = {"name": org}

        # Strip empty dicts
        for s_name in list(data.keys()):
            if len(data[s_name]) == 0:
                del data[s_name]

        pconfig = {
            "id": "fastq_screen_plot",
            "title": "FastQ Screen: Mapped Reads",
            "cpswitch_c_active": False,
            "ylab": "Percentages",
        }
        cats["Multiple Genomes"] = {"name": "Multiple Genomes", "color": "#820000"}
        cats["No hits"] = {"name": "No hits", "color": "#cccccc"}

        return bargraph.plot(data, cats, pconfig)

    def fqscreen_bisulfite_plot(self):
        """Make a stacked barplot for the bisulfite data, if we have any"""

        pconfig = {
            "id": "fastq_screen_bisulfite_plot",
            "title": "FastQ Screen: Bisulfite Mapping Strand Orientation",
            "hide_zero_cats": False,
            "ylab": "Percentages",
            "data_labels": [],
        }

        cats = OrderedDict()
        cats["original_top_strand"] = {"name": "Original top strand", "color": "#80cdc1"}
        cats["complementary_to_original_top_strand"] = {
            "name": "Complementary to original top strand",
            "color": "#018571",
        }
        cats["complementary_to_original_bottom_strand"] = {
            "name": "Complementary to original bottom strand",
            "color": "#a6611a",
        }
        cats["original_bottom_strand"] = {"name": "Original bottom strand", "color": "#dfc27d"}

        # Pull out the data that we want
        pdata_unsorted = {}
        org_counts = {}
        max_count = 0
        for s_name, d in self.fq_screen_data.items():
            for org, dd in d.items():
                if org not in pdata_unsorted:
                    pdata_unsorted[org] = {}
                    org_counts[org] = 0
                try:
                    pdata_unsorted[org][s_name] = dd["bisulfite_counts"]
                    total_counts = sum([c for c in dd["bisulfite_counts"].values()])
                    org_counts[org] += total_counts
                    max_count = max(max_count, total_counts)
                except (TypeError, KeyError):
                    pass

        # Consistent y-max across all genomes
        pconfig["ymax"] = max_count

        # Show tabbed bar plots, sorted by total read count
        pdata = []
        pcats = []
        for org in sorted(org_counts, key=org_counts.get, reverse=True):
            if org_counts[org] > 0:
                pconfig["data_labels"].append(org)
                pdata.append(pdata_unsorted[org])
                pcats.append(cats)

        if len(pdata) == 0:
            return None

        self.add_section(
            name="Bisulfite Reads", anchor="fastq_screen_bisulfite", plot=bargraph.plot(pdata, pcats, pconfig)
        )
