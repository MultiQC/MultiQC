""" MultiQC module to parse log output from Xenome Classify """

from collections import defaultdict

import logging
import spectra
from typing import Dict, Union

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Xenome",
            anchor="xenome",
            href="https://github.com/data61/gossamer/blob/master/docs/xenome.md",
            info="is a tool for classifying reads from xenograft sources.",
            doi="doi.org/10.1093/bioinformatics/bts236",
        )

        # Find and load any Xenome reports
        self.detail_percents = dict()
        self.detail_counts = dict()
        self.summary_percents = dict()
        self.summary_counts = dict()
        for f in self.find_log_files("xenome"):
            self._parse_xenome_logs(f)

        # Filter to strip out ignored sample names
        self.detail_percents = self.ignore_samples(self.detail_percents)
        self.detail_counts = self.ignore_samples(self.detail_counts)
        self.summary_percents = self.ignore_samples(self.summary_percents)
        self.summary_counts = self.ignore_samples(self.summary_counts)
        if len(self.detail_percents) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.detail_percents)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.detail_percents, "multiqc_xenome_detail_percent")
        self.write_data_file(self.detail_counts, "multiqc_xenome_detail_counts")
        self.write_data_file(self.summary_percents, "multiqc_xenome_summary_percent")
        self.write_data_file(self.summary_counts, "multiqc_xenome_summary_counts")

        self.all_species = self._collect_all_species(self.summary_counts)
        self._build_table()
        self._xenome_stats_plot()

    @staticmethod
    def _collect_all_species(summary_data):
        grafts = []
        hosts = []
        for sn, cnt_by_cls in summary_data.items():
            s = list(cnt_by_cls)[0]
            if s not in grafts:
                grafts.append(s)
            s = list(cnt_by_cls)[1]
            if s not in hosts:
                hosts.append(s)
        return grafts + [s for s in hosts if s not in grafts]

    def _parse_xenome_logs(self, f):
        """
        Expected two sections: "Statistics" and "Summary", one following the other after a line break
        """
        s_name = f["s_name"]

        lines = iter(f["contents_lines"])
        try:
            detail_percents, detail_counts = self._parse_xenome_section(lines, "Statistics")
            summary_percents, summary_counts = self._parse_xenome_section(lines, "Summary")
        except (AssertionError, StopIteration) as e:
            log.error(f"Error parsing Xenome log file '{f['fn']}' for sample '{s_name}': {e}")
            return

        if s_name in self.detail_counts:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.add_data_source(f, s_name)
        self.detail_counts[s_name] = detail_counts
        self.detail_percents[s_name] = detail_percents
        self.summary_counts[s_name] = summary_counts
        self.summary_percents[s_name] = summary_percents

    @staticmethod
    def _parse_xenome_section(lines, expected_title: str):
        """
        Parse one section, e.g.:

        Statistics
        B	G	H	M	count	percent	class
        0	0	0	0	73938	0.231509	"neither"
        ...

        OR

        Summary
        count	percent	class
        30220658	94.6246	human
        ...
        """
        title = next(lines).strip()
        assert title == expected_title, f"Expected '{expected_title}', got '{title}'"
        header = next(lines).strip().split("\t")
        expected_fields = ["count", "percent", "class"]
        assert all(
            field in header for field in expected_fields
        ), f"Expected the header to contain '{expected_fields}', got {header}"

        cnt_by_class = defaultdict(int)
        pct_by_class = defaultdict(float)
        for line in lines:
            if line.strip() == "":
                break
            fields = line.strip().split("\t")
            assert len(fields) == len(header), f"Expected {len(header)} fields, got {len(fields)}"
            data = dict(zip(header, fields))
            cls = data["class"].strip('"')
            cnt_by_class[cls] += int(data["count"])
            pct_by_class[cls] += float(data["percent"])
        return pct_by_class, cnt_by_class

    # to support multiple species, when running on logs for many unrelated runs
    COLORS = [
        ("#377eb8", "Blues"),  # blue
        ("#e41a1c", "Reds"),  # red
        ("#4daf4a", "Greens"),  # green
        ("#ff7f00", "Oranges"),  # orange
        ("#ffff33", "YlOrRd"),  # yellow
        ("#a65628", "YlOrBr"),  # brown
        ("#f781bf", "PuRd"),  # pink
    ]

    BOTH = ("#984ea3", "Purples")  # purple
    AMBIGUOUS = ("#616161", "Greys")  # grey
    NEITHER = ("#b3b3b3", "Greys")  # light grey

    @staticmethod
    def _lighten_color(code: str, lighten=1.0):
        def _rgb_converter(x):
            return max(0, min(1, 1 + ((x - 1) * lighten)))

        c = spectra.html(code)
        c = spectra.rgb(*[_rgb_converter(v) for v in c.rgb])
        return c.hexcode

    def _get_color(self, cls: Union[str, int], lighten=1.0, return_scale=False) -> str:
        if not isinstance(cls, int) and cls in self.all_species:
            cls = self.all_species.index(cls)
        if isinstance(cls, int):
            code, scale = MultiqcModule.COLORS[cls % len(MultiqcModule.COLORS)]
        elif cls == "both":
            code, scale = MultiqcModule.BOTH
        elif cls == "ambiguous":
            code, scale = MultiqcModule.AMBIGUOUS
        else:
            code, scale = MultiqcModule.NEITHER

        if return_scale:
            return scale
        else:
            return MultiqcModule._lighten_color(code, lighten)

    def _build_table(self):
        """
        Prepare headers and data for a table. Add a section with a table,
        and add a few columns into the general stats.
        """
        headers: Dict[str, Dict] = {}
        table_data = defaultdict(dict)

        for sn, data in self.summary_percents.items():
            for cls, val in data.items():
                table_data[sn][f"{cls}_reads_pct"] = val
                if cls == "human":
                    headers[f"{cls}_reads_pct"] = {
                        "rid": f"{self.anchor}_{cls}_reads_pct",  # to make the ID unique from xengsort
                        "title": "Human reads",
                        "description": "share of human reads in the sample",
                        "min": 0,
                        "suffix": "%",
                        "scale": self._get_color(cls, return_scale=True),
                        "format": "{:,.1f}",
                    }
                else:
                    headers[f"{cls}_reads_pct"] = {
                        "rid": f"{self.anchor}_{cls}_reads_pct",  # to make the ID unique from xengsort
                        "title": f"{cls.capitalize()} reads",
                        "description": f"share of {cls} reads in the sample",
                        "min": 0,
                        "suffix": "%",
                        "scale": self._get_color(cls, return_scale=True),
                        "format": "{:,.1f}",
                        "hidden": cls in ["both", "neither", "ambiguous"],
                    }
        self.general_stats_addcols(table_data, headers)
        detail_headers = headers.copy()
        for metric in headers:
            detail_headers[metric]["hidden"] = False

        for sn, data in self.summary_counts.items():
            for cls, val in data.items():
                table_data[sn][f"{cls}_reads_cnt"] = val
                if cls == "human":
                    detail_headers[f"{cls}_reads_cnt"] = {
                        "title": "Human reads",
                        "description": "number of human reads in the sample",
                        "min": 0,
                        "scale": self._get_color(cls, return_scale=True),
                        "format": "{:,d}",
                        "hidden": True,
                    }
                else:
                    detail_headers[f"{cls}_reads_cnt"] = {
                        "title": f"{cls.capitalize()} reads",
                        "description": f"number of {cls} reads in the sample",
                        "min": 0,
                        "scale": self._get_color(cls, return_scale=True),
                        "format": "{:,d}",
                        "hidden": True,
                    }

        self.add_section(
            name="Summary table",
            anchor="xenome-summary-table",
            plot=table.plot(table_data, detail_headers, pconfig={"id": "xenome-table"}),
        )

    def _xenome_stats_plot(self):
        """
        Create two bar plots: based on summary and detail data.
        """
        colors = {s: self._get_color(i, lighten=0.6) for i, s in enumerate(self.all_species)}
        definitely_colors = {s: self._get_color(i, lighten=1.0) for i, s in enumerate(self.all_species)}
        probably_colors = {s: self._get_color(i, lighten=0.3) for i, s in enumerate(self.all_species)}

        summary_cats = {}
        detail_cats = {}
        for s in self.all_species:
            summary_cats[s] = {"name": s.capitalize(), "color": colors[s]}
            detail_cats[f"definitely {s}"] = {"name": f"Definitely {s}", "color": definitely_colors[s]}
            detail_cats[f"probably {s}"] = {"name": f"Probably {s}", "color": probably_colors[s]}

        summary_cats["both"] = {"name": "Both", "color": self._get_color("both", lighten=0.6)}
        summary_cats["ambiguous"] = {"name": "Ambiguous", "color": self._get_color("ambiguous")}
        summary_cats["neither"] = {"name": "Neither", "color": self._get_color("neither")}
        detail_cats["both"] = summary_cats["both"]
        detail_cats["probably both"] = {
            "name": "Probably both",
            "color": self._get_color("both", 0.3),
        }
        detail_cats["ambiguous"] = summary_cats["ambiguous"]
        detail_cats["neither"] = summary_cats["neither"]

        self.add_section(
            description="This plot shows the number of reads classified by Xenome",
            helptext="""
            There are 5 possible categories:  
            * Reads found in graft species, e.g. **Human**
            * Reads found in host species, e.g. **Mouse**  
            * **Both**: read was found in either of the species
            * **Neither**: Read was found in neither of the species
            * **Ambiguous**: Read origin could not be adequately determined.  
            """,
            name="Summary classification",
            anchor="xenome_summary_bar_plot_section",
            plot=bargraph.plot(
                self.summary_counts,
                summary_cats,
                {
                    "id": "xenome_summary_bar_plot",
                    "title": "Xenome: summary classification",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of reads",
                    "cpswitch_c_active": False,
                },
            ),
        )

        self.add_section(
            description="This plot shows the number of reads classified by Xenome, "
            "with a more detail certainty split than the plot above",
            name="Detailed classification",
            anchor="xenome_detail_bar_plot_section",
            plot=bargraph.plot(
                self.detail_counts,
                detail_cats,
                {
                    "id": "xenome_detail_bar_plot",
                    "title": "Xenome: detailed classification",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of reads",
                    "cpswitch_c_active": False,
                },
            ),
        )
