""" MultiQC module to parse log output from Xenome Classify """

import logging
from collections import defaultdict
from typing import Dict, Union, List

import spectra

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.utils import config

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

        self.show_pct = not getattr(config, "xenome", {}).get("show_read_counts", False)

        # Find and load any Xenome reports
        self.detail_data = dict()
        self.summary_data = dict()
        for f in self.find_log_files("xenome"):
            self._parse_xenome_logs(f)

        # Filter to strip out ignored sample names
        self.detail_data = self.ignore_samples(self.detail_data)
        self.summary_data = self.ignore_samples(self.summary_data)
        if len(self.detail_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.detail_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.detail_data, "multiqc_xenome_detail")
        self.write_data_file(self.summary_data, "multiqc_xenome_summary")

        self.all_species = self._collect_all_species(self.summary_data)

        self._xenome_general_stats_table()
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
            detail_by_class = self._parse_xenome_table(lines, "Statistics")
            summary_by_class = self._parse_xenome_table(lines, "Summary")
        except (AssertionError, StopIteration) as e:
            log.error(f"Error parsing Xenome log file '{f['fn']}' for sample '{s_name}': {e}")
            return
        if detail_by_class is None or summary_by_class is None:
            return

        if s_name in self.detail_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.add_data_source(f, s_name)
        self.detail_data[s_name] = detail_by_class
        self.summary_data[s_name] = summary_by_class

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None, s_name)

    def _parse_xenome_table(self, lines, expected_title: str):
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
        return pct_by_class if self.show_pct else cnt_by_class

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

    def _xenome_general_stats_table(self):
        """
        Add the numbers of reads classified as one of the species into the general stats table.
        """
        headers: Dict[str, Dict] = {}
        table_data = defaultdict(dict)
        for sn, data in self.summary_data.items():
            for cls, cnt in data.items():
                table_data[sn][f"{cls}_reads"] = cnt
                if cls == "human":
                    headers[f"{cls}_reads"] = {
                        "title": f"Human reads",
                        "description": f"{'share' if self.show_pct else 'number'} of human reads in the sample",
                        "min": 0,
                        "suffix": "%" if self.show_pct else "",
                        "scale": self._get_color(cls, return_scale=True),
                        "format": "{:,.1f}" if self.show_pct else "{:,d}",
                    }
                else:
                    headers[f"{cls}_reads"] = {
                        "title": f"{cls.capitalize()} reads",
                        "description": f"{'share' if self.show_pct else 'number'} of {cls} reads in the sample",
                        "min": 0,
                        "suffix": "%" if self.show_pct else "",
                        "scale": self._get_color(cls, return_scale=True),
                        "format": "{:,.1f}" if self.show_pct else "{:,d}",
                        "hidden": cls in ["both", "neither", "ambiguous"],
                    }
        self.general_stats_addcols(table_data, headers)

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
            plot=bargraph.plot(
                self.summary_data,
                summary_cats,
                {
                    "id": "xenome_stats_summary",
                    "title": "Xenome: summary classification",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of reads",
                },
            ),
        )

        self.add_section(
            description="This plot shows the number of reads classified by Xenome, "
            "with a more detail certainty split than the plot above",
            plot=bargraph.plot(
                self.detail_data,
                detail_cats,
                {
                    "id": "xenome_stats_detail",
                    "title": "Xenome: detailed classification",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of reads",
                },
            ),
        )
