""" MultiQC module to parse log output from Xenome Classify """

import logging
from collections import defaultdict
from typing import Dict

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

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self._xenome_general_stats_table(self.summary_data)

        # Alignment Rate Plot
        self._xenome_stats_plot(self.summary_data, self.detail_data)

    def _parse_xenome_logs(self, f):
        """
        Expected two sections: "Statistics" and "Summary", one following the other after a line break
        """
        s_name = f["s_name"]

        lines = iter(f["contents_lines"])
        try:
            detail_cnt_by_class = self._parse_xenome_table(lines, "Statistics")
            summary_cnt_by_class = self._parse_xenome_table(lines, "Summary")
        except (AssertionError, StopIteration) as e:
            log.error(f"Error parsing Xenome log file '{f['fn']}' for sample '{s_name}': {e}")
            return
        if detail_cnt_by_class is None or summary_cnt_by_class is None:
            return

        if s_name in self.detail_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.add_data_source(f, s_name)
        self.detail_data[s_name] = detail_cnt_by_class
        self.summary_data[s_name] = summary_cnt_by_class

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None, s_name)

    @staticmethod
    def _parse_xenome_table(lines, expected_title: str):
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
        return cnt_by_class

    def _xenome_general_stats_table(self, summary_data):
        """
        Add the numbers of reads classified as one of the species into the general stats table.
        """
        headers: Dict[str, Dict] = {}
        table_data = defaultdict(dict)
        for sn, data in summary_data.items():
            for cls, cnt in data.items():
                table_data[sn][f"{cls}_reads"] = cnt
                if cls == "human":
                    headers[f"{cls}_reads"] = {
                        "title": f"Human Reads, {config.read_count_prefix}",
                        "description": f"The number of human reads in the sample ({config.read_count_desc})",
                        "min": 0,
                        "scale": "Blues",
                        "modify": lambda x: x * config.read_count_multiplier,
                        "shared_key": "read_count",
                    }
                else:
                    headers[f"{cls}_reads"] = {
                        "title": f"{cls.capitalize()} Reads",
                        "description": "The number of {cls} reads in the sample",
                        "min": 0,
                        "scale": "Greys" if cls in ["both", "neither", "ambiguous"] else "Reds",
                        "format": "{:,.0f}",
                        "hidden": cls in ["both", "neither", "ambiguous"],
                    }
        self.general_stats_addcols(table_data, headers)

    def _xenome_stats_plot(self, summary_data, detail_data):
        """
        Create two bar plots: based on summary and detail data.
        """
        COLORS = [
            "#377eb8",
            "#e41a1c",
            "#4daf4a",
            "#984ea3",
            "#ff7f00",
            "#ffff33",
            "#a65628",
            "#f781bf",
        ]  # to support multiple species, when running on logs for many unrelated runs
        BOTH = "#984ea3"  # purple
        AMBIGUOUS = "#616161"  # grey
        NEITHER = "#b3b3b3"  # light grey

        def _lighten_color(code, lighten=0.0):
            def _rgb_converter(x):
                return max(0, min(1, 1 + ((x - 1) * lighten)))

            c = spectra.html(code)
            c = spectra.rgb(*[_rgb_converter(v) for v in c.rgb])
            return c.hexcode

        def _get_color(i, lighten=0.0):
            code = COLORS[i % len(COLORS)]
            return _lighten_color(code, lighten)

        grafts = []
        hosts = []
        for sn, cnt_by_cls in summary_data.items():
            s = list(cnt_by_cls)[0]
            if s not in grafts:
                grafts.append(s)
            s = list(cnt_by_cls)[1]
            if s not in hosts:
                hosts.append(s)
        all_species = grafts + [s for s in hosts if s not in grafts]
        colors = {s: _get_color(i, lighten=0.6) for i, s in enumerate(all_species)}
        definitely_colors = {s: _get_color(i, lighten=1.0) for i, s in enumerate(all_species)}
        probably_colors = {s: _get_color(i, lighten=0.3) for i, s in enumerate(all_species)}

        summary_cats = {}
        detail_cats = {}
        for s in all_species:
            summary_cats[s] = {"name": s.capitalize(), "color": colors[s]}
            detail_cats[f"definitely {s}"] = {"name": f"Definitely {s}", "color": definitely_colors[s]}
            detail_cats[f"probably {s}"] = {"name": f"Probably {s}", "color": probably_colors[s]}

        summary_cats["both"] = {"name": "Both", "color": _lighten_color(BOTH, 0.6)}
        summary_cats["ambiguous"] = {"name": "Ambiguous", "color": AMBIGUOUS}
        summary_cats["neither"] = {"name": "Neither", "color": NEITHER}
        detail_cats["both"] = summary_cats["both"]
        detail_cats["probably both"] = {"name": "Probably both", "color": _lighten_color(BOTH, 0.3)}
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
                summary_data,
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
                detail_data,
                detail_cats,
                {
                    "id": "xenome_stats_detail",
                    "title": "Xenome: detailed classification",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of reads",
                },
            ),
        )
