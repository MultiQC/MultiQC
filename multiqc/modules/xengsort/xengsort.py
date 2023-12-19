""" MultiQC module to parse log output from xengsort classify """

from collections import defaultdict

import logging
from typing import Dict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="xengsort",
            anchor="xengsort",
            href="https://gitlab.com/genomeinformatics/xengsort",
            info="is a fast xenograft read sorter based on space-efficient k-mer hashing",
            doi="doi.org/10.4230/LIPIcs.WABI.2020.4",
        )

        # Find and load any Xenome reports
        self.percents = dict()
        self.counts = dict()
        for f in self.find_log_files("xengsort"):
            self._parse_log(f)

        # Filter to strip out ignored sample names
        self.percents = self.ignore_samples(self.percents)
        self.counts = self.ignore_samples(self.counts)
        if len(self.percents) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.counts)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.percents, f"multiqc_{self.anchor}_percents")
        self.write_data_file(self.counts, f"multiqc_{self.anchor}_counts")

        self._build_table()
        self._build_plot()

    def _parse_log(self, f):
        lines = iter(f["contents_lines"])
        for line in iter(lines):
            if "\t" in line:
                fields = line.strip().split("\t")
                if set(fields) == {"prefix", "host", "graft", "ambiguous", "both", "neither"}:
                    values = next(lines).strip().split("\t")
                    data = dict(zip(fields, values))
                    s_name = data.pop("prefix")
                    f["s_name"] = s_name
                    data = {k: int(v) for k, v in data.items()}
                    percents = {k: v / sum(data.values()) * 100 for k, v in data.items()}

                    if s_name in self.counts:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    self.counts[s_name] = data
                    self.percents[s_name] = percents
                    break

    def _build_table(self):
        """
        Prepare headers and data for a table. Add a section with a table,
        and add a few columns into the general stats.
        """
        headers: Dict[str, Dict] = {}
        table_data = defaultdict(dict)

        scale_by_cls = {
            "graft": "Blues",
            "host": "Reds",
            "both": "Purples",
            "ambiguous": "Greys",
            "neither": "Greys",
        }
        for sn, data in self.percents.items():
            for cls in ["graft", "host", "ambiguous", "both", "neither"]:
                table_data[sn][f"{cls}_reads_pct"] = data.get(cls)
                headers[f"{cls}_reads_pct"] = {
                    "rid": f"{self.anchor}_{cls}_reads_pct",  # to make the ID unique from xenome
                    "title": f"{cls.capitalize()} reads",
                    "description": f"share of {cls} reads in the sample",
                    "min": 0,
                    "suffix": "%",
                    "scale": scale_by_cls[cls],
                    "format": "{:,.2f}",
                    "hidden": cls in ["both", "neither", "ambiguous"],
                }
        self.general_stats_addcols(table_data, headers)
        detail_headers = headers.copy()
        for metric in headers:
            detail_headers[metric]["hidden"] = False

        for sn, data in self.counts.items():
            for cls in ["graft", "host", "ambiguous", "both", "neither"]:
                table_data[sn][f"{cls}_reads_cnt"] = data.get(cls)
                detail_headers[f"{cls}_reads_cnt"] = {
                    "rid": f"{self.anchor}_{cls}_reads_cnt",  # to make the ID unique from xenome
                    "title": f"{cls.capitalize()} reads",
                    "description": f"number of {cls} reads in the sample",
                    "min": 0,
                    "scale": scale_by_cls.get(cls),
                    "format": "{:,d}",
                    "hidden": True,
                }

        self.add_section(
            name="Summary table",
            anchor=f"{self.anchor}-summary-table-section",
            plot=table.plot(table_data, detail_headers, {"id": f"{self.anchor}-summary-table"}),
        )

    def _build_plot(self):
        """
        Create two bar plots: based on summary and detail data.
        """
        cats = {
            "graft": {"name": "Graft", "color": "#377eb8"},  # blue
            "host": {"name": "Host", "color": "#e41a1c"},  # red
            "both": {"name": "Both", "color": "#984ea3"},  # purple
            "ambiguous": {"name": "Ambiguous", "color": "#616161"},  # grey
            "neither": {"name": "Neither", "color": "b3b3b3"},  # light grey
        }
        self.add_section(
            description=f"This plot shows the number of reads classified by {self.name}",
            helptext="""
            There are 5 possible categories:  
            * **Graft**: reads found in graft species, e.g. human
            * **Host**: reads found in host species, e.g. mouse
            * **Both**: reads found in either of the species
            * **Neither**: reads was found in neither of the species
            * **Ambiguous**: reads origin could not be adequately determined.  
            """,
            name="Summary classification",
            anchor=f"{self.anchor}_summary_bar_plot_section",
            plot=bargraph.plot(
                self.counts,
                cats,
                {
                    "id": f"{self.anchor}_summary_bar_plot",
                    "title": f"{self.name}: summary classification",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of reads",
                    "cpswitch_c_active": False,
                },
            ),
        )
