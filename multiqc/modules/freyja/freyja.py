import csv
import logging
import re
from ast import literal_eval
from typing import Dict

from natsort import natsorted

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.plots.table_object import ColumnDict
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)

# Pathogen reported by legacy Freyja versions that do not emit a "pathogen" line.
DEFAULT_PATHOGEN = "SARS-CoV-2"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Freyja",
            anchor="freyja",
            href="https://github.com/andersen-lab/Freyja",
            info="Recovers relative lineage abundances from mixed samples.",
            extra="""
            Freyja is a tool to recover relative lineage abundances from mixed samples from a
            sequencing dataset and uses lineage-determining mutational "barcodes" derived from the UShER global
            phylogenetic tree to solve the constrained (unit sum, non-negative) de-mixing problem.
            """,
            doi="10.1038/s41586-022-05049-6",
            license="BSD 2-Clause License",
            license_url="https://github.com/andersen-lab/Freyja/blob/main/LICENSE",
        )

        # Store the summary data grouped by pathogen:
        #   data_by_pathogen[pathogen][s_name] = {lineage: abundance}
        data_by_pathogen: Dict[str, Dict[str, Dict[str, float]]] = {}

        for f in self.find_log_files("freyja", filehandles=True):
            # The demix report is a TSV with one "key\tvalue" row per field. We need the
            # "summarized" row, and, for newer Freyja versions, the "pathogen" row that names
            # the analysed organism:
            # ...
            # summarized	[('BQ.1*', 0.983), ('Omicron', 0.011), ('key', value)]
            # ...
            # pathogen	SARS-CoV-2
            # ...
            s_name: str = f["s_name"]

            # Collect the value of each field, keyed by the field name in the first column.
            fields: Dict[str, str] = {}
            for row in csv.reader(f["f"], delimiter="\t"):
                if len(row) >= 2:
                    fields[row[0]] = row[1]

            # This will be present because the search pattern requires the "summarized" row to be there:
            summarized = fields["summarized"]

            # Legacy Freyja reports do not include a pathogen field; they only ever covered SARS-CoV-2.
            pathogen = fields.get("pathogen") or DEFAULT_PATHOGEN

            try:
                sample_dict: Dict[str, float] = dict(literal_eval(summarized))
            except (ValueError, SyntaxError):
                log.error(f"Error parsing 'summarized' line for '{s_name}': {summarized}, skipping sample")
                continue
            if not sample_dict:
                log.debug(f"No data in the 'summarized' line for '{s_name}': {summarized}, skipping sample")
                continue

            samples = data_by_pathogen.setdefault(pathogen, {})
            if s_name in samples:
                log.debug(f"Duplicate sample name found for pathogen '{pathogen}'! Overwriting: {s_name}")
            samples[s_name] = sample_dict
            self.add_data_source(f, s_name)

        # Remove filtered samples
        for pathogen in list(data_by_pathogen.keys()):
            data_by_pathogen[pathogen] = self.ignore_samples(data_by_pathogen[pathogen])
            if not data_by_pathogen[pathogen]:
                del data_by_pathogen[pathogen]

        # Let MultiQC know this module found no data
        if not data_by_pathogen:
            raise ModuleNoSamplesFound

        total_samples = sum(len(samples) for samples in data_by_pathogen.values())
        log.info(f"Found {total_samples} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.scale = mqc_colour.mqc_colour_scale("plot_defaults")

        for pathogen in natsorted(data_by_pathogen.keys()):
            # Sort samples to keep a reproducible order - files are discovered non-deterministically
            samples = data_by_pathogen[pathogen]
            data_by_sample = {s_name: samples[s_name] for s_name in natsorted(samples.keys())}

            suffix = f" ({pathogen})"
            # Slugs disambiguate section anchors and column keys per pathogen.
            base_slug = re.sub(r"[^a-zA-Z0-9]+", "-", pathogen).strip("-").lower()
            anchor_slug = f"-{base_slug}"  # e.g. "freyja-summary-sars-cov-2"
            col_slug = f"_{base_slug.replace('-', '_')}"  # e.g. "Top_lineage_freyja_sars_cov_2"

            top_lineages_dict: Dict[str, Dict] = {}
            all_lineages: set = set()
            for s_name, sample_data in data_by_sample.items():
                top_lineage, top_lineage_value = max(sample_data.items(), key=lambda xv: xv[1])
                top_lineages_dict[s_name] = {
                    f"Top_lineage_freyja{col_slug}": top_lineage,
                    f"Top_lineage_freyja{col_slug}_percentage": top_lineage_value,
                }
                all_lineages.add(top_lineage)
            for s_name, sample_data in data_by_sample.items():
                for lineage, val in sorted(sample_data.items(), key=lambda xv: xv[1], reverse=True):
                    if lineage not in all_lineages:
                        all_lineages.add(lineage)

            self.general_stats_cols(top_lineages_dict, all_lineages, suffix, col_slug)
            self.add_freyja_section(all_lineages, data_by_sample, suffix, anchor_slug)

        # Write the data file at the very end, after all sections are added.
        # Disambiguate sample names with the pathogen.
        combined: Dict[str, Dict[str, float]] = {}
        for pathogen, samples in data_by_pathogen.items():
            for s_name, sample_data in samples.items():
                combined[f"{s_name} ({pathogen})"] = sample_data
        self.write_data_file(combined, "multiqc_freyja")

    def general_stats_cols(self, top_lineages_dict, all_lineages, suffix, slug):
        """Add a column displaying the most abundant lineage to the General Statistics table"""
        headers: Dict[str, ColumnDict] = {
            f"Top_lineage_freyja{slug}": {
                "title": f"Top lineage{suffix}",
                "description": f"The most abundant lineage in the sample{suffix}",
                "bgcols": {x: self.scale.get_colour(i) for i, x in enumerate(all_lineages)},
            },
            f"Top_lineage_freyja{slug}_percentage": {
                "title": f"Top lineage %{suffix}",
                "description": f"The percentage of the most abundant lineage in the sample{suffix}",
                "max": 100,
                "min": 0,
                "scale": "Blues",
                "modify": lambda x: x * 100,
                "suffix": "%",
            },
        }
        self.general_stats_addcols(top_lineages_dict, headers)

    def add_freyja_section(self, lineages, data_by_sample, suffix, slug):
        pconfig = {
            "id": f"Freyja_plot{slug}",
            "title": f"Freyja: Top lineages{suffix}",
            "ylab": "relative abundance",
            "y_clipmax": 1,
            "cpswitch": False,
            "cpswitch_c_active": False,
        }
        cats = {x: {"name": x, "color": self.scale.get_colour(i, lighten=1)} for i, x in enumerate(lineages)}

        self.add_section(
            name=f"Freyja Summary{suffix}",
            anchor=f"freyja-summary{slug}",
            description="""
                Relative lineage abundances from mixed samples. Hover over the column headers for descriptions and click _Help_ for more in-depth documentation.
                """,
            helptext="""
                The graph denotes a sum of all lineage abundances in a particular WHO designation, otherwise they are grouped into "Other".
                Lineages abundances are calculated as the number of reads that are assigned to a particular lineage.
                Lineages and their corresponding abundances are summarized by constellation.

                > **Note**: Lineage designation is based on the used WHO nomenclature, which could vary over time.
                """,
            plot=bargraph.plot(data_by_sample, cats, pconfig),
        )
