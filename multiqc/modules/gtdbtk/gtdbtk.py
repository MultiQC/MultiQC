import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.plots.table_object import TableConfig

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses `summary.tsv` outputs from GTDB-Tk's `classify.py` and `classify_wf.py`.

    The module only works for version >= 2.4.0 because column names changed.

    `classify.py` and `classify_wf.py` are used to determine the taxonomic classification of input genomes.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="GTDB-Tk",
            anchor="gtdbtk",
            href="https://ecogenomics.github.io/GTDBTk/index.html",
            info="Assigns objective taxonomic classifications to bacterial and archaeal genomes.",
            doi=["10.1093/bioinformatics/btac672"],
        )

        data_by_sample = {}
        for f in self.find_log_files("gtdbtk"):
            self.parse_file(f, data_by_sample)
            self.add_data_source(f)

        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version()

        self.write_data_file(data_by_sample, "multiqc_gtdbtk")

        self.closest_taxa_table(data_by_sample)

        # Add important columns to the general stats table
        headers = {
            "classification": {
                "title": "Classification",
                "description": "GTDB taxonomy string inferred by the GTDB-Tk.",
            },
            "ANI": {
                "title": "ANI to closest genome",
                "description": "Depending on the classification method, either the 'closest_genome_ani' or 'closest_placement_ani'.",
                "min": 0,
                "max": 100,
                "hidden": True,
            },
            "AF": {
                "title": "AF to closest genome",
                "description": "Depending on the classification method, either the 'closest_genome_af' or 'closest_placement_af'.",
                "min": 0,
                "max": 1,
                "scale": "Purples",
                "hidden": True,
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

    def parse_file(self, f, data_by_sample):
        """Parse the summary.tsv outputs."""
        column_names = (
            "user_genome",
            "classification",
            "closest_genome_reference",
            "closest_genome_reference_radius",
            "closest_genome_taxonomy",
            "closest_genome_ani",
            "closest_genome_af",
            "closest_placement_reference",
            "closest_placement_radius",
            "closest_placement_taxonomy",
            "closest_placement_ani",
            "closest_placement_af",
            "pplacer_taxonomy",
            "classification_method",
            "note",
            "other_related_references",
            "msa_percent",
            "translation_table",
            "red_value",
            "warnings",
        )
        lines = f["f"].splitlines()
        if len(lines) <= 1:
            log.warning(f"Skipping file {f['fn']} because it has no data")
            return
        for line in lines[1:]:
            row = line.rstrip("\n").split("\t")
            if len(row) != len(column_names):
                log.warning(f"Skipping line {line} because it has {len(row)} columns instead of {len(column_names)}")
                continue
            sname = row[0]
            if sname in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {sname}")
            data_by_sample[sname] = {k: v for k, v in zip(column_names[1:], row[1:]) if v}

    def closest_taxa_table(self, data_by_sample):
        """Add a table showing the closest taxa for each query genome."""
        classication_method_translate_dict = {
            "taxonomic classification defined by topology and ANI": ["closest_genome_ani", "closest_genome_af"],
            "ani_screen": ["closest_genome_ani", "closest_genome_af"],
            "taxonomic classification fully defined by topology": [
                "closest_placement_ani",
                "closest_placement_af",
            ],
            "taxonomic novelty determined using RED": [None, None],
        }
        table_data = {}
        for sname, data in data_by_sample.items():
            classification_value_columns = classication_method_translate_dict.get(
                data.get("classification_method"), [None, None]
            )
            classification = data.get("classification", "").split(";")
            last_classification = classification[-1]
            table_data[sname] = {
                "classification": last_classification,
                "full_classification": "; ".join(classification),
                "classification_method": data.get("classification_method", None),
                "ANI": data.get(classification_value_columns[0], None),
                "AF": data.get(classification_value_columns[1], None),
                "red_value": data.get("red_value", None),
                "note": data.get("note"),
                "warnings": data.get("warnings"),
            }
        headers = {
            "classification": {
                "title": "Classification",
                "description": "GTDB taxonomy string inferred by the GTDB-Tk.",
            },
            "full_classification": {
                "title": "Full classification",
                "description": "Full GTDB taxonomy string inferred by the GTDB-Tk.",
                "hidden": True,
            },
            "classification_method": {
                "title": "Classification method",
                "description": "Indicates the rule used to classify the genome.",
                "hidden": True,
            },
            "ANI": {
                "title": "ANI to closest genome",
                "description": "Depending on the classification method, either the 'closest_genome_ani' or 'closest_placement_ani'.",
                "min": 0,
                "max": 100,
            },
            "AF": {
                "title": "AF to closest genome",
                "description": "Depending on the classification method, either the 'closest_genome_af' or 'closest_placement_af'.",
                "min": 0,
                "max": 1,
                "scale": "Purples",
            },
            "red_value": {
                "title": "RED",
                "description": "Indicates the relative evolutionary divergence (RED) for a query genome. RED is not calculated when a query genome can be classified based on ANI.",
                "min": 0,
                "max": 1,
            },
            "warnings": {
                "title": "Warnings",
                "description": "Indicates unusual characteristics of the query genome that may impact the taxonomic assignment.",
            },
            "note": {
                "title": "Notes",
                "description": "Provides additional information regarding the classification of the genome.",
            },
        }
        pconfig = TableConfig(
            title="Taxonomy classifications",
            id="gtdbtk-first-table",
            col1_header="User genome",
        )
        self.add_section(
            name="MAG taxonomy",
            anchor="gtdbtk-taxonomy",
            description="The taxonomy of a MAG as found by GTDB.",
            helptext="GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy GTDB. It is designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate and single-cell genomes. ",
            plot=table.plot(data=table_data, headers=headers, pconfig=pconfig),
        )
