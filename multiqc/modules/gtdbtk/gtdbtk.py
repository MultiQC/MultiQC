import logging

from multiqc.base_module import BaseMultiqcModule
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
            info="Toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes.",
            doi=["10.1093/bioinformatics/btac672"],
        )

        self.gtdbtk_data = {}
        for f in self.find_log_files("gtdbtk", filehandles=True):
            self.parse_file(f)
            self.add_data_source(f)
        self.gtdbtk_data = self.ignore_samples(self.gtdbtk_data)
        if len(self.gtdbtk_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.gtdbtk_data)} reports")
        self.add_software_version(None)  # may need to fix this?

        self.closest_taxa_table()

    def parse_file(self, f):
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
        parsed_data = {}
        for line in f["f"]:
            if line.startswith("user_genome\t"):
                continue
            column_values = line.rstrip().split("\t")
            parsed_data[column_values[0]] = dict(zip(column_names, column_values))
        self.gtdbtk_data.update(parsed_data)

    def closest_taxa_table(self):
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
        for sample_name in self.gtdbtk_data:
            sample = self.gtdbtk_data[sample_name]
            classification_value_columns = classication_method_translate_dict.get(
                sample.get("classification_method"), [None, None]
            )
            table_data[sample.get("user_genome")] = {
                "classification": sample.get("classification", None),
                "classification_method": sample.get("classification_method", None),
                "ANI": sample.get(classification_value_columns[0], None),
                "AF": sample.get(classification_value_columns[1], None),
                "red_value": sample.get("red_value", None),
                "note": sample.get("note"),
                "warnings": sample.get("warnings"),
            }
        headers = {
            "classification": {
                "title": "Classification",
                "description": "GTDB taxonomy string inferred by the GTDB-Tk.",
            },
            "classification_method": {
                "title": "Classification Method",
                "description": "Indicates the rule used to classify the genome.",
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
                "title": "Relative Evolutionary Divergence (RED)",
                "description": "Indicates the relative evolutionary divergence (RED) for a query genome. RED is not calculated when a query genome can be classified based on ANI.",
                "min": 0,
                "max": 1,
            },
            "note": {
                "title": "Notes",
                "description": "Provides additional information regarding the classification of the genome.",
            },
            "warnings": {
                "title": "Warnings",
                "description": "Indicates unusual characteristics of the query genome that may impact the taxonomic assignment.",
            },
        }
        pconfig = TableConfig(
            title="Taxonomy classifications",
            id="gtdbtk-first-table",
        )
        self.add_section(
            name="MAG taxonomy",
            anchor="gtdbtk-taxonomy",
            description="The taxonomy of a MAG as found by GTDB.",
            helptext="GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy GTDB. It is designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate and single-cell genomes. ",
            plot=table.plot(data=table_data, headers=headers, pconfig=pconfig),
        )
