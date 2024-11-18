import logging

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses summary.tsv outputs from GTDB-Tk's classify.py and classify_wf.py .

    Only works for >= 2.4.0 because column names changed.
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
            raise ModuleNotFoundError
        log.info(f"Found {len(self.gtdbtk_data)} reports")
        self.add_software_version(None)  # may need to fix this?

        self.closest_taxa_table()

    def parse_file(self, f):
        # loop through every line
        # split on tabs
        # grab the first two columns
        # grab the 13th column to say how it classified
        # classification_method: indicates the rule used to classify the genome.
        # This field will be one of:
        # i) ANI, indicating a species assignement was based solely on the calculated ANI and AF with a reference genome;
        # ii) ANI/Placement, indicating a species assignment was made based on both ANI and the placement of the genome in the reference tree;
        # iii) taxonomic classification fully defined by topology, indicating that the classification could be determine based solely on the genome’s position in the reference tree; or
        # iv) taxonomic novelty determined using RED, indicating that the relative evolutionary divergence (RED) and placement of the genome in the reference tree were used to determine the classification.
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
        self.gtdbtk_data = parsed_data

    def closest_taxa_table(self):
        """Add a table showing the closest taxa for each query genome."""
        # generate a table with :
        # rows: a MAG
        # columns: MAG name, classifcation, classification_method, classification_value, note, warnings
        # classifacation_value is just the following:
        # if classification_method:
        classication_method_translate_dict = {
            "taxonomic classification defined by topology and ANI": ["closest_genome_ani", "closest_genome_af"],
            "ani_screen": ["closest_genome_ani", "closest_genome_af"],
            "taxonomic classification fully defined by topology": [
                "closest_placement_ani",
                "closest_placement_af",
            ],
            "taxonomic novelty determined using RED": ["red_value", ""],
        }
        table_data = {}
        for sample_name in self.gtdbtk_data:
            sample = self.gtdbtk_data[sample_name]
            classification_value_columns = classication_method_translate_dict.get(sample.get("classification_method"), [None, None])
            table_data[sample.get("user_genome")] = {
                "classification": sample.get("classification", None),
                "classification_method": sample.get("classification_method", None),
                "classification_value_1": sample.get(classification_value_columns[0], None),
                "classification_value_1_unit": classification_value_columns[0],
                "classification_value_2": sample.get(classification_value_columns[1], None),
                "classification_value_2_unit": classification_value_columns[1],
                "note": sample.get("note"),
                "warnings": sample.get("warnings"),
            }
        self.add_section(
            name="MAG taxonomy",
            anchor="gtdbtk-taxonomy",
            description="The taxonomy of a MAG as found by GTDB.",
            plot=table.plot(data=table_data),
        )
