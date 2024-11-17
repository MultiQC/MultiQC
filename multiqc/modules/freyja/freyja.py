import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Freyja",
            anchor="freyja",
            href="https://github.com/andersen-lab/Freyja",
            info="Recovers relative lineage abundances from mixed SARS-CoV-2 samples.",
            extra="""
            Freyja is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a 
            sequencing dataset and uses lineage-determining mutational "barcodes" derived from the UShER global 
            phylogenetic tree to solve the constrained (unit sum, non-negative) de-mixing problem.
            """,
            doi="10.1038/s41586-022-05049-6",
        )

        # To store the summary data
        data_by_sample: Dict[str, Dict] = dict()

        for f in self.find_log_files("freyja", filehandles=True):
            # Freyja has multiple summary files, but we only need to parse the one from the demix command.
            # More specifically, we only need the line that starts with "summarized".
            # ...
            # summarized	[('BQ.1*', 0.983), ('Omicron', 0.011), ('key', value)]
            # ...

            s_name: str = f["s_name"]

            # This will not raise because the search pattern requires it to be there:
            summarized_line = next(line for line in f["f"] if line.startswith("summarized\t"))
            dict_str = summarized_line.split("\t")[1].strip()
            try:
                sample_dict: Dict[str, float] = dict(eval(dict_str))
            except ValueError:
                log.error(f"Error parsing 'summarized' line for '{s_name}': {dict_str}, skipping sample")
                continue
            if not sample_dict:
                log.debug(f"No data in the 'summarized' line for '{s_name}': {dict_str}, skipping sample")
                continue

            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            data_by_sample[s_name] = sample_dict
            self.add_data_source(f, s_name)

        # Remove filtered samples
        data_by_sample = self.ignore_samples(data_by_sample)

        # Let MultiQC know this module found no data
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.write_data_file(data_by_sample, "multiqc_freyja")

        # sort data_by_sample to keep the reproducible order - as the files are discovered in a non-deterministic order
        data_by_sample = dict(sorted(data_by_sample.items()))

        top_lineages_dict = {}
        all_lineages = set()
        for s_name, sample_data in data_by_sample.items():
            top_lineage, top_lineage_value = max(sample_data.items(), key=lambda xv: xv[1])
            top_lineages_dict[s_name] = {
                "Top_lineage_freyja": top_lineage,
                "Top_lineage_freyja_percentage": top_lineage_value,
            }
            all_lineages.add(top_lineage)
        for s_name, sample_data in data_by_sample.items():
            for lineage, val in sorted(sample_data.items(), key=lambda xv: xv[1], reverse=True):
                if lineage not in all_lineages:
                    all_lineages.add(lineage)

        self.scale = mqc_colour.mqc_colour_scale("plot_defaults")
        self.general_stats_cols(top_lineages_dict, all_lineages)
        self.add_freyja_section(all_lineages, data_by_sample)

    def general_stats_cols(self, top_lineages_dict, all_lineages):
        """Add a single column displaying the most abundant lineage to the General Statistics table"""
        headers = {
            "Top_lineage_freyja": {
                "title": "Top lineage",
                "description": "The most abundant lineage in the sample",
                "bgcols": {x: self.scale.get_colour(i) for i, x in enumerate(all_lineages)},
            },
            "Top_lineage_freyja_percentage": {
                "title": "Top lineage %",
                "description": "The percentage of the most abundant lineage in the sample",
                "max": 100,
                "min": 0,
                "scale": "Blues",
                "modify": lambda x: x * 100,
                "suffix": "%",
            },
        }

        self.general_stats_addcols(top_lineages_dict, headers)

    def add_freyja_section(self, lineages, data_by_sample):
        pconfig = {
            "id": "Freyja_plot",
            "title": "Freyja: Top lineages",
            "ylab": "relative abundance",
            "y_clipmax": 1,
            "cpswitch": False,
            "cpswitch_c_active": False,
        }
        cats = {x: {"name": x, "color": self.scale.get_colour(i, lighten=1)} for i, x in enumerate(lineages)}

        self.add_section(
            name="Freyja Summary",
            anchor="freyja-summary",
            description="""
                Relative lineage abundances from mixed SARS-CoV-2 samples. Hover over the column headers for descriptions and click _Help_ for more in-depth documentation. 
                """,
            helptext="""
                The graph denotes a sum of all lineage abundances in a particular WHO designation , otherwise they are grouped into "Other".
                Lineages abundances are calculated as the number of reads that are assigned to a particular lineage. 
                Lineages and their corresponding abundances are summarized by constellation. 

                > **Note**: Lineage designation is based on the used WHO nomenclature, which could vary over time. 
                """,
            plot=bargraph.plot(data_by_sample, cats, pconfig),
        )
