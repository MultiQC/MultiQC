import csv
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap
from multiqc.types import SectionAlert

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module produces summary statistics from AMRFinder and AMRFinderPlus reports.

    The sample name is parsed from the cleaned file name.

    The software version information for AMRFinder can be added to the MultiQC report using MultiQC configuration, or by creating a stand-alone YAML file (https://docs.seqera.io/multiqc/reports/customisation#listing-software-versions).
    """

    def __init__(self):
        super().__init__(
            name="AMRFinder",
            anchor="amrfinder",
            href="https://github.com/ncbi/amr",
            info="Finds acquired antimicrobial resistance genes and point mutations in protein and/or assembled nucleotide sequences",
            doi="10.1038/s41598-021-91456-0",
        )

        self.plot_data = []
        self.general_stats_data = {}
        self.samples_w_no_elements = []

        # Parse AMRFinder log files
        for f in self.find_log_files("amrfinder", filehandles=True):
            self.parse_amrfinder_log(f)

        # Raise error if no samples are found
        if len(self.general_stats_data) == 0:
            raise ModuleNoSamplesFound

        # Adding to general data at top of report
        self.amrfinder_general_stats_table()

        # Adding section for bar graph of element types
        self.add_section(
            name="Element Type Bar Graph",
            anchor="element_type_bar_graph",
            description="A bar graph showing the count of each element type or subtype found by AMRFinder for each sample",
            plot=self.amrfinder_bar_graph(),
            alerts=SectionAlert(
                message=f"**{len(self.samples_w_no_elements)} samples** with no elements hidden from this plot.",
                level="warning",
                affected_samples=self.samples_w_no_elements,
            ),
        )

        # Adding section for heatmap
        self.add_section(
            name="Element Coverage Heatmap",
            anchor="element_coverage_heatmap",
            description="A heatmap showing the coverage for all elements found by AMRFinder",
            plot=self.amrfinder_heatmap(),
            alerts=SectionAlert(
                message=f"**{len(self.samples_w_no_elements)} samples** with no elements hidden from this plot.",
                level="warning",
                affected_samples=self.samples_w_no_elements,
            ),
        )

    def parse_amrfinder_log(self, f):
        sample_name = f["s_name"]
        row_count = 0
        coverage_sum = 0.0
        identity_sum = 0.0
        reader = csv.DictReader(f["f"], delimiter="\t")

        # if the sample is not in ignore samples, then add it to the data frames
        if not self.is_ignore_sample(sample_name):
            # Warning for overwriting duplicate samples
            if sample_name in self.general_stats_data:
                log.debug(f"Duplicate sample name found! Overwriting: {sample_name}")

            # Adding to data sources file
            self.add_data_source(f)

            # Iterating through each row in the file
            for row in reader:
                row["Sample"] = sample_name  # Add sample ID

                self.plot_data.append(row)
                row_count += 1
                coverage_sum += float(row["% Coverage of reference"])
                identity_sum += float(row["% Identity to reference"])
            # adding a row to general stats data for this sample that has the count of elements, average coverage, and average identity
            self.general_stats_data[sample_name] = {
                "elements": row_count,
                "avg_coverage": coverage_sum / row_count if row_count > 0 else 0,
                "avg_identity": identity_sum / row_count if row_count > 0 else 0,
            }

            if row_count == 0:
                self.samples_w_no_elements.append(sample_name)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, sample_name)

    def amrfinder_general_stats_table(self):
        headers = {
            "elements": {
                "title": "# AMRFinder Elements Identified",
                "description": "Count of elements found in each sample",
                "format": "{:,}",
                "scale": "Blues",
            },
            "avg_coverage": {
                "title": "AMRFinder Elements Avg Cov",
                "description": "Mean coverage across elements found in each sample",
                "format": "{:.2f}%",
                "scale": "Greens",
                "hidden": True,
            },
            "avg_identity": {
                "title": "AMRFinder Elements Avg Identity",
                "description": "Mean identity across elements found in each sample",
                "format": "{:.2f}%",
                "scale": "Greens",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.general_stats_data, headers)

    def amrfinder_bar_graph(self):
        # creating a dictionary of dictionaries that has a key for each sample, and each sample has a dictionary with the count of each element type: {sample: {type: count}}
        element_type_counts = {}
        for row in self.plot_data:
            sample = row["Sample"]
            if sample not in element_type_counts:
                element_type_counts[sample] = {}
            element_type = row["Type"]
            if element_type not in element_type_counts[sample]:
                element_type_counts[sample][element_type] = 1
            else:
                element_type_counts[sample][element_type] += 1

        # similar to above but for element subtypes {sample: {subtype: count}}
        element_subtype_counts = {}
        for row in self.plot_data:
            sample = row["Sample"]
            if sample not in element_subtype_counts:
                element_subtype_counts[sample] = {}
            element_subtype = row["Subtype"]
            if element_subtype not in element_subtype_counts[sample]:
                element_subtype_counts[sample][element_subtype] = 1
            else:
                element_subtype_counts[sample][element_subtype] += 1

        # Config for the plot
        pconfig = {
            "data_labels": ["Type", "Subtype"],
            "id": "amrfinder_bar_graph",
            "title": "AMRFinder: Element Type",
            "xlab": "Count",
            "ylab": "Sample",
            "tt_decimals": 0,
        }

        return bargraph.plot([element_type_counts, element_subtype_counts], pconfig=pconfig)

    def amrfinder_heatmap(self):
        # Pivot plot_data to: {element: {sample: coverage}}
        heatmap_data = {}
        for row in self.plot_data:
            element = row["Element symbol"]
            sample = row["Sample"]
            value = float(row["% Coverage of reference"])
            if element not in heatmap_data:
                heatmap_data[element] = {}
            heatmap_data[element][sample] = value

        # sorting heatmap data so that the elements are in alphabetical order
        heatmap_data = dict(sorted(heatmap_data.items()))

        # Config for the plot
        pconfig = {
            "id": "amrfinder_heatmap",
            "title": "AMRFinder: Element Coverage",
            "xlab": "Sample",
            "ylab": "Element",
            "zlab": "% Coverage",
        }

        return heatmap.plot(data=heatmap_data, pconfig=pconfig)
