import logging
from collections import defaultdict
from statistics import median

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module assumes that the Percolator output file is named `*percolator_feature_weights.tsv`.
    Make sure to run it using:

    ```
    percolator ... > samples.percolator_feature_weights.tsv
    ```

    The module accepts one configuration option:
     - `group_to_feature`: A dictionary mapping group names to feature names (empty per default), e.g. in `multiqc_config.yaml`:

    ```yaml
    percolator:
      group_to_feature:
        psm_file_combined: [MS:1002255, MS:1002252]
        ms2pip: [ionb_min_abs_diff, iony_min_abs_diff]
        deeplc: [rt_diff]
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Percolator",
            anchor="percolator",
            href="https://github.com/percolator/percolator",
            info="Semi-supervised learning for peptide identification from shotgun proteomics datasets.",
            doi="10.1007/s13361-016-1460-7",
        )

        # Parse logs
        self.data_by_sample = dict()
        for f in self.find_log_files("percolator", filehandles=True):
            self.parse_percolator(f)

        # Filter to strip out ignored sample names
        self.data_by_sample = self.ignore_samples(self.data_by_sample)

        if len(self.data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.data_by_sample)} logs")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with an actual version if it is available
        self.add_software_version(None)

        # Summarize the data from all samples by computing the median across all samples for each feature
        values_by_feature = defaultdict(list)
        for sample_dict in self.data_by_sample.values():
            for feature, value in sample_dict.items():
                values_by_feature[feature].extend(value)
        # Since our bars are features and not samples, computing the median of each feature across samples
        median_weights_by_feature = dict()
        for feature, values in values_by_feature.items():
            median_value = median(values)
            median_weights_by_feature[feature] = median_value
        # Sort the features by the median weight
        median_weights_by_feature = dict(sorted(median_weights_by_feature.items(), key=lambda x: -x[1], reverse=True))

        # Color feature groups if provided in the config
        group_to_features = getattr(config, "percolator", {}).get("group_to_feature", {})
        feature_to_group = {feature: group for group, features in group_to_features.items() for feature in features}
        # If the feature is not assigned to a group, assign it to 'other'
        barplot_data = {
            feature: {f"{feature_to_group.get(feature, 'other') if feature_to_group else 'weight'}": median_value}
            for feature, median_value in median_weights_by_feature.items()
        }

        # Plot the bar plot of the median weight for each feature
        self.add_section(
            anchor="percolator_median_feature_weights",
            helptext="""The bar plot illustrates the median weights for each feature, 
        which have been assigned to one of three categories: psm_file_combined, ms2pip and deeplc.""",
            description="""For each feature, the associated weight is the median value calculated over all input samples.""",
            plot=bargraph.plot(
                barplot_data,
                pconfig={
                    "id": "percolator_plot",
                    "title": "Percolator: Median Feature Weights",
                    "sort_samples": False,
                    "hide_zero_cats": False,
                },
            ),
        )

        self.write_data_file(self.data_by_sample, "multiqc_percolator_barplot")

    def parse_percolator(self, file):
        """Extract the normalized weights for each feature from the percolator logs."""
        s_name = file["s_name"]
        lines = [line for line in file["f"] if not line.startswith("#")]
        header = lines[0].strip().split("\t")
        if s_name in self.data_by_sample:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.data_by_sample[s_name] = {feature: list() for feature in header}
        # Select every third line starting from the first data line after the header (normalized values)
        for i in range(1, len(lines), 3):
            line = lines[i].strip().split("\t")
            # Take the absolute value of each element
            line = [float(x) for x in line]
            assert len(line) == len(header)
            # For each sample, collect all values for each feature in a list
            for elem in range(len(line)):
                self.data_by_sample[s_name][header[elem]].append(line[elem])
        self.add_data_source(file, s_name)
