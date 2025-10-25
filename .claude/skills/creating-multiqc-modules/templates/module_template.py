"""
MultiQC Module Template

Replace {TOOLNAME}, {Display Name}, etc. with actual values.
"""

import logging
from typing import Dict, Union

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    MultiQC module to parse output from **{Display Name}**.

    {Display Name} is a bioinformatics tool that {brief description}.

    The module can summarize data from the following {Display Name} outputs:
    - `{command}` - {command description}

    ### Example Output

    The module generates:
    - **General Statistics**: {list key metrics}
    - **{Plot Name}**: {plot description}

    For more information, see the
    [{Display Name} documentation]({homepage_url}).
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="{Display Name}",
            anchor="{toolname}",
            href="{homepage_url}",
            info="{One-line description starting with capital letter}",
            # doi="{doi}" if available
        )

        # Dictionary to store parsed data by sample
        data_by_sample: Dict[str, Dict[str, Union[int, float, str]]] = {}

        # Find and parse log files
        for f in self.find_log_files("{toolname}"):
            parsed_data = self.parse_{toolname}_log(f["f"])

            if parsed_data:
                s_name = f["s_name"]

                # Handle duplicate sample names
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

                data_by_sample[s_name] = parsed_data
                self.add_data_source(f)

        # Filter out ignored samples
        data_by_sample = self.ignore_samples(data_by_sample)

        # Raise error if no samples found
        if not data_by_sample:
            raise ModuleNoSamplesFound

        # Add software version (required even if None)
        version = self.extract_version(data_by_sample)
        self.add_software_version(version)

        log.info(f"Found {len(data_by_sample)} reports")

        # Add data to general statistics table
        self.add_general_stats(data_by_sample)

        # Add plots and sections
        self.add_{toolname}_plots(data_by_sample)

        # Write data file (MUST BE LAST!)
        self.write_data_file(data_by_sample, "multiqc_{toolname}")

    def parse_{toolname}_log(self, f: str) -> Dict[str, Union[int, float, str]]:
        """
        Parse {toolname} log file and extract metrics.

        Args:
            f: File contents as string

        Returns:
            Dictionary of metrics, or None if parsing failed
        """
        data = {}

        # TODO: Implement parsing logic based on file format
        # Examples:
        #
        # For TSV files:
        # for line in f.splitlines():
        #     parts = line.strip().split("\t")
        #     if len(parts) >= 2:
        #         data[parts[0]] = parts[1]
        #
        # For JSON files:
        # import json
        # data = json.loads(f)
        #
        # For key-value pairs:
        # import re
        # for line in f.splitlines():
        #     match = re.match(r"(\w+):\s*(.+)", line)
        #     if match:
        #         data[match.group(1)] = match.group(2)

        return data if data else None

    def extract_version(self, data_by_sample: Dict) -> Union[str, None]:
        """
        Extract software version from parsed data.

        Args:
            data_by_sample: Dictionary of parsed data by sample

        Returns:
            Version string or None if not found
        """
        # TODO: Implement version extraction
        # Example:
        # for sample_data in data_by_sample.values():
        #     if "version" in sample_data:
        #         return sample_data["version"]
        return None

    def add_general_stats(self, data_by_sample: Dict):
        """
        Add key metrics to the general statistics table.

        Args:
            data_by_sample: Dictionary of parsed data by sample
        """
        headers = {
            "metric_name": {
                "title": "Metric Display Name",
                "description": "Detailed description of this metric",
                "scale": "RdYlGn",  # Color scale (RdYlGn, Blues, etc.) or False
                "format": "{:,.0f}",  # Number formatting
                "suffix": " reads",  # Optional suffix
                "min": 0,  # Optional minimum value for scale
            }
        }

        self.general_stats_addcols(data_by_sample, headers)

    def add_{toolname}_plots(self, data_by_sample: Dict):
        """
        Create plot sections for the module.

        Args:
            data_by_sample: Dictionary of parsed data by sample
        """
        # Example: Bar graph
        # self.add_section(
        #     name="Results Summary",
        #     anchor="{toolname}-results",
        #     description="Overview of {toolname} results",
        #     plot=bargraph.plot(
        #         data_by_sample,
        #         pconfig={
        #             "id": "{toolname}_bargraph",
        #             "title": "{Display Name}: Results",
        #             "ylab": "Count",
        #         },
        #     ),
        # )

        # Example: Line graph
        # line_data = {}
        # for s_name, data in data_by_sample.items():
        #     line_data[s_name] = data.get("distribution", {})
        #
        # self.add_section(
        #     name="Quality Distribution",
        #     anchor="{toolname}-quality",
        #     description="Quality score distribution",
        #     plot=linegraph.plot(
        #         line_data,
        #         pconfig={
        #             "id": "{toolname}_linegraph",
        #             "title": "{Display Name}: Quality Scores",
        #             "xlab": "Position",
        #             "ylab": "Quality Score",
        #         },
        #     ),
        # )

        # Example: Table
        # headers = {
        #     "metric1": {"title": "Metric 1", "format": "{:,.0f}"},
        #     "metric2": {"title": "Metric 2", "format": "{:,.2f}%"},
        # }
        #
        # self.add_section(
        #     name="Detailed Statistics",
        #     anchor="{toolname}-stats",
        #     description="Detailed per-sample statistics",
        #     plot=table.plot(data_by_sample, headers),
        # )

        pass
