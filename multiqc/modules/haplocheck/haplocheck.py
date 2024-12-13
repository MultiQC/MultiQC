from typing import Dict, Union
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Haplocheck",
            anchor="haplocheck",
            href="https://github.com/genepi/haplocheck/",
            info="Haplocheck detects in-sample contamination in mtDNA or WGS sequencing studies by analyzing the mitchondrial content.",
            doi=["10.1101/gr.256545.119"],
        )
        haplocheck_data: Dict = dict()

        for f in self.find_log_files("haplocheck"):
            haplocheck_data = self.parse_logs(f)

        haplocheck_data = self.ignore_samples(haplocheck_data)

        if len(haplocheck_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(haplocheck_data)} Haplocheck reports")

        # Write data to file
        self.write_data_file(haplocheck_data, "multiqc_haplocheck")

        self.add_software_version(None)

        # Headers dictionary for Contamination Data
        headers = {
            "Contamination Status": {
                "title": "Contamination Status",
                "description": "Indicates whether contamination was detected in the sample.",
                "scale": False,
            },
            "Contamination Level": {
                "title": "Contamination Level",
                "description": "Estimated level of contamination (if applicable).",
                "min": 0,
                "scale": "Oranges",
                "format": "{:,.2f}",
            },
            "Distance": {
                "title": "Distance",
                "description": "Genomic distance value associated with contamination.",
                "min": 0,
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "Sample Coverage": {
                "title": "Sample Coverage",
                "description": "The total coverage of the sample sequence.",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.0f}",
            },
        }

        self.add_section(
            name="Haplocheck",
            anchor="haplocheck-section",
            description="Haplocheck detects in-sample contamination in mtDNA or WGS sequencing studies by analyzing the mitchondrial content.",
            plot=table.plot(
                data=haplocheck_data,
                headers=headers,
                pconfig={
                    "id": "haplocheck-table",
                    "title": "Haplocheck",
                },
            ),
        )

        for header in headers.values():
            header["hidden"] = True
        headers["Contamination Status"]["hidden"] = False

        self.general_stats_addcols(haplocheck_data, headers)

    def parse_logs(self, f: str) -> Dict[str, Dict[str, Union[float, str]]]:
        parsed_data: Dict[str, Dict[str, Union[float, str]]] = {}
        file_content = f["f"]
        lines = file_content.strip().splitlines()

        # Assuming the first line contains headers
        headers = [header.strip().replace('"', "") for header in lines[0].strip().split("\t")[1:]]

        for line in lines[1:]:
            columns = line.strip().split("\t")
            sample_name = columns[0]  # First column is the sample name
            s_name = sample_name.replace("_haplocheck", "").replace('"', "")  # Clean sample name and remove quotes
            values = columns[1:]
            self.add_data_source(f, s_name)

            # Map headers to values for this sample
            parsed_data[s_name] = {
                key: (
                    float(value.strip().replace('"', ""))
                    if value.strip().replace('"', "").replace(".", "", 1).isdigit()
                    else value.strip().replace('"', "")
                )
                for key, value in zip(headers, values)
            }

        return parsed_data
