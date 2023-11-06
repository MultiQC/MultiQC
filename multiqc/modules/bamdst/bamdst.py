""" MultiQC module to parse output from Bamdst """


import logging
from typing import Dict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils import config

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Bamdst",
            anchor="bamdst",
            href="https://https://github.com/shiquan/bamdst",
            info="""is a lightweight tool to stat the depth coverage of target regions of bam file(s)""",
            # doi="", # No DOI
        )

        # Find and load reports
        data_by_sample = dict()
        for f in self.find_log_files("bamdst/coverage"):
            data = self._parse_coverage_report(f)
            if data:
                data_by_sample[f["s_name"]] = data

        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} reports")

        self.general_stats(data_by_sample)

    def _parse_coverage_report(self, f: Dict) -> Dict:
        """
        Parse one coverage report.
        """
        data = dict()
        s_name = None
        version = None
        for line in f["f"].split("\n"):
            line = line.strip()
            if line.startswith("## Files : "):
                # Split a line like "## Files : example/test1.bam example/test2.bam "
                # Cannot split by spaces because the file names can contain spaces as well
                path_line = line.split(":")[1].strip().replace(".bam", "<EXT>").replace(".cram", "<EXT>")
                names = [self.clean_s_name(sn) for sn in path_line.split("<EXT>") if sn]
                f["s_name"] = "-".join(names)  # concatenating file names
                continue
            if line.startswith("## Version : "):
                version = line.split(":")[1].strip()
                continue
            fields = line.split("\t")
            if not len(fields) == 2:
                continue
            metric, value = fields
            if value.endswith("%"):
                value = value.strip("%")
                try:
                    value = float(value)
                except ValueError:
                    continue
                value = value / 100.0
            else:
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        continue
            data[metric] = value

        if s_name and version and data:
            self.add_software_version(version, sample=s_name)
            self.add_data_source(f, s_name)
        return data

    def general_stats(self, data_by_sample):
        headers = {
            "[Total] Mapped Reads": {
                "title": f"Mapped ({config.read_count_prefix})",
                "description": f"Mapped reads ({config.read_count_desc})",
                "scale": "RdYlGn",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "hidden": True,
            },
            "[Total] Fraction of Mapped Reads": {
                "title": "Mapped",
                "description": f"Fraction of mapped reads in all reads ({config.read_count_desc})",
                "scale": "PuBu",
                "min": 0.0,
                "max": 100.0,
                "suffix": "%",
                "modify": lambda x: x * 100.0,
                "hidden": True,
            },
            "[Target] Target Reads": {
                "title": f"Target reads ({config.read_count_prefix})",
                "description": f"Reads mapped on target ({config.read_count_desc})",
                "scale": "RdYlGn",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
            "[Target] Fraction of Target Reads in all reads": {
                "title": "% in all",
                "description": "Fraction of target reads in all reads",
                "scale": "PuBu",
                "min": 0.0,
                "max": 100.0,
                "suffix": "%",
                "modify": lambda x: x * 100.0,
                "hidden": True,
            },
            "[Target] Fraction of Target Reads in mapped reads": {
                "title": "% in mapped",
                "description": "Fraction of target reads in mapped reads",
                "scale": "PuBu",
                "min": 0.0,
                "max": 100.0,
                "suffix": "%",
                "modify": lambda x: x * 100.0,
                "hidden": True,
            },
            "[Target] Average depth": {
                "title": "Avg depth",
                "description": "Average depth, on target",
                "min": 0,
                "suffix": "X",
                "scale": "BuPu",
            },
            "[Target] Len of region": {
                "title": "Region len",
                "description": "Length of target region",
                "min": 0,
                "format": "{:,d}",
                "suffix": "&nbsp;bp",
                "scale": "Greys",
                "hidden": True,
            },
        }

        coverage_metrics = []
        for s_name, d in data_by_sample.items():
            for k in d.keys():
                if k.startswith("[Target] Coverage (") and k not in coverage_metrics:
                    coverage_metrics.append(k)
        for m in coverage_metrics:
            short = m.replace("[Target] Coverage (", "").replace(")", "")
            headers[m] = {
                "title": f"{short}",
                "description": f"Target coverage ({short})",
                "min": 0,
                "suffix": "%",
                "format": "{:.2f}",
                "scale": "RdYlGn",
                "modify": lambda x: x * 100.0,
                "hidden": short != ">0x",
            }

        self.general_stats_addcols(data_by_sample, headers)
