""" MultiQC module to parse output from MEGAHIT """


import logging
import re


from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """MEGAHIT module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="MEGAHIT",
            anchor="megahit",
            href="https://github.com/voutcn/megahit",
            info="is an ultra-fast and memory-efficient NGS assembler",
            doi="10.1093/bioinformatics/btv033",
        )

        data = dict()
        for f in self.find_log_files("megahit", filehandles=True):
            for line in f["f"]:
                self.add_data_source(f, f["s_name"])
                if " - MEGAHIT v" in line:
                    version = line.split(" - MEGAHIT v")[1].strip()
                    self.add_software_version(version, sample=f["s_name"])
                elif " - Memory used: " in line:
                    mem = line.split(" - Memory used: ")[1].strip()
                    try:
                        mem_mb = float(mem) / 1024 / 1024
                    except ValueError:
                        mem_mb = None
                    data[f["s_name"]] = {"megahit_mem": mem_mb}
                elif " - ALL DONE. Time elapsed: " in line:
                    time = line.split(" - ALL DONE. Time elapsed: ")[1].strip().replace(" seconds", "")
                    data[f["s_name"]]["megahit_time"] = time
                elif " contigs, total " in line:
                    # Parse line like this:
                    # 2022-06-13 14:22:47 - 64408 contigs, total 46322050 bp, min 200 bp, max 94571 bp, avg 719 bp, N50 831 bp
                    regex = re.compile(
                        r"(\d+) contigs, total (\d+) bp, min (\d+) bp, max (\d+) bp, avg (\d+) bp, N50 (\d+) bp"
                    )
                    match = regex.search(line)
                    if match:
                        data[f["s_name"]]["megahit_contigs"] = int(match.group(1))
                        data[f["s_name"]]["megahit_bases"] = int(match.group(2))
                        data[f["s_name"]]["megahit_min_contig"] = int(match.group(3))
                        data[f["s_name"]]["megahit_max_contig"] = int(match.group(4))
                        data[f["s_name"]]["megahit_avg_contig"] = int(match.group(5))
                        data[f["s_name"]]["megahit_n50"] = int(match.group(6))

        # Filter to strip out ignored sample names
        data = self.ignore_samples(data)
        if len(data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data)} reports")

        # Write parsed report data to a file
        self.write_data_file(data, f"multiqc_{self.anchor}")

        headers = {
            "megahit_contigs": {
                "title": "Contigs",
                "description": "Number of contigs",
                "min": 0,
                "scale": "Greens",
                "format": "{:,d}",
            },
            "megahit_bases": {
                "title": f"Bases, {config.base_count_prefix}",
                "description": f"Total number of bases ({config.base_count_desc})",
                "shared_key": "base_count",
                "min": 0,
                "scale": "Blues",
                "hidden": True,
            },
            "megahit_min_contig": {
                "title": "Min contig",
                "description": "Minimum contig length",
                "min": 0,
                "scale": "YlGn",
                "suffix": "&nbsp;bp",
                "format": "{:,d}",
                "hidden": True,
            },
            "megahit_max_contig": {
                "title": "Max contig",
                "description": "Maximum contig length",
                "min": 0,
                "scale": "YlGn",
                "suffix": "&nbsp;bp",
                "format": "{:,d}",
                "hidden": True,
            },
            "megahit_avg_contig": {
                "title": "Avg",
                "description": "Average contig length",
                "min": 0,
                "scale": "YlGn",
                "suffix": "&nbsp;bp",
                "format": "{:,d}",
            },
            "megahit_n50": {
                "title": "N50",
                "description": "N50 contig length",
                "min": 0,
                "scale": "Greens",
                "suffix": "&nbsp;bp",
                "format": "{:,d}",
            },
            "megahit_mem": {
                "title": "Memory",
                "description": "Memory used, MB",
                "min": 0,
                "suffix": "&nbsp;MB",
                "scale": "Reds",
                "hidden": True,
            },
            "megahit_time": {
                "title": "Time",
                "description": "Time elapsed",
                "min": 0,
                "suffix": "&nbsp;s",
                "scale": "Greys",
                "hidden": True,
            },
        }

        # Basic Stats Table
        self.general_stats_addcols(data, headers)

        # Separate table with all metrics visible
        for k, v in headers.items():
            v["hidden"] = False
        self.add_section(
            name="Run statistics",
            anchor="megahit-stats",
            description="Metrics and run statistics from MEGAHIT run logs.",
            plot=table.plot(
                data,
                headers,
                pconfig={"id": "megahit-stats-table", "title": "Run statistics"},
            ),
        )
