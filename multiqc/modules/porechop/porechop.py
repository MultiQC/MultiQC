import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Porechop",
            anchor="porechop",
            href="https://github.com/rrwick/Porechop",
            info="a tool for finding and removing adapters from Oxford Nanopore reads.",
            # doi="",  # No DOI available
        )

        # Find and load reports
        self.porechop_data = dict()

        # Find all files for porechop
        for f in self.find_log_files("porechop", filehandles=True):
            self.parse_logs(f)

        self.porechop_data = self.ignore_samples(self.porechop_data)

        if len(self.porechop_data) == 0:
            raise ModuleNoSamplesFound

        log.info("Found {} reports".format(len(self.porechop_data)))

        # Write data to file
        self.write_data_file(self.porechop_data, "porechop")

        self.porechop_general_stats()
        if max(len(v) for v in self.porechop_data.values()) > 1:
            self.start_trim_barplot()
            self.end_trim_barplot()
            self.middle_split_barplot()
        self.no_adapters_found()

    def parse_logs(self, logfile):
        """Parsing Logs. Note: careful of ANSI formatting log"""

        def get_float(val):
            """Get float from string"""
            val = val.replace(",", "")
            try:
                return float(val)
            except ValueError:
                return val

        file_content = logfile["f"]
        for l in file_content:
            ## Find line after loading reads, and remove suffixes for sample name
            if "Loading reads" in l:
                s_name = next(file_content).rstrip()
                s_name = self.clean_s_name(s_name, logfile)
                self.add_data_source(logfile, s_name=s_name)
                if s_name in self.porechop_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.porechop_data[s_name] = {}

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, s_name)

            ## Find each valid metric, clean up for plain integer
            # 10,000 reads loaded
            if "reads loaded" in l:
                reads_loaded = re.search(r"([\d,]+)\s*reads loaded", l)
                if reads_loaded:
                    self.porechop_data[s_name]["Input Reads"] = get_float(reads_loaded.group(1))

            # 7,100 / 10,000 reads had adapters trimmed from their start (425,196 bp removed)
            if "reads had adapters trimmed from their start" in l:
                adapter_start = re.search(
                    r"([\d,]+)\s*/\s*([\d,]+)\s*reads had adapters trimmed from their start \(([\d,]+) bp removed\)", l
                )
                if adapter_start:
                    self.porechop_data[s_name]["Start Trimmed"] = get_float(adapter_start.group(1))
                    self.porechop_data[s_name]["Start Trimmed Total"] = get_float(adapter_start.group(2))
                    self.porechop_data[s_name]["Start Trimmed (bp)"] = get_float(adapter_start.group(3))
                    self.porechop_data[s_name]["Start Untrimmed"] = (
                        self.porechop_data[s_name]["Start Trimmed Total"] - self.porechop_data[s_name]["Start Trimmed"]
                    )
                    try:
                        self.porechop_data[s_name]["Start Trimmed Percent"] = (
                            self.porechop_data[s_name]["Start Trimmed"]
                            / self.porechop_data[s_name]["Start Trimmed Total"]
                        ) * 100
                    except ZeroDivisionError:
                        pass

            # 4,849 / 10,000 reads had adapters trimmed from their end (283,192 bp removed)'
            if "reads had adapters trimmed from their end" in l:
                end_trimmed = re.search(
                    r"([\d,]+)\s*/\s*([\d,]+)\s*reads had adapters trimmed from their end \(([\d,]+) bp removed\)", l
                )
                if end_trimmed:
                    self.porechop_data[s_name]["End Trimmed"] = get_float(end_trimmed.group(1))
                    self.porechop_data[s_name]["End Trimmed Total"] = get_float(end_trimmed.group(2))
                    self.porechop_data[s_name]["End Trimmed (bp)"] = get_float(end_trimmed.group(3))
                    self.porechop_data[s_name]["End Untrimmed"] = (
                        self.porechop_data[s_name]["End Trimmed Total"] - self.porechop_data[s_name]["End Trimmed"]
                    )
                    try:
                        self.porechop_data[s_name]["End Trimmed Percent"] = (
                            self.porechop_data[s_name]["End Trimmed"]
                            / self.porechop_data[s_name]["Start Trimmed Total"]
                            * 100
                        )
                    except ZeroDivisionError:
                        pass

            # 7 / 10,000 reads were split based on middle adapters
            if "reads were split based on middle adapters" in l:
                split_stats = re.search(r"([\d,]+)\s*/\s*([\d,]+)\s*reads were split based on middle adapters", l)
                if split_stats:
                    self.porechop_data[s_name]["Middle Split"] = get_float(split_stats.group(1))
                    self.porechop_data[s_name]["Middle Split Total"] = get_float(split_stats.group(2))
                    self.porechop_data[s_name]["Middle Not-Split"] = (
                        self.porechop_data[s_name]["Middle Split Total"] - self.porechop_data[s_name]["Middle Split"]
                    )
                    try:
                        self.porechop_data[s_name]["Middle Split Percent"] = (
                            self.porechop_data[s_name]["Middle Split"]
                            / self.porechop_data[s_name]["Middle Split Total"]
                            * 100
                        )
                    except ZeroDivisionError:
                        pass

    def porechop_general_stats(self):
        """Porechop General Stats Table"""
        headers = OrderedDict()
        headers["Input Reads"] = {
            "title": "Input Reads ({})".format(config.read_count_prefix),
            "description": "Number of reads loaded into Porechop ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["Start Trimmed"] = {
            "title": "Start Trimmed ({})".format(config.read_count_prefix),
            "description": "Number of reads that had adapters trimmed from the start ({})".format(
                config.read_count_prefix
            ),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Start Trimmed Percent"] = {
            "title": "Start Trimmed",
            "description": "Percent of reads that had adapters trimmed from the start",
            "suffix": "%",
            "max": 100,
            "scale": "RdYlGn",
        }
        headers["End Trimmed"] = {
            "title": "End Trimmed ({})".format(config.read_count_prefix),
            "description": "Number of reads that had adapters trimmed from the end ({})".format(
                config.read_count_prefix
            ),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["End Trimmed Percent"] = {
            "title": "End Trimmed",
            "description": "Percent of reads that had adapters trimmed from the end",
            "suffix": "%",
            "max": 100,
            "scale": "RdYlGn",
        }
        headers["Middle Split"] = {
            "title": "Middle Split ({})".format(config.read_count_prefix),
            "description": "Number of reads split based on middle adapters ({})".format(config.read_count_prefix),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }
        headers["Middle Split Percent"] = {
            "title": "Middle Split",
            "description": "Percent of reads that were split based on middle adapters",
            "suffix": "%",
            "max": 100,
            "scale": "RdYlGn",
        }

        self.general_stats_addcols(self.porechop_data, headers)

    def start_trim_barplot(self):
        """Barplot of number of reads adapter trimmed at read start"""
        cats = OrderedDict()
        cats["Start Trimmed"] = {"name": "Start Trimmed", "color": "#7cb5ec"}
        cats["Start Untrimmed"] = {"name": "Start Untrimmed", "color": "#f7a35c"}
        config = {
            "id": "porechop-starttrim-barplot",
            "title": "Porechop: Read Start Adapter Timmed",
            "ylab": "Read Counts",
        }
        self.add_section(
            name="Reads adapter-trimmed read start",
            anchor="porechop-starttrim",
            description="Shows the number of reads that had adapters removed from read start.",
            plot=bargraph.plot(self.porechop_data, cats, config),
        )

    def end_trim_barplot(self):
        """Barplot of number of reads adapter trimmed at read end"""
        cats = OrderedDict()
        cats["End Trimmed"] = {"name": "End Trimmed", "color": "#7cb5ec"}
        cats["End Untrimmed"] = {"name": "End Untrimmed", "color": "#f7a35c"}
        config = {
            "id": "porechop-endtrim-barplot",
            "title": "Porechop: Read End Adapter Timmed",
            "ylab": "Read Counts",
        }
        self.add_section(
            name="Reads adapter-trimmed read end",
            anchor="porechop-endtrim",
            description="Shows the number of reads that had adapters removed from read end.",
            plot=bargraph.plot(self.porechop_data, cats, config),
        )

    def middle_split_barplot(self):
        """Barplot of number of reads adapter trimmed at read end"""
        cats = OrderedDict()
        cats["Middle Split"] = {"name": "Split Reads", "color": "#7cb5ec"}
        cats["Middle Not-Split"] = {"name": "Unsplit Reads", "color": "#f7a35c"}
        config = {
            "id": "porechop-middlesplit-barplot",
            "title": "Porechop: Middle Split",
            "ylab": "Read Counts",
        }
        self.add_section(
            name="Middle split reads",
            anchor="porechop-middlesplit",
            description="Shows the number of reads that were split due to adapter being present in middle of read.",
            plot=bargraph.plot(self.porechop_data, cats, config),
        )

    def no_adapters_found(self):
        """Show any samples that did not have any trimming"""
        no_adapters = []
        for s_name in self.porechop_data:
            if len(self.porechop_data[s_name]) == 1:
                no_adapters.append(s_name)
        if len(no_adapters):
            self.add_section(
                name="No adapters found",
                anchor="porechop-noadapters",
                description="The following samples did not have any adapters found - output reads were unchanged from input reads:",
                content=f"""
                    <ul>
                        <li><code>{'</code></li><li><code>'.join(no_adapters)}</code></li>
                    </ul>
                """,
            )
