import logging

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module takes summary statistics of number of long reads filtered and displays them in the General Stats table.

    #### Bases Kept

    Sometimes, the Filtlong log message contains this:

    ```
    Filtering long reads
      target: 123456789 bp
      reads already fall below target after filtering

    Outputting passed long reads
    ```

    In these cases we cannot say for sure how many bases were kept. As such, this field is left blank.
    If you have a better solution, please suggest in an issue or pull request.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Filtlong",
            anchor="filtlong",
            href="https://github.com/rrwick/Filtlong",
            info="Filters long reads by quality.",
            extra="It can take a set of long reads and produce a smaller, better subset. It uses both read length "
            "(longer is better) and read identity (higher is better) when choosing which reads pass the filter.",
            # doi="", # No DOI
        )

        # Find and load reports
        self.filtlong_data = dict()

        # Find all files for filtlong
        for f in self.find_log_files("filtlong", filehandles=True):
            self.parse_logs(f)

        self.filtlong_data = self.ignore_samples(self.filtlong_data)

        if len(self.filtlong_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.filtlong_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write data to file
        self.write_data_file(self.filtlong_data, "filtlong")
        self.filtlong_general_stats()

    def parse_logs(self, f):
        """Parsing Logs. Note: careful of ANSI formatting log"""
        for line in f["f"]:
            # Find the valid metric
            if "target:" in line:
                self.add_data_source(f)
                if f["s_name"] in self.filtlong_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                target_bases = line.lstrip().split(" ")[1]
                # Remove . thousand separators - see MultiQC/MultiQC#1843
                target_bases = float(target_bases.replace(".", ""))
                self.filtlong_data[f["s_name"]] = {"Target bases": target_bases}

            elif "keeping" in line and f["s_name"] in self.filtlong_data:
                bases_kept = line.lstrip().split(" ")[1]
                # Remove . thousand separators - see MultiQC/MultiQC#1843
                bases_kept = float(bases_kept.replace(".", ""))
                self.filtlong_data[f["s_name"]]["Bases kept"] = bases_kept

            elif "fall below" in line and f["s_name"] in self.filtlong_data:
                log.debug(f"{f['s_name']}: reads already fall below target after filtering")

    def filtlong_general_stats(self):
        """Filtlong General Stats Table"""
        headers = {
            "Target bases": {
                "title": f"Target bases ({config.read_count_prefix})",
                "description": "Keep only the best reads up to this many total bases ({})".format(
                    config.read_count_desc
                ),
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "Bases kept": {
                "title": f"Bases kept ({config.read_count_prefix})",
                "description": f"Bases kept ({config.read_count_desc})",
                "scale": "Purples",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
        }

        self.general_stats_addcols(self.filtlong_data, headers)
