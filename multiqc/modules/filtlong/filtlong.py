import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Filtlong",
            anchor="filtlong",
            href="https://github.com/rrwick/Filtlong",
            info="A tool for filtering long reads by quality.",
            # doi="", # No DOI
        )

        # Find and load reports
        self.filtlong_data = dict()

        # Find all files for filtlong
        for f in self.find_log_files("filtlong", filehandles=True):
            self.parse_logs(f)

        self.filtlong_data = self.ignore_samples(self.filtlong_data)

        if len(self.filtlong_data) == 0:
            raise UserWarning

        log.info(f"Found {len(self.filtlong_data)} reports")

        # Write data to file
        self.write_data_file(self.filtlong_data, "filtlong")
        self.filtlong_general_stats()

    def parse_logs(self, f):
        """Parsing Logs. Note: careful of ANSI formatting log"""
        for l in f["f"]:
            # Find the valid metric
            if "target:" in l:
                self.add_data_source(f)
                if f["s_name"] in self.filtlong_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                target_bases = l.lstrip().split(" ")[1]
                # Remove . thousand separators - see ewels/MultiQC#1843
                target_bases = float(target_bases.replace(".", ""))
                self.filtlong_data[f["s_name"]] = {"Target bases": target_bases}

            elif "keeping" in l and f["s_name"] in self.filtlong_data:
                bases_kept = l.lstrip().split(" ")[1]
                # Remove . thousand separators - see ewels/MultiQC#1843
                bases_kept = float(bases_kept.replace(".", ""))
                self.filtlong_data[f["s_name"]]["Bases kept"] = bases_kept

            elif "fall below" in l and f["s_name"] in self.filtlong_data:
                log.debug(f"{f['s_name']}: reads already fall below target after filtering")

    def filtlong_general_stats(self):
        """Filtlong General Stats Table"""
        headers = OrderedDict()
        headers["Target bases"] = {
            "title": "Target bases ({})".format(config.read_count_prefix),
            "description": "Keep only the best reads up to this many total bases ({})".format(config.read_count_desc),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["Bases kept"] = {
            "title": "Bases kept ({})".format(config.read_count_prefix),
            "description": "Bases kept ({})".format(config.read_count_desc),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        }

        self.general_stats_addcols(self.filtlong_data, headers)
