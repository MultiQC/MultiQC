import logging

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Sargasso",
            anchor="sargasso",
            href="http://biomedicalinformaticsgroup.github.io/Sargasso/",
            info="Separates mixed-species RNA-seq reads according to their species of origin.",
            doi="10.1038/s41596-018-0029-2",
        )

        # Find and load any Sargasso reports
        self.sargasso_data = dict()
        self.sargasso_files = list()
        self.sargasso_keys = list()  # header keys

        for f in self.find_log_files("sargasso"):
            self.parse_sargasso_logs(f)
            self.sargasso_files.append(f)
            self.add_data_source(f)

        # log.info('Removing ignored samples...')
        self.sargasso_data = self.ignore_samples(self.sargasso_data)

        if len(self.sargasso_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.sargasso_files)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.sargasso_data, "multiqc_sargasso")

        # Basic Stats Table
        self.sargasso_stats_table()

        # Assignment bar plot
        self.add_section(plot=self.sargasso_chart())

    def parse_sargasso_logs(self, f):
        """Parse the sargasso log file."""
        species_name = list()
        items = list()
        is_first_line = True
        for line in f["f"].splitlines():
            s = line.split(",")
            # Check that this actually is a Sargasso file
            if is_first_line and s[0] != "Sample":
                return None

            if len(s) < 7:
                continue
            if is_first_line:
                # prepare header
                is_first_line = False
                header = s
                for i in header[1:]:
                    # find out what species included
                    sname = i.split("-")[-1]
                    if sname not in species_name:
                        species_name.append(sname)
                    # find out what is being counted
                    kname = "-".join(i.split("-")[-3:-1])
                    if kname not in items:
                        items.append(kname)
            else:
                # start sample lines.
                sample_name = s.pop(0)

                chunk_by_species = [s[i : i + len(items)] for i in range(0, len(s), len(items))]
                for idx, v in enumerate(chunk_by_species):
                    # adding species name to the sample name for easy interpretation
                    new_sample_name = "_".join([sample_name, species_name[idx]])

                    # Clean up sample name
                    new_sample_name = self.clean_s_name(new_sample_name, f)

                    if new_sample_name in self.sargasso_data.keys():
                        log.debug(f"Duplicate sample name found! Overwriting: {new_sample_name}")

                    try:
                        self.sargasso_data[new_sample_name] = dict(zip(items, map(int, v)))
                    except ValueError:
                        pass

        self.sargasso_keys = items

        for idx, f_name in enumerate(self.sargasso_data.keys()):
            # Reorganised parsed data for this sample
            # Collect total READ count number
            self.sargasso_data[f_name]["Total"] = 0
            for key, value in list(self.sargasso_data[f_name].items()):  # iter on both keys and values
                if key.endswith("Reads"):
                    self.sargasso_data[f_name]["Total"] += value

            # Calculate the percent aligned if we can
            try:
                self.sargasso_data[f_name]["sargasso_percent_assigned"] = (
                    float(self.sargasso_data[f_name]["Assigned-Reads"]) / float(self.sargasso_data[f_name]["Total"])
                ) * 100.0
            except (KeyError, ZeroDivisionError):
                pass

    def sargasso_stats_table(self):
        """Take the parsed stats from the sargasso report and add them to the
        basic stats table at the top of the report"""

        headers = {
            "sargasso_percent_assigned": {
                "title": "% Assigned",
                "description": "Sargasso % Assigned reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "Assigned-Reads": {
                "title": f"{config.read_count_prefix} Assigned",
                "description": f"Sargasso Assigned reads ({config.read_count_desc})",
                "min": 0,
                "scale": "PuBu",
                "modify": lambda x: float(x) * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(self.sargasso_data, headers)

    def sargasso_chart(self):
        """Make the sargasso plot"""

        # Config for the plot
        config = {
            "id": "sargasso_assignment_plot",
            "title": "Sargasso: Assigned Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        # We only want to plot the READs at the moment
        return bargraph.plot(self.sargasso_data, [name for name in self.sargasso_keys if "Reads" in name], config)
