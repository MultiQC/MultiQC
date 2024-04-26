#!/usr/bin/env python

import json
import logging
import os
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Hostile Module"""

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Hostile",
            anchor="hostile",
            href="https://github.com/bede/hostile",
            info="Short and Long host reads removal tools.",
            doi="10.1093/bioinformatics/btad728",
        )

        self.parse_data = dict()
        for f in self.find_log_files("hostile", filehandles=True):
            self.parse_logs(f)
        
        # add version
        self.add_software_version(self.version,f)
        
        # Ignore samples
        self.parse_data = self.ignore_samples(self.parse_data)

        # Check if no matching log files were found
        if len(self.parse_data) == 0:
            raise ModuleNoSamplesFound
        
        # Number of files found
        log.info(f"Found {len(self.parse_data)} reports")
        
        # Write hostile data.
        self.write_data_file(self.parse_data, "multiqc_hostile")

        # Plot
        self.hostile_plot()

    def parse_logs(self, json_file):
        """
        Parsing json_file
        """
        try:
            parse_data = json.load(json_file["f"])

        except json.JSONDecodeError:
            log.warning(f"Could not parse JSON file {json_file['f']}")
            return

        if len(parse_data) > 0:
            self.add_data_source(json_file)
        
        s_name = self.clean_s_name(json_file["fn"])
        self.parse_data[s_name] = parse_data

        #print(self.parse_data[s_name][0])
        if "version" in self.parse_data[s_name][0]:
            self.version = self.parse_data[s_name][0]['version']
        else:
            self.version = "0.4.0"

    def hostile_plot(self):
        """
        Extract the data
        Plot the barplot
        """
        data = {}

        for f_name, values in self.parse_data.items():
            s_name = values[0]["fastq1_in_name"].split(".")[0]
            database = os.path.basename(values[0]["index"])
            data[s_name] = {"Cleaned reads": values[0]["reads_out"], "Host reads": values[0]["reads_removed"]}
        
        ## categories
        cats = ["Cleaned reads", "Host reads"]

        pconfig = {
            "title": "Hostile: Reads Filtered",
            "id": "he_reads_plots",
            "ylab": "# Reads",
            "plot_type": "bargraph"
        }

        self.add_section(
            name="Reads Filtering",
            anchor="hostile-reads",
            description=f"This plot shows the number of cleaned reads vs host-reads per sample (database index: {database}).",
            plot=bargraph.plot(data, cats, pconfig),
        )