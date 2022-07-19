from multiqc.modules.base_module import BaseMultiqcModule
import logging
import re

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Porechop",
            anchor="porechop",
            href="https://github.com/rrwick/Porechop",
            info="A tool for finding and removing adapters from Oxford Nanopore reads.",
        )

        # Find and load reports
        self.porechop_data = dict()

        # Find all files for porechop
        for f in self.find_log_files("porechop", filehandles=True):
            self.parse_logs(f)

        self.porechop_data = self.ignore_samples(self.porechop_data)

        if len(self.porechop_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.porechop_data)))

        # Write data to file
        self.write_data_file(self.porechop_data, "porechop")

        print(self.porechop_data)

    def parse_logs(self, logfile):

        file_content = logfile["f"]
        for l in file_content:
            ## Find line after loading reads, and remove suffixes for sample name
            if "Loading reads" in l:
                s_name = next(file_content).rstrip()
                s_name = self.clean_s_name(s_name, logfile)
                self.add_data_source(logfile, s_name=s_name)
                self.porechop_data[s_name] = {}
            ## Find each valid metric, clean up for plain integer
            elif "reads loaded" in l:
                self.porechop_data[s_name]["input_reads"] = {}
                self.porechop_data[s_name]["input_reads"]["total"] = l.split(" ")[0]
            elif "from their start" in l:
                self.porechop_data[s_name]["trimmed_from_start"] = {}
                self.porechop_data[s_name]["trimmed_from_start"]["trimmed"] = l.split(" ")[0]
                self.porechop_data[s_name]["trimmed_from_start"]["total"] = l.split(" ")[2]
                self.porechop_data[s_name]["trimmed_from_start"]["bp"] = l.split(" ")[10].strip("(").replace(",", "")
            elif "from their end" in l:
                self.porechop_data[s_name]["trimmed_from_end"] = {}
                self.porechop_data[s_name]["trimmed_from_end"]["trimmed"] = l.split(" ")[0]
                self.porechop_data[s_name]["trimmed_from_end"]["total"] = l.split(" ")[2]
                self.porechop_data[s_name]["trimmed_from_end"]["bp"] = l.split(" ")[10].strip("(").replace(",", "")
            elif "split based on" in l:
                self.porechop_data[s_name]["split_middle"] = {}
                self.porechop_data[s_name]["split_middle"]["trimmed"] = l.split(" ")[0]
                self.porechop_data[s_name]["split_middle"]["total"] = l.split(" ")[2]
