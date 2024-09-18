import logging
from multiqc.plots import table

log = logging.getLogger(__name__)


from multiqc.base_module import BaseMultiqcModule
 
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='ngs-bits SampleGender',
            anchor='samplegender',
            href="https://github.com/imgag/ngs-bits",
            info="SampleGender determines the gender of a sample from the BAM/CRAM file.",
            doi=[""],
        )

    def parse_reports(self):
        """Find and parse ngsbits SampleGender TSV output files."""

        self.samplegender_data = dict()
        for f in self.find_log_files("ngsbits_samplegender", filehandles=True):
            parsed_data = self.samplegender_parse_reports(f["f"])
            if parsed_data:
                if f["s_name"] in self.samplegender_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="samplegender")
                self.samplegender_data[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.samplegender_data = self.ignore_samples(self.samplegender_data)

        n_reports_found = len(self.samplegender_data)
        if n_reports_found == 0:
            return 0

        log.info(f"Found {n_reports_found} SampleGender reports")

        # Add software version information if available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.samplegender_data, "multiqc_ngsbits_samplegender.txt")

        # Add SampleGender Table
        self.samplegender_table(self.samplegender_data)

    def samplegender_parse_reports(self, f):
        """Parse the ngsbits SampleGender TSV file for all samples."""
        
        data = dict()
        headers = None
        
        for line in f:
            if line.startswith("#"):  # Skip header or comment lines
                continue
            parts = line.strip().split("\t")
            if headers is None:
                headers = parts  # Use the first non-comment line as headers
                continue
            if len(parts) >= 5:  # Ensure there are at least five columns
                sample_data = {
                    "file": parts[0],
                    "gender": parts[1],
                    "reads_chry": int(parts[2]),
                    "reads_chrx": int(parts[3]),
                    "ratio_chry_chrx": float(parts[4]),
                }
                data[parts[0]] = sample_data  # Use the file name as the key
    
        return data


    def samplegender_table(self, data):
        """Build a table from the parsed SampleGender output."""
    
        headers = {
            "file": {
                "title": "File Name",
                "description": "The name of the BAM file analyzed.",
                "namespace": "ngsbits",
            },
            "gender": {
                "title": "Predicted Gender",
                "description": "The predicted gender based on chromosome read ratios.",
                "namespace": "ngsbits",
            },
            "reads_chry": {
                "title": "Reads on ChrY",
                "description": "The number of reads mapped to the Y chromosome.",
                "namespace": "ngsbits",
                "min": 0,
            },
            "reads_chrx": {
                "title": "Reads on ChrX",
                "description": "The number of reads mapped to the X chromosome.",
                "namespace": "ngsbits",
                "min": 0,
            },
            "ratio_chry_chrx": {
                "title": "ChrY/ChrX Ratio",
                "description": "The ratio of reads mapped to ChrY vs ChrX.",
                "namespace": "ngsbits",
                "min": 0,
                "format": "{:.4f}",
            },
        }
        table_html = table.plot(data, headers, {"id": "ngsbits_samplegender_table", "title": "ngsbits - SampleGender"})
        return table_html
