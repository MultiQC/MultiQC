import csv
from typing import Dict, Union
import logging
from multiqc.base_module import ModuleNoSamplesFound

from multiqc.base_module import BaseMultiqcModule


log = logging.getLogger(__name__)

def parse_file(fh) -> Dict[str, Union[float, int]]:
    data = {}
    reader = csv.DictReader(fh)
    for row in reader:
        metric_name = row.get("Metric Name", "").strip()
        value = row.get("Metric Value", "").strip()

        if not metric_name or not value:
            continue

        value = value.replace(",", "").replace("%", "")
        try:
            data[metric_name] = float(value)
        except ValueError:
            log.warning(f"Couldn't convert value '{value}' for metric '{metric_name}'")
    return data


# def parse_file(fh) -> Dict[str, Union[float, int]]:
#         """Parses an open CSV filehandle and returns a dict of metrics"""
#         data = {}
#         reader = csv.DictReader(fh)
#         for row in reader:
#             for key, val in row.items():
#                 if key.lower() == "sample":
#                     continue  # skip sample name column
#                 try:
#                     data[key] = float(val)
#                 except (ValueError, TypeError):
#                     log.warning(f"Couldn't convert value '{val}' for key '{key}'")
#         return data

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
          name="CSVModule",
          anchor="custom_csv",
          info="parses in files of the type .csv and extracts metrics",
        )
        log.info("Hello World!")

        data_by_sample: Dict[str, Dict[str, Union[float, int]]] = {}

        for f in self.find_log_files("csv", filehandles=True):
            parsed = parse_file(f["f"])

            if parsed:
                s_name = f["s_name"]
                log.info(parsed)
                log.info(s_name)
                if self.is_ignore_sample(s_name):
                    continue
                
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                data_by_sample[s_name] = parsed
                self.add_data_source(f)

        self.add_software_version(None)
    
            
        data_by_sample = self.ignore_samples(data_by_sample)
        
        
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")
        
        self.write_data_file(data_by_sample, "multiqc_csv")

        self.add_section(
            name="CSV Results Summary",
            anchor="csv_summary",
            description="Summary of parsed results from CSV files."
        )



    