# import csv
# from typing import Dict, Union
# import logging
# from multiqc.base_module import ModuleNoSamplesFound

# from multiqc.base_module import BaseMultiqcModule


# log = logging.getLogger(__name__)

# def parse_file(fh) -> Dict[str, Union[float, int]]:
#     data = {}
#     reader = csv.DictReader(fh)
#     for row in reader:
#         print("Row:", row)  # shows all columns in the row
#         metric_name = row.get("Metric Name", "").strip()
#         value = row.get("Metric Value", "").strip()
#         print("Parsed metric:", metric_name, "| Value:", value)


#         if not metric_name or not value:
#             continue

#         value = value.replace(",", "").replace("%", "")
#         try:
#             data[metric_name] = float(value)
#         except ValueError:
#             log.warning(f"Couldn't convert value '{value}' for metric '{metric_name}'")
#     return data


# # def parse_file(fh) -> Dict[str, Union[float, int]]:
# #         """Parses an open CSV filehandle and returns a dict of metrics"""
# #         data = {}
# #         reader = csv.DictReader(fh)
# #         for row in reader:
# #             for key, val in row.items():
# #                 if key.lower() == "sample":
# #                     continue  # skip sample name column
# #                 try:
# #                     data[key] = float(val)
# #                 except (ValueError, TypeError):
# #                     log.warning(f"Couldn't convert value '{val}' for key '{key}'")
# #         return data

# class MultiqcModule(BaseMultiqcModule):
#     def __init__(self):
#         super(MultiqcModule, self).__init__(
#           name="CSVFlex",
#           anchor="custom_csv",
#           info="parses in files of the type .csv and extracts metrics",
#         )
#         log.info("Hello World!")

#         data_by_sample: Dict[str, Dict[str, Union[float, int]]] = {}

#         for f in self.find_log_files("csv", filehandles=True):
#             parsed = parse_file(f["f"])
#             print(parsed)

#             if parsed:
#                 s_name = f["s_name"]
#                 print("Sample name detected:", f["s_name"])
#                 log.info(parsed)
#                 log.info(s_name)
#                 if self.is_ignore_sample(s_name):
#                     continue
                
#                 if s_name in data_by_sample:
#                     log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
#                 data_by_sample[s_name] = parsed
#                 self.add_data_source(f)

#         self.add_software_version(None)
    
            
#         data_by_sample = self.ignore_samples(data_by_sample)
        
        
#         if len(data_by_sample) == 0:
#             raise ModuleNoSamplesFound
#         log.info(f"Found {len(data_by_sample)} reports")
        
#         self.write_data_file(data_by_sample, "multiqc_csv")

#         self.add_section(
#             name="CSV Results Summary",
#             anchor="csv_summary",
#             description="Summary of parsed results from CSV files."
#         )


# import csv
# from typing import Dict, Union
# import logging
# from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# log = logging.getLogger(__name__)

# def parse_file(fh) -> Dict[str, Union[float, int]]:
#     data = {}
#     reader = csv.DictReader(fh)
#     for row in reader:
#         log.debug(f"Row: {row}")

#         metric_name = row.get("Metric Name", "").strip()
#         value = row.get("Metric Value", "").strip()

#         if not metric_name or not value:
#             continue

#         # Remove commas and percentage signs
#         value = value.replace(",", "").replace("%", "")
#         try:
#             data[metric_name] = float(value)
#         except ValueError:
#             log.warning(f"Couldn't convert value '{value}' for metric '{metric_name}'")
#     return data


# class MultiqcModule(BaseMultiqcModule):
#     def __init__(self):
#         super().__init__(
#             name="CSVFlex",
#             anchor="custom_csv",
#             info="Parses .csv files and extracts metrics for summary reporting."
#         )
#         log.info("Hello World!")

#         data_by_sample: Dict[str, Dict[str, Union[float, int]]] = {}

#         for f in self.find_log_files("csv", filehandles=True):
#             log.info(f"Processing file: {f}")
#             log.info(f"Found log file: {f['fn']} (sample name: {f['s_name']})")
#             # parsed = parse_file(f["f"])
#             parsed = parse_file(f)
#             log.debug(f"Parsed data: {parsed}")


#             if parsed:
#                 s_name = f.get("s_name") or f.get("fn") or "unknown_sample"
#                 log.info(f"Sample name: {s_name}")

#                 if s_name

#                 if self.is_ignore_sample(s_name):
#                     continue

#                 if s_name in data_by_sample:
#                     log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

#                 data_by_sample[s_name] = parsed
#                 self.add_data_source(f)

#         self.add_software_version(None)
#         data_by_sample = self.ignore_samples(data_by_sample)

#         if len(data_by_sample) == 0:
#             log.warning("No sample data found after processing CSV files.")
#             raise ModuleNoSamplesFound

#         self.write_data_file(data_by_sample, "multiqc_csv")

#         self.add_section(
#             name="CSV Results Summary",
#             anchor="csv_summary",
#             description="Summary of parsed results from CSV files."
#         )

import csv
from typing import Dict, Union
import logging
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


def parse_file(fh) -> Dict[str, Dict[str, Union[float, int]]]:
    """
    Parses a CSV file and returns a dictionary mapping sample names (Group Name)
    to their associated metrics.
    """
    data: Dict[str, Dict[str, Union[float, int]]] = {}

    reader = csv.DictReader(fh)

    for row in reader:
        metric_name = row.get("Metric Name", "").strip()
        value = row.get("Metric Value", "").strip()
        s_name = row.get("Group Name", "").strip() or "global"

        log.debug(f"Parsing row - Sample: {s_name}, Metric: {metric_name}, Value: {value}")

        if not metric_name or not value:
            continue

        # Clean numeric values
        value = value.replace(",", "").replace("%", "")

        try:
            val = float(value)
        except ValueError:
            log.warning(f"Couldn't convert value '{value}' for metric '{metric_name}'")
            continue

        # Initialize sample if not already present
        if s_name not in data:
            data[s_name] = {}

        data[s_name][metric_name] = val

    return data


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="CSVFlex",
            anchor="custom_csv",
            info="Parses .csv files and extracts summary metrics grouped by sample.",
        )
        log.info("Custom CSV Module initialized!")

        data_by_sample: Dict[str, Dict[str, Union[float, int]]] = {}

        for f in self.find_log_files("custom_csv", filehandles=True):
            log.info(f"Processing file: {f['fn']} (initial sample guess: {f['s_name']})")

            parsed_samples = parse_file(f["f"])

            if not parsed_samples:
                log.warning(f"No data parsed from: {f['fn']}")
                continue

            for s_name, metrics in parsed_samples.items():
                log.info(f"Parsed sample: {s_name} with {len(metrics)} metrics")
                if self.is_ignore_sample(s_name):
                    continue

                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

                data_by_sample[s_name] = metrics
                self.add_data_source(f)

        self.add_software_version(None)

        data_by_sample = self.ignore_samples(data_by_sample)

        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} sample(s) with metrics.")

        self.write_data_file(data_by_sample, "multiqc_csv")

        self.add_section(
            name="CSV Results Summary",
            anchor="csv_summary",
            description="Summary of parsed results from CSV files."
        )
