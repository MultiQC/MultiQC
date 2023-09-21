"""
MultiQC submodule to parse output from Iso-Seq refine logs
"""

import csv
import json
import logging
import math
from collections import defaultdict
from pprint import pprint

from multiqc.plots import table

log = logging.getLogger(__name__)


class RefineMixin:
    def parse_refine_log(self):
        data = dict()

        for f in self.find_log_files("isoseq/refine-json", filehandles=True):
            data[f["s_name"]] = json.load(f["f"])
            self.add_data_source(f)

        for f in self.find_log_files("isoseq/refine-csv", filehandles=True):
            counts = defaultdict(int)
            sums = defaultdict(int)
            squared_sums = defaultdict(int)
            mins = defaultdict(lambda: float("inf"))
            maxs = defaultdict(lambda: float("-inf"))
            strand_counts = defaultdict(int)
            primer_counts = defaultdict(int)

            # Read the CSV file
            reader = csv.DictReader(f["f"])
            for row in reader:
                # Process columns of interest
                for column in ["fivelen", "threelen", "polyAlen", "insertlen"]:
                    value = float(row[column])
                    counts[column] += 1
                    sums[column] += value
                    squared_sums[column] += value**2
                    mins[column] = min(mins[column], value)
                    maxs[column] = max(maxs[column], value)

                # Count occurrences of 'strand' and 'primer' values
                strand_counts[row["strand"]] += 1
                primer_counts[row["primer"]] += 1

            d = dict()
            if counts:
                # Compute mean, standard deviation, min, and max and store them in data
                for column in ["fivelen", "threelen", "polyAlen", "insertlen"]:
                    mean = sums[column] / counts[column]
                    variance = (squared_sums[column] / counts[column]) - mean**2
                    std_dev = math.sqrt(variance) if variance > 0 else 0.0

                    d[f"min_{column}"] = mins[column]
                    d[f"mean_{column}"] = mean
                    d[f"std_{column}"] = std_dev
                    d[f"max_{column}"] = maxs[column]
                d["strand_counts"] = dict(strand_counts)
                d["primer_counts"] = dict(primer_counts)
                if d:
                    if f["s_name"] in data:
                        data[f["s_name"]].update(d)
                    else:
                        data[f["s_name"]] = d
                    self.add_data_source(f)
        self.write_data_file(data, "multiqc_isoseq_refine_report")
        return data

    def add_general_stats_refine(self, data):
        gstats_data = {}
        for s_name, attrs in data.items():
            gstats_data[s_name] = {}
            if "num_reads_fl" not in attrs:
                pprint(data)
            gstats_data[s_name]["num_reads_fl"] = attrs["num_reads_fl"]
            gstats_data[s_name]["num_reads_flnc"] = attrs["num_reads_flnc"]
            gstats_data[s_name]["num_reads_flnc_polya"] = attrs["num_reads_flnc_polya"]

        headers = dict()
        headers["num_reads_fl"] = {
            "title": "# full-length reads",
            "description": "Number of CCS where both primers have been detected",
            "suffix": " fl",
            "scale": "spectral",
        }
        headers["num_reads_flnc"] = {
            "title": "# non-chimeric full-length reads",
            "description": "Number of non-chimeric CCS where both primers have been detected",
            "suffix": " flnc",
            "scale": "RdYlGn",
        }
        headers["num_reads_flnc_polya"] = {
            "title": "# poly(A) free non-chimeric full-length reads",
            "description": "Number of non-chimeric CCS where both primers have been detected and the poly(A) tail has been removed",
            "suffix": " flnc_polya",
            "scale": "spectral",
        }
        self.general_stats_addcols(gstats_data, headers, namespace="refine")

    def add_table_refine(self, data_by_sample):
        headers = dict()
        plot_data = {s_name: dict() for s_name in data_by_sample}

        # 5' primer length
        for s_name, data in data_by_sample.items():
            plot_data[s_name]["min_fivelen"] = data["min_fivelen"]
            plot_data[s_name]["mean_fivelen"] = data["mean_fivelen"]
            plot_data[s_name]["std_fivelen"] = data["std_fivelen"]
            plot_data[s_name]["max_fivelen"] = data["max_fivelen"]

        headers["min_fivelen"] = {
            "title": "Min 5' primer length",
            "description": "The minimum 5' primer length in base pair",
            "scale": "spectral",
        }
        headers["mean_fivelen"] = {
            "title": "Mean 5' primer length",
            "description": "The mean 5' primer length in base pair",
            "scale": "RdYlGn",
        }
        headers["std_fivelen"] = {
            "title": "Std of 5' primer length",
            "description": "The standard deviation of 5' primer length in base pair",
            "scale": "spectral",
        }
        headers["max_fivelen"] = {
            "title": "Max 5' primer length",
            "description": "The maximum 5' primer length in base pair",
            "scale": "RdYlGn",
        }

        # 3' primer length
        for s_name, data in data_by_sample.items():
            plot_data[s_name]["min_threelen"] = data["min_threelen"]
            plot_data[s_name]["mean_threelen"] = data["mean_threelen"]
            plot_data[s_name]["std_threelen"] = data["std_threelen"]
            plot_data[s_name]["max_threelen"] = data["max_threelen"]

        headers["min_threelen"] = {
            "title": "Min 3' primer length",
            "description": "The minimum 3' primer length in base pair",
            "scale": "spectral",
        }
        headers["mean_threelen"] = {
            "title": "Mean 3' primer length",
            "description": "The mean 3' primer length in base pair",
            "scale": "RdYlGn",
        }
        headers["std_threelen"] = {
            "title": "Std of 3' primer length",
            "description": "The standard deviation of 3' primer length in base pair",
            "scale": "spectral",
        }
        headers["max_threelen"] = {
            "title": "Max 3' primer length",
            "description": "The maximum 3' primer length in base pair",
            "scale": "RdYlGn",
        }

        # Poly-A tail length
        for s_name, data in data_by_sample.items():
            plot_data[s_name]["min_polyAlen"] = data["min_polyAlen"]
            plot_data[s_name]["mean_polyAlen"] = data["mean_polyAlen"]
            plot_data[s_name]["std_polyAlen"] = data["std_polyAlen"]
            plot_data[s_name]["max_polyAlen"] = data["max_polyAlen"]

        headers["min_polyAlen"] = {
            "title": "Min polyA tail length",
            "description": "The minimum polyA tail length in base pair",
            "scale": "spectral",
        }
        headers["mean_polyAlen"] = {
            "title": "Mean polyA tail length",
            "description": "The mean polyA tail length in base pair",
            "scale": "RdYlGn",
        }
        headers["std_polyAlen"] = {
            "title": "Std of polyA tail length",
            "description": "The standard deviation of polyA tail length in base pair",
            "scale": "spectral",
        }
        headers["max_polyAlen"] = {
            "title": "Max polyA tail length",
            "description": "The maximum polyA tail length in base pair",
            "scale": "RdYlGn",
        }

        # Insert length
        for s_name, data in data_by_sample.items():
            plot_data[s_name]["min_insertlen"] = data["min_insertlen"]
            plot_data[s_name]["mean_insertlen"] = data["mean_insertlen"]
            plot_data[s_name]["std_insertlen"] = data["std_insertlen"]
            plot_data[s_name]["max_insertlen"] = data["max_insertlen"]

        headers["min_insertlen"] = {
            "title": "Min insert length",
            "description": "The minimum insert length in base pair",
            "scale": "spectral",
        }
        headers["mean_insertlen"] = {
            "title": "Mean insert length",
            "description": "The mean insert length in base pair",
            "scale": "RdYlGn",
        }
        headers["std_insertlen"] = {
            "title": "Std of insert length",
            "description": "The standard deviation of insert length in base pair",
            "scale": "spectral",
        }
        headers["max_insertlen"] = {
            "title": "Max insert length",
            "description": "The maximum insert length in base pair",
            "scale": "RdYlGn",
        }

        config = {
            "id": "isoseq_refine_table",
            "title": "Iso-Seq refine",
        }

        self.add_section(
            name="Iso-Seq refine",
            anchor="insert-refine-stats",
            description="Iso-Seq refine statistics",
            helptext="""
            The .report.csv file contains information about 5' prime and 3' primers length, insert length,
            poly(A) length, and couple of primers detected for each CCS.
            The table presents min, max, mean, standard deviation for each parameter.
            """,
            plot=table.plot(plot_data, headers, config),
        )
