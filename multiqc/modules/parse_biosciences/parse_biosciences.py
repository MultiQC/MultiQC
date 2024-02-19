""" MultiQC module to parse output from Parse Biosciences pipelines """
from typing import Dict, Union

import logging

from multiqc.plots import table
from multiqc.utils import config

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Parse Biosciences",
            anchor="parse_biosciences",
            href="https://www.parsebiosciences.com/",
            info="provides single cell RNA-Seq pipelines",
            # doi="",
        )

        split_pipe_samples = self.parse_split_pipe()
        if split_pipe_samples > 0:
            log.info(f"Found {split_pipe_samples} samples in Parse Biosciences Split pipeline outputs")

        if split_pipe_samples == 0:
            raise ModuleNoSamplesFound

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    CSV_HEADERS = {
        "sample_well_count": {
            "title": "Sample wells",
            "description": "Number of wells with cells",
            "format": "{:,d}",
            "min": 0,
        },
        "number_of_cells": {
            "rid": "bpsplit_cells",
            "title": "Cells",
            "description": "Number of cells detected",
            "format": "{:,d}",
            "min": 0,
            "scale": "RdPu",
        },
        "mean_reads_per_cell": {
            "title": "Reads per cell",
            "description": "Mean reads per cell",
            "format": "{:,.0f}",
            "min": 0,
            "scale": "Blues",
        },
        "number_of_reads": {
            "rid": "bpsplit_reads",
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Number of reads ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "scale": "YlGn",
        },
        "number_of_tscps": {
            "title": f"{config.read_count_prefix} TSCP",
            "description": f"Number of TSCP ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
        },
        "sequencing_saturation": {
            "title": "Saturation",
            "description": "Sequencing saturation",
            "format": "{:.2%}",
            "min": 0,
            "max": 1,
            "scale": "YlOrRd",
        },
        "valid_barcode_fraction": {
            "title": "Valid BC",
            "description": "Valid barcode fraction",
            "format": "{:.2%}",
            "min": 0,
            "max": 1,
            "scale": "Spectral",
        },
        "bc1_Q30": {
            "title": "BC1 Q30",
            "description": "BC1 Q30",
            "format": "{:.2%}",
            "min": 0,
            "max": 1,
            "hidden": True,
        },
        "bc2_Q30": {
            "title": "BC2 Q30",
            "description": "BC2 Q30",
            "format": "{:.2%}",
            "min": 0,
            "max": 1,
        },
        "bc3_Q30": {
            "title": "BC3 Q30",
            "description": "BC3 Q30",
            "format": "{:.2%}",
            "min": 0,
            "max": 1,
            "hidden": True,
        },
        "cDNA_Q30": {
            "title": "cDNA Q30",
            "description": "cDNA Q30",
            "format": "{:.2%}",
            "min": 0,
            "max": 1,
            "hidden": True,
        },
        "polyN_Q30": {
            "title": "PolyN Q30",
            "description": "PolyN Q30",
            "format": "{:.2%}",
            "min": 0,
            "max": 1,
            "hidden": True,
        },
    }

    def parse_split_pipe(self) -> int:
        data_by_sample = dict()
        for f in self.find_log_files("parse_biosciences/split_pipe_csv", filehandles=True):
            d_by_sn = self._parse_split_pipe_csv(f["f"])
            if not d_by_sn:
                return 0
            self.add_data_source(f)
            if any(sn in data_by_sample for sn in d_by_sn):
                log.debug(f"Duplicate sample names found in {f['fn']}! Overwriting: {', '.join(d_by_sn.keys())}")
            data_by_sample.update(d_by_sn)

        self.add_section(
            name="Split pipe",
            anchor="parse-biosciences-split-pipe-stats",
            description="Summary QC metrics from Parse Biosciences Split pipeline",
            plot=table.plot(
                data_by_sample,
                MultiqcModule.CSV_HEADERS,
                pconfig={
                    "namespace": "Split pipe",
                    "id": "parse-biosciences-split-pipe-stats-table",
                },
            ),
        )
        data_by_sample = self.ignore_samples(data_by_sample)
        self.write_data_file(data_by_sample, "multiqc_bp_split")
        return len(data_by_sample)

    @staticmethod
    def _parse_split_pipe_csv(fh) -> Dict[str, Dict[str, Union[float, int]]]:
        header = next(fh).strip().split(",")
        if header[0] != "statistic":
            return {}
        sample_names = header[1:]
        data_by_sample = {sn: {} for sn in sample_names}
        for line in fh:
            row = line.strip().split(",")
            metric = row[0]
            if metric not in MultiqcModule.CSV_HEADERS:
                continue
            for si, val in enumerate(row[1:]):
                try:
                    val = float(val)
                except ValueError:
                    log.debug(f"Could not parse {val} as float on line {line.strip()} in {fh.name}")
                    continue
                if MultiqcModule.CSV_HEADERS[metric].get("format") == "{:,d}":
                    val = int(val)
                data_by_sample[sample_names[si]][metric] = val
        return data_by_sample
