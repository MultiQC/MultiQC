from __future__ import annotations

import csv
import logging
from typing import Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """mlst scans contig files against PubMLST typing schemes to determine sequence types.

    To generate files compatible with MultiQC, run `mlst` with the `--full` option:

    ```bash
    mlst --full contigs.fasta > mlst_report.tsv
    # Or CSV format:
    mlst --full --csv contigs.fasta > mlst_report.csv
    ```

    The default output without `--full` omits the header row and status/score metrics
    required for MultiQC parsing.
    """

    def __init__(self) -> None:
        super().__init__(
            name="MLST",
            anchor="mlst",
            href="https://github.com/tseemann/mlst",
            info="Scans contig files against PubMLST typing schemes.",
            doi="",
            license="GPL-2.0",
            license_url="https://github.com/tseemann/mlst/blob/master/LICENSE",
        )

        self.mlst_data: dict[str, dict[str, Any]] = {}

        for f in self.find_log_files("mlst"):
            self.parse_mlst_log(f)

        # Filter out ignored samples via MultiQC core configuration
        self.mlst_data = self.ignore_samples(self.mlst_data)

        if not self.mlst_data:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.mlst_data)} MLST reports")

        # Track tool version (required by MultiQC linting even if unknown)
        self.add_software_version(None)

        # General Statistics table
        self.add_general_stats()

        # Dedicated module section with detailed table
        self.add_mlst_section()

        # Export parsed data files (JSON/TSV/YAML)
        self.write_data_file(self.mlst_data, "multiqc_mlst")

    def parse_mlst_log(self, f: dict[str, Any]) -> None:
        """Parse mlst TSV or CSV output."""
        lines = f["f"].splitlines()
        if not lines:
            return

        # Determine delimiter (TSV default or CSV)
        first_line = lines[0]
        delimiter = "," if ("," in first_line and "\t" not in first_line) else "\t"

        reader = csv.DictReader(lines, delimiter=delimiter)
        for row in reader:
            raw_filename = row.get("FILE", "").strip()
            if not raw_filename or raw_filename == "-":
                sample_name = f["s_name"]
            else:
                sample_name = self.clean_s_name(raw_filename, f)

            scheme = row.get("SCHEME", "").strip()
            st = row.get("ST", "").strip()
            status = row.get("STATUS", "").strip()
            raw_score = row.get("SCORE", "").strip()
            alleles = row.get("ALLELES", "").strip()

            try:
                score: float | None = float(raw_score)
            except (ValueError, TypeError):
                score = None

            # Convert "-" or empty values to None so MultiQC renders them
            # as blank, uncolored cells rather than assigning a category color.
            clean_st = st if (st and st != "-") else None
            clean_scheme = scheme if (scheme and scheme != "-") else None
            clean_status = status if (status and status != "-") else None
            clean_alleles = alleles if (alleles and alleles != "-") else None

            if sample_name in self.mlst_data:
                log.debug(f"Duplicate sample name found! Overwriting: {sample_name}")

            self.mlst_data[sample_name] = {
                "scheme": clean_scheme,
                "st": clean_st,
                "status": clean_status,
                "score": score,
                "alleles": clean_alleles,
            }

    def add_general_stats(self) -> None:
        """Add MLST summary columns to the General Statistics table."""
        headers: dict[str, dict[str, Any]] = {
            "scheme": {
                "title": "Scheme",
                "description": "PubMLST scheme name",
                "scale": "Set1",
            },
            "st": {
                "title": "ST",
                "description": "PubMLST Sequence Type (ST)",
                "scale": "Set2",
            },
            "status": {
                "title": "MLST Status",
                "description": "Prediction confidence status (PERFECT, OK, NOVEL, etc.)",
                "scale": False,
            },
            "score": {
                "title": "MLST Score",
                "description": "MLST prediction score (0-100)",
                "min": 0,
                "max": 100,
                "scale": "Blues",
                "format": "{:,.0f}",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.mlst_data, headers)

    def add_mlst_section(self) -> None:
        """Add dedicated section containing full prediction details."""
        headers: dict[str, dict[str, Any]] = {
            "scheme": {
                "title": "Scheme",
                "description": "Auto-detected PubMLST scheme",
                "scale": "Set1",
            },
            "st": {
                "title": "Sequence Type (ST)",
                "description": "Assigned Sequence Type identifier",
                "scale": "Set2",
            },
            "status": {
                "title": "Status",
                "description": "Quality assessment of the MLST call",
                "scale": False,
            },
            "score": {
                "title": "Score",
                "description": "Confidence score out of 100",
                "min": 0,
                "max": 100,
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "alleles": {
                "title": "Alleles",
                "description": "Identified allele call profiles",
                "scale": False,
            },
        }

        table_config = {
            "id": "mlst_predictions_table",
            "title": "MLST Sequence Typing Results",
            "save_file": False,
        }

        self.add_section(
            name="MLST Typings",
            anchor="mlst-typings",
            description=(
                "Scanned contig files evaluated against PubMLST schemes. "
                "Requires <code>mlst</code> to be run with the <code>--full</code> option."
            ),
            helptext="""
            The `mlst` module parses outputs created with the `--full` option (`FILE`, `SCHEME`, `ST`, `STATUS`, `SCORE`, `ALLELES`).
            Standard `mlst` output without `--full` lacks the header row and completeness metrics required for MultiQC parsing.
            """,
            plot=table.plot(self.mlst_data, headers=headers, pconfig=table_config),
        )
        