import os
import pandas as pd
from multiqc.utils import report
from multiqc.utils.config import config
from multiqc.modules.base_module import BaseMultiqcModule


def parse_reports(self):
    yield_files = self.find_log_files(
        fn_match=r'_collapsed_grouped\.duplex_yield_metrics\.txt$',
        contents_match='fraction',
        filecontents_res=[r"fraction"]
    )

    if not yield_files:
        self.log.warning("❌ No yield metric files matched.")
        return 0
    else:
        self.log.info(f"✅ Found {len(yield_files)} yield metric files.")

    parsed_data = {}

    for f in yield_files:
        self.log.debug(f"📄 Matched file: {f['fn']} at path: {f['f']}")

        with open(f['f'], 'r') as fh:
            lines = fh.readlines()

        self.log.debug(f"🧪 First line: {lines[0].strip()}")
        if len(lines) > 1:
            self.log.debug(f"🧪 Second line: {lines[1].strip()}")

        df = pd.read_csv(f['f'], sep="\t")
        total_yield_fraction = df["fraction"].sum()

        sample_name = self.clean_s_name(f["s_name"])
        parsed_data[sample_name] = {"Total Yield Fraction": total_yield_fraction}

    # Add data table section
    self.add_section(
        name="Fgbio Duplex Yield Metrics",
        anchor="collapsed_duplex_metrics",
        description="Summed yield fractions from fgbio CollectDuplexSeqMetrics files.",
        plot_type="table",
        data=parsed_data,
        headers={
            "Total Yield Fraction": {
                "title": "Total Yield Fraction",
                "description": "Sum of 'fraction' column from the duplex_yield_metrics file",
                "format": "{:,.4f}"
            }
        }
    )

    return len(parsed_data)
