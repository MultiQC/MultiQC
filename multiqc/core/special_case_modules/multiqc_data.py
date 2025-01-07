"""Special case MultiQC module to load multiqc_data.json"""

import json
import logging
from pathlib import Path
from typing import Union

from multiqc import report
from multiqc.base_module import BaseMultiqcModule
from multiqc.types import Anchor

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="MultiQC Data",
            anchor=Anchor("multiqc_data"),
            info="loads multiqc_data.json",
        )

        for f in self.find_log_files("multiqc_data"):
            parse_data_json(Path(f["root"]) / f["fn"])


def parse_data_json(path: Union[str, Path]):
    """
    Try find multiqc_data.json in the given directory, and load it into the report.

    @param path: Path to the directory containing multiqc_data.json or the path to the file itself.
    """
    path = Path(path)
    assert path.suffix == ".json"
    log.info(f"Loading data from {path}")
    try:
        with path.open("r") as f:
            data = json.load(f)

        for mod, sections in data["report_data_sources"].items():
            log.info(f"Loaded module {mod}")
            for section, sources in sections.items():
                for sname, source in sources.items():
                    report.data_sources[mod][section][sname] = source
        for id, plot_dump in data["report_plot_data"].items():
            log.info(f"Loaded plot {id}")
            report.plot_data[id] = plot_dump
    except (json.JSONDecodeError, KeyError) as e:
        log.error(f"Error loading data from multiqc_data.json: {e}")
