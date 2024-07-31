import logging
import xml.etree.cElementTree
from typing import Dict

from multiqc import config
from multiqc.modules.ngsbits.utils import parse_qcml_by
from multiqc.plots import table

log = logging.getLogger(__name__)


def check_paired_end(qcml_contents):
    """Check if both R1 and R2 are present in input files."""
    found_r1 = False
    found_r2 = False
    root = xml.etree.cElementTree.fromstring(qcml_contents)

    for el in root.findall(".//{http://www.prime-xs.eu/ms/qcml}metaDataParameter"):
        if el.attrib["name"] == "source file":
            if "R1" in el.attrib["value"]:
                found_r1 = True
            if "R2" in el.attrib["value"]:
                found_r2 = True

    return found_r1 and found_r2


def parse_reports(self):
    """Find ngs-bits ReadQC reports and parse their data"""

    readqc: Dict = dict()
    readqc_keys: Dict = dict()

    for f in self.find_log_files("ngsbits/readqc"):
        values, params = parse_qcml_by(f["f"], "qualityParameter")
        is_pe = check_paired_end(f["f"])

        if len(values) > 0:
            if f["s_name"] in readqc:
                log.debug(f'Duplicate sample name found! Overwriting: {f["s_name"]}')
            self.add_data_source(f, section="readqc")
            readqc[f["s_name"]] = values
            readqc_keys.update(params)
            readqc[f["s_name"]].update(
                {
                    "cluster count": readqc[f["s_name"]]["read count"] / (1 + is_pe),
                    "paired-end": "yes" if is_pe else "no",
                }
            )

    # Filter to strip out ignored sample names
    readqc = self.ignore_samples(readqc)

    if len(readqc) == 0:
        return 0

    # Write to file
    self.write_data_file(readqc, "multiqc_ngsbits_readqc")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Convert numbers given in megabases to bases
    readqc_keys["bases sequenced"] = ("Bases sequenced in total.", "")
    for _, kv in readqc.items():
        kv["bases sequenced"] = kv["bases sequenced (MB)"] * 1e6
        kv.pop("bases sequenced (MB)")

    # Add cluster count
    readqc_keys["cluster count"] = ("Clusters sequenced in total.", "")
    readqc_keys["paired-end"] = ("Whether input files were paired-end sequences.", "")

    # Improve table headers
    readqc_keys_table = {key: {"description": value[0]} for key, value in readqc_keys.items()}
    readqc_keys_table["read count"].update(
        {
            "title": "Reads",
            "format": "{:,.2f}",
            "shared_key": "read_count",
            "scale": "Purples",
            "placement": 10,
        }
    )
    readqc_keys_table["cluster count"].update(
        {
            "title": "Clusters",
            "format": "{:,.2f}",
            "shared_key": "read_count",
            "scale": "Purples",
            "placement": 20,
        }
    )
    readqc_keys_table["bases sequenced"].update(
        {
            "title": "Bases",
            "format": "{:,.2f}",
            "shared_key": "base_count",
            "scale": "Blues",
            "placement": 30,
        }
    )
    readqc_keys_table["gc content %"].update(
        {
            "title": "GC content",
            "suffix": "%",
            "format": "{:,.2f}",
            "max": 100,
            "scale": "Spectral",
            "placement": 40,
        }
    )
    readqc_keys_table["Q20 read %"].update(
        {
            "title": "Q20",
            "suffix": "%",
            "format": "{:,.2f}",
            "max": 100,
            "scale": "Reds",
            "placement": 50,
        }
    )
    readqc_keys_table["Q30 base %"].update(
        {
            "title": "Q30",
            "suffix": "%",
            "format": "{:,.2f}",
            "max": 100,
            "scale": "Oranges",
            "placement": 60,
        }
    )
    readqc_keys_table["read length"].update(
        {
            "title": "Read length",
            "suffix": "bp",
            "format": "{:,.0f}",
            "scale": "Greens",
            "placement": 70,
        },
    )
    readqc_keys_table["no base call %"].update(
        {
            "title": "No base call",
            "suffix": "%",
            "format": "{:,.2f}",
            "floor": 1,
            "scale": "BuGn",
            "placement": 80,
        }
    )

    # overview table with all values
    self.add_section(
        name="ReadQC",
        anchor="ngsbits-readqc",
        description='<a href="https://github.com/imgag/ngs-bits/blob/master/doc/tools/ReadQC.md" target="_blank">ReadQC</a>'
        " calculates QC metrics on unprocessed NGS reads.",
        plot=table.plot(
            readqc,
            readqc_keys_table,
            pconfig={
                "namespace": "ReadQC",
                "id": "ngsbits_readqc_table",
                "title": "ngs-bits: ReadQC Summary",
            },
        ),
    )

    return len(readqc)
