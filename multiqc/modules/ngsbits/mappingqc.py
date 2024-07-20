import logging
from copy import copy
from typing import Dict

from multiqc import config, BaseMultiqcModule
from multiqc.modules.ngsbits.utils import parse_qcml_by
from multiqc.plots import table


log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find ngs-bits MappingQC reports and parse their data"""

    mappingqc: Dict = dict()
    mappingqc_keys: Dict = dict()

    for f in module.find_log_files("ngsbits/mappingqc"):
        values, params = parse_qcml_by(f["f"], "qualityParameter")

        if len(values) > 0:
            if f["s_name"] in mappingqc:
                log.debug(f'Duplicate sample name found! Overwriting: {f["s_name"]}')
            module.add_data_source(f, section="mappingqc")
            mappingqc[f["s_name"]] = values
            mappingqc_keys.update(params)

    # Filter to strip out ignored sample names
    mappingqc = module.ignore_samples(mappingqc)
    if len(mappingqc) == 0:
        return 0

    # Write to file
    module.write_data_file(mappingqc, "multiqc_ngsbits_mappingqc")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Convert numbers given in megabases to bases
    mappingqc_keys["bases usable"] = ("Bases usable in total.", "")
    for _, kv in mappingqc.items():
        kv["bases usable"] = kv["bases usable (MB)"] * 1e6
        kv.pop("bases usable (MB)")

    headers: Dict = dict()
    headers["bases usable"] = {
        "title": "Usable",
        "description": mappingqc_keys["bases usable"][0],
        "format": "{:,.2f}",
        "scale": "Blues",
        "shared_key": "base_count",
    }
    headers["mapped read %"] = {
        "title": "Mapped",
        "description": mappingqc_keys["mapped read %"][0],
        "suffix": "%",
        "format": "{:,.2f}",
        "max": 100,
        "scale": "Reds",
    }
    # always available, even without target file
    headers["on-target read %"] = {
        "title": "On-target",
        "description": mappingqc_keys["on-target read %"][0],
        "suffix": "%",
        "format": "{:,.2f}",
        "max": 100,
        "scale": "Purples",
    }

    # only available if duplicates marked
    if "duplicate read %" in mappingqc_keys:
        headers["duplicate read %"] = {
            "title": "Duplicates",
            "description": mappingqc_keys["duplicate read %"][0],
            "suffix": "%",
            "format": "{:,.2f}",
            "max": 100,
            "scale": "YlOrRd",
        }

    # only available if paired-end
    try:
        headers["properly-paired read %"] = {
            "title": "Properly paired",
            "description": mappingqc_keys["properly-paired read %"][0],
            "suffix": "%",
            "format": "{:,.2f}",
            "max": 100,
            "scale": "GnBu",
        }
        headers["insert size"] = {
            "title": "Insert size",
            "description": mappingqc_keys["insert size"][0],
            "suffix": "bp",
            "format": "{:,.2f}",
            "scale": "RdYlGn",
        }
    except KeyError:
        pass

    # only available if target file provided
    all_covs = (10, 20, 30, 50, 100, 200, 500)
    table_covs = (30, 100, 500)
    try:
        headers["target region read depth"] = {
            "title": "Target depth",
            "description": mappingqc_keys["target region read depth"][0],
            "suffix": "x",
            "format": "{:,.2f}",
        }
        for x in all_covs:
            headers[f"target region {x:d}x %"] = {
                "title": f"Target {x:d}x",
                "description": mappingqc_keys[f"target region {x:d}x %"][0],
                "suffix": "%",
                "format": "{:,.2f}",
                "max": 100,
                "scale": "YlGn",
                "hidden": x not in table_covs,
            }
    except KeyError:
        pass

    headers["trimmed base %"] = {
        "title": "Trimmed",
        "description": mappingqc_keys["trimmed base %"][0],
        "suffix": "%",
        "format": "{:,.2f}",
        "floor": 1,
        "scale": "PuBu",
    }
    headers["clipped base %"] = {
        "title": "Clipped",
        "description": mappingqc_keys["clipped base %"][0],
        "suffix": "%",
        "format": "{:,.2f}",
        "floor": 1,
        "scale": "PuRd",
        "hidden": True,
    }

    # only available if human
    if "SNV allele frequency deviation" in mappingqc_keys:
        headers["SNV allele frequency deviation"] = {
            "title": "SNV AF deviation",
            "description": mappingqc_keys["SNV allele frequency deviation"][0],
            "suffix": "",
            "format": "{:,.2f}",
            "floor": 0,
            "ceiling": 10,
            "minRange": 10,
            "scale": "Greys",
            "hidden": True,
        }

    # overview table with all values
    module.add_section(
        name="MappingQC",
        anchor="ngsbits-mappingqc",
        description='<a href="https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC.md" target="_blank">MappingQC</a>'
        " calculates QC metrics on mapped NGS reads.",
        plot=table.plot(
            mappingqc,
            headers,
            pconfig={
                "namespace": "MappingQC",
                "id": "ngsbits_mappingqc_table",
                "title": "ngs-bits: MappingQC Summary",
            },
        ),
    )

    gen_stats_headers = {
        "bases usable": copy(headers["bases usable"]),
        "mapped read %": copy(headers["mapped read %"]),
        "on-target read %": copy(headers["on-target read %"]),
        "duplicate read %": copy(headers["duplicate read %"]),
    }
    for x in table_covs:
        gen_stats_headers[f"target region {x:d}x %"] = copy(headers[f"target region {x:d}x %"])

    for k in gen_stats_headers:
        gen_stats_headers[k]["hidden"] = True
    gen_stats_headers["bases usable"]["hidden"] = False
    gen_stats_headers["target region 30x %"]["hidden"] = False

    module.general_stats_addcols(
        mappingqc,
        gen_stats_headers,
        namespace="MappingQC",
    )

    return len(mappingqc)
