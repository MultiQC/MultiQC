import logging
from typing import Dict

from multiqc import config
from multiqc.plots import table


log = logging.getLogger(__name__)


def parse_reports(self):
    """Find ngs-bits MappingQC reports and parse their data"""

    mappingqc: Dict = dict()
    mappingqc_keys: Dict = dict()

    for f in self.find_log_files("ngsbits/mappingqc"):
        values, params = self.parse_qcml_by(f["f"], "qualityParameter")

        if len(values) > 0:
            if f["s_name"] in mappingqc:
                log.debug(f'Duplicate sample name found! Overwriting: {f["s_name"]}')
            self.add_data_source(f, section="mappingqc")
            mappingqc[f["s_name"]] = values
            mappingqc_keys.update(params)

    # Filter to strip out ignored sample names
    mappingqc = self.ignore_samples(mappingqc)

    if len(mappingqc) == 0:
        return 0

    # Write to file
    self.write_data_file(mappingqc, "multiqc_ngsbits_mappingqc")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Convert numbers given in megabases to bases
    mappingqc_keys["bases usable"] = ("Bases usable in total.", "")
    for _, kv in mappingqc.items():
        kv["bases usable"] = kv["bases usable (MB)"] * 1e6
        kv.pop("bases usable (MB)")

    # Improve table headers
    mappingqc_keys_table = {key: {"title": key, "description": value[0]} for key, value in mappingqc_keys.items()}

    mappingqc_keys_table["trimmed base %"].update({"suffix": "%", "format": "{:,.2f}", "floor": 1, "scale": "PuBu"})
    mappingqc_keys_table["clipped base %"].update({"suffix": "%", "format": "{:,.2f}", "floor": 1, "scale": "PuRd"})
    mappingqc_keys_table["mapped read %"].update({"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "Reds"})
    mappingqc_keys_table["bases usable"].update(
        {
            "suffix": config.base_count_prefix,
            "format": "{:,.2f}",
            "modify": lambda x: x * config.base_count_multiplier,
            "scale": "Greens",
        }
    )
    # always available, even without target file
    mappingqc_keys_table["on-target read %"].update(
        {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "Purples"}
    )

    # only available if duplicates marked
    try:
        mappingqc_keys_table["duplicate read %"].update(
            {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "YlOrRd"}
        )
    except KeyError:
        pass

    # only available if paired-end
    try:
        mappingqc_keys_table["properly-paired read %"].update(
            {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "GnBu"}
        )
        mappingqc_keys_table["insert size"].update({"suffix": "bp", "format": "{:,.2f}", "scale": "RdYlGn"})
    except KeyError:
        pass

    # only available if human
    try:
        mappingqc_keys_table["SNV allele frequency deviation"].update(
            {"suffix": "", "format": "{:,.2f}", "floor": 0, "ceiling": 10, "minRange": 10, "scale": "Greys"}
        )
    except KeyError:
        pass

    # only available if target file provided
    coverage_values = (10, 20, 30, 50, 100, 200, 500)
    try:
        mappingqc_keys_table["target region read depth"].update({"suffix": "x", "format": "{:,.2f}"})
        for x in coverage_values:
            mappingqc_keys_table["target region {:d}x %".format(x)].update(
                {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "YlGn"}
            )
    except KeyError:
        pass

    # overview table with all values
    self.add_section(
        name="MappingQC",
        anchor="ngsbits-mappingqc",
        description='<a href="https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC.md" target="_blank">MappingQC</a>'
        " calculates QC metrics on mapped NGS reads.",
        plot=table.plot(
            mappingqc,
            mappingqc_keys_table,
            pconfig={
                "namespace": "ngsbits_mappingqc",
                "id": "ngsbits_mappingqc_table",
                "title": "ngs-bits: MappingQC Summary",
            },
        ),
    )

    return len(mappingqc)
