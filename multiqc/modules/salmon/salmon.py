import json
import logging
import os

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph
from multiqc import config

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The Salmon module parses `meta_info.json`, `lib_format_counts.json` and `flenDist.txt` files, if found.

    :::note
    Note that `meta_info.json` must be within a directory called either `aux_info` or `aux` and will be ignored
    otherwise.
    :::
    """

    def __init__(self):
        super().__init__(
            name="Salmon",
            anchor="salmon",
            href="https://combine-lab.github.io/salmon/",
            info="Quantifies expression of transcripts using RNA-seq data.",
            doi="10.1038/nmeth.4197",
        )

        # Parse meta information. JSON win!
        self.salmon_meta = dict()
        for f in self.find_log_files("salmon/meta"):
            # Get the s_name from the parent directory
            if os.path.basename(f["root"]) in ["aux_info", "aux"]:
                s_name = os.path.basename(os.path.dirname(f["root"]))
                s_name = self.clean_s_name(s_name, f)
                self.salmon_meta[s_name] = {
                    metric: val
                    for metric, val in json.loads(f["f"]).items()
                    if isinstance(val, (int, float, str, list))
                }
                self.add_software_version(self.salmon_meta[s_name]["salmon_version"], s_name)

        # Parse library insert size distribution logs
        self.salmon_fld = dict()
        for f in self.find_log_files("salmon/fld"):
            # Get the s_name from the parent directory
            if os.path.basename(f["root"]) == "libParams":
                s_name = os.path.basename(os.path.dirname(f["root"]))
                s_name = self.clean_s_name(s_name, f)
                parsed = dict()
                for i, v in enumerate(f["f"].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed

        # Parse Library Format Counts information. JSON file expected
        self.salmon_lfc = dict()
        for f in self.find_log_files("salmon/lfc"):
            s_name = os.path.basename(f["root"])  # lfc file located at root folder
            s_name = self.clean_s_name(s_name, f)
            self.salmon_lfc[s_name] = json.loads(f["f"])

        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)
        self.salmon_lfc = self.ignore_samples(self.salmon_lfc)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0 and len(self.salmon_lfc) == 0:
            raise ModuleNoSamplesFound

        if len(self.salmon_meta) > 0:
            log.info(f"Found {len(self.salmon_meta)} meta reports")
            self.write_data_file(self.salmon_meta, "multiqc_salmon")
        if len(self.salmon_fld) > 0:
            log.info(f"Found {len(self.salmon_fld)} library insert size distributions")
        if len(self.salmon_lfc) > 0:
            log.info(f"Found {len(self.salmon_lfc)} library format counts reports")

        if self.salmon_meta:
            # Add alignment rate to the general stats table
            # Convert library types to string:
            for d in self.salmon_meta.values():
                if "library_types" in d:
                    d["library_types"] = ", ".join(d["library_types"])

            headers = {
                "percent_mapped": {
                    "title": "% Aligned",
                    "description": "% Mapped reads",
                    "max": 100,
                    "min": 0,
                    "suffix": "%",
                    "scale": "YlGn",
                },
                "num_mapped": {
                    "title": "M Aligned",
                    "description": "Mapped reads (millions)",
                    "min": 0,
                    "scale": "PuRd",
                    "modify": lambda x: float(x) * config.read_count_multiplier,
                    "suffix": config.read_count_prefix,
                    "shared_key": "read_count",
                },
                "library_types": {
                    "title": "Library types",
                    "description": "Library types",
                    "scale": False,
                    # Hide if all samples have the same value
                    "hidden": len(set(d.get("library_types") for d in self.salmon_meta.values())) == 1,
                },
            }
            self.general_stats_addcols(self.salmon_meta, headers)

        if self.salmon_lfc:
            # Compatible fragments ratios data
            lfc_headers = {
                "compatible_fragment_ratio": {
                    "title": "CFR",
                    "description": "Compatible fragment ratio",
                    "min": 0.0,
                    "max": 100.0,
                    "modify": lambda x: x * 100.0,
                    "suffix": "%",
                    "scale": "YlGn",
                },
                "strand_mapping_bias": {
                    "title": "M Bias",
                    "description": "Strand mapping bias",
                    "scale": "BuGn",
                    "max": 1.0,
                },
            }
            # add strand mapping bias data
            self.general_stats_addcols(self.salmon_lfc, lfc_headers)

        if self.salmon_fld:
            salmon_fld_pct = {
                sample: {x: y * 100.0 for x, y in series.items()} for sample, series in self.salmon_fld.items()
            }

            # Library insert size distribution plot
            pconfig = {
                "smooth_points": 500,
                "id": "salmon_plot",
                "title": "Inferred Library Insert Size Distribution",
                "ylab": "Percentage",
                "ysuffix": "%",
                "xlab": "Library Insert Size (bp)",
                "ymin": 0,
                "xmin": 0,
                "tt_label": "<b>{point.x:,.0f} bp</b>: {point.y:.2f}%",
            }
            self.add_section(
                name="Library Insert Size Distribution",
                description="Inferred library insert size density estimated by Salmon.",
                plot=linegraph.plot(salmon_fld_pct, pconfig),
            )
