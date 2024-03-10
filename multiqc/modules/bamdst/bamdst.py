""" MultiQC module to parse output from Bamdst """


import logging
from collections import defaultdict
from typing import Dict, Union
import csv
import math
import os
import fnmatch

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table, bargraph, linegraph
from multiqc.utils import config

# Initialise the logger
log = logging.getLogger(__name__)


default_exclude_contigs = [
    "*_alt",
    "*_decoy",
    "*_random",
    "*_fix",
    "HLA*",
    "chrUn*",
    "chrEBV",
    "chrM",
]


def _read_config():
    cfg = getattr(config, "bamdst", dict())
    if not isinstance(cfg, dict):
        return {}

    cfg["exclude_include_are_set_by_user"] = "exclude_contigs" in cfg or "include_contigs" in cfg

    cfg["include_contigs"] = cfg.get("include_contigs", [])
    if not isinstance(cfg["include_contigs"], list):
        cfg["include_contigs"] = []

    cfg["exclude_contigs"] = cfg.get("exclude_contigs", default_exclude_contigs)
    if not isinstance(cfg["exclude_contigs"], list):
        cfg["exclude_contigs"] = []

    if cfg["include_contigs"]:
        log.debug(f"Configuration: would try to include only these contigs: {', '.join(cfg['include_contigs'])}")
    if cfg["exclude_contigs"]:
        log.debug(
            f"Configuration: would exclude these contigs, unless all would be filtered: {', '.join(cfg['exclude_contigs'])}"
        )

    cutoff = cfg.get("perchrom_fraction_cutoff", 0.0)
    try:
        cutoff = float(cutoff)
    except ValueError:
        cutoff = 0.0
    if cutoff != 0.0:
        log.debug(
            f"Configuration: coverage cutoff of contigs is {cutoff * 100.0}%, unless everything would be filtered"
        )
    cfg["perchrom_fraction_cutoff"] = cutoff

    return cfg


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Bamdst",
            anchor="bamdst",
            href="https://https://github.com/shiquan/bamdst",
            info="""is a lightweight tool to stat the depth coverage of target regions of bam file(s)""",
            # doi="", # No DOI
        )

        self.cfg = _read_config()
        data_by_sample = dict()
        data_by_chromosome_by_sample = dict()
        for f in self.find_log_files("bamdst/coverage"):
            if data := self._parse_coverage_report(f):
                data_by_sample[f["s_name"]] = data

                # Parse "chromosomes.report" sitting next to the "coverage.report" file:
                if os.path.exists(chroms_path := os.path.join(f["root"], "chromosomes.report")):
                    if data := self._parse_chromosomes_report(chroms_path, f["s_name"]):
                        data_by_chromosome_by_sample[f["s_name"]] = data

        data_by_sample = self.ignore_samples(data_by_sample)
        data_by_chromosome_by_sample = self.ignore_samples(data_by_chromosome_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} reports")
        self.write_data_file(data_by_sample, "multiqc_bamdst_coverage")
        self.write_data_file(data_by_chromosome_by_sample, "multiqc_bamdst_chromosomes")

        self.build_tables(data_by_sample)
        self._build_per_chrom_plot(data_by_chromosome_by_sample)

    def _parse_coverage_report(self, f: Dict) -> Dict:
        """
        Parse one coverage report.
        """
        data = dict()
        s_name = None
        version = None
        for line in f["f"].split("\n"):
            line = line.strip()
            if line.startswith("## Files : "):
                # Split a line like "## Files : example/test1.bam example/test2.bam "
                # Cannot split by spaces because the file names can contain spaces as well
                path_line = line.split(":")[1].strip()
                for ext in [".bam", ".cram", ".sam"]:
                    path_line = path_line.replace(ext, "<EXT>")
                names = [
                    self.clean_s_name(sn.strip())
                    for sn in path_line.split("<EXT>")
                    if sn
                    if sn.strip() not in ["-", "/dev/stdin"]
                ]
                if names:
                    f["s_name"] = "-".join(names)  # concatenating file names
                continue
            if line.startswith("## Version : "):
                version = line.split(":")[1].strip()
                continue
            fields = line.split("\t")
            if not len(fields) == 2:
                continue
            metric, value = fields
            if value.endswith("%"):
                value = value.strip("%")
                try:
                    value = float(value)
                except ValueError:
                    continue
                value = value / 100.0
            else:
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        continue
            data[metric] = value

        if s_name and version and data:
            self.add_software_version(version, sample=s_name)
            self.add_data_source(f, s_name=s_name, section="coverage table")
        return data

    def _parse_chromosomes_report(self, path: str, s_name: str) -> Dict[str, Dict[str, Union[float, str]]]:
        data_by_contig = defaultdict(dict)
        with open(path) as fh:
            reader: csv.DictReader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                contig = row["#Chromosome"].strip()
                for k, v in row.items():
                    if k == "#Chromosome":
                        continue
                    v = v.strip()
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                    data_by_contig[contig][k.strip()] = v

        if data_by_contig:
            self.add_data_source(source=path, s_name=s_name, section="per-chromosome")
        return data_by_contig

    def build_tables(self, data_by_sample):
        headers = {
            "[Total] Mapped Reads": {
                "title": f"Mapped ({config.read_count_prefix})",
                "description": f"Mapped reads ({config.read_count_desc})",
                "scale": "RdYlGn",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
            "[Total] Fraction of Mapped Reads": {
                "title": "Mapped",
                "description": f"Fraction of mapped reads in all reads ({config.read_count_desc})",
                "scale": "PuBu",
                "min": 0.0,
                "max": 100.0,
                "suffix": "%",
                "modify": lambda x: x * 100.0,
            },
            "[Target] Target Reads": {
                "title": f"Target reads ({config.read_count_prefix})",
                "description": f"Reads mapped on target ({config.read_count_desc})",
                "scale": "RdYlGn",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "[Target] Fraction of Target Reads in all reads": {
                "title": "% in all",
                "description": "Fraction of target reads in all reads",
                "scale": "PuBu",
                "min": 0.0,
                "max": 100.0,
                "suffix": "%",
                "modify": lambda x: x * 100.0,
            },
            "[Target] Fraction of Target Reads in mapped reads": {
                "title": "% in mapped",
                "description": "Fraction of target reads in mapped reads",
                "scale": "PuBu",
                "min": 0.0,
                "max": 100.0,
                "suffix": "%",
                "modify": lambda x: x * 100.0,
            },
            "[Target] Average depth": {
                "title": "Avg depth",
                "description": "Average depth, on target",
                "min": 0,
                "suffix": "X",
                "scale": "BuPu",
            },
            "[Target] Len of region": {
                "title": "Region len",
                "description": "Length of target region",
                "min": 0,
                "format": "{:,d}",
                "suffix": " bp",
                "scale": "Greys",
            },
        }

        coverage_metrics = []
        for s_name, d in data_by_sample.items():
            for k in d.keys():
                if k.startswith("[Target] Coverage (") and k not in coverage_metrics:
                    coverage_metrics.append(k)
        for m in coverage_metrics:
            short = m.replace("[Target] Coverage (", "").replace(")", "")
            headers[m] = {
                "title": short,
                "description": f"Target coverage ({short})",
                "min": 0,
                "suffix": "%",
                "format": "{:.2f}",
                "scale": "RdYlGn",
                "modify": lambda x: x * 100.0,
            }

        self.add_section(
            name="Coverage metrics",
            anchor="bamdst-coverage",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig={
                    "id": "bamdst-coverage-table",
                    "title": "Bamdst: coverage metrics",
                },
            ),
        )

        genstats_headers = {
            k: v
            for k, v in headers.items()
            if k
            in [
                "[Target] Average depth",
                "[Target] Coverage (>0x)",
            ]
        }
        self.general_stats_addcols(data_by_sample, genstats_headers)

    def _filter_contigs(self, data_by_chrom_by_sample):
        # Filter contigs based on include/exclude patterns
        filtered_data_by_chrom_by_sample = defaultdict(dict)
        for s_name, data_by_chrom in data_by_chrom_by_sample.items():
            for contig, data in data_by_chrom.items():
                if any(fnmatch.fnmatch(contig, str(pattern)) for pattern in self.cfg["exclude_contigs"]):
                    try:
                        if self.cfg.get("show_excluded_debug_logs") is True:
                            log.debug(f"Skipping excluded contig '{contig}'")
                    except (AttributeError, KeyError):
                        pass
                    continue

                # filter out contigs based on inclusion patterns
                if len(self.cfg["include_contigs"]) > 0 and not any(
                    fnmatch.fnmatch(contig, pattern) for pattern in self.cfg["include_contigs"]
                ):
                    # Commented out since this could be many thousands of contigs!
                    # log.debug(f"Skipping not included contig '{contig}'")
                    continue
                filtered_data_by_chrom_by_sample[s_name][contig] = data
        if not filtered_data_by_chrom_by_sample:
            filtered_data_by_chrom_by_sample = data_by_chrom_by_sample
            (log.warning if self.cfg.get("exclude_include_are_set_by_user", False) else log.debug)(
                "All contigs would be filtered out by the contigs inclusion/exclusion patterns. "
                "Consider changing the contig.include_contigs/contig.exclude_contigs parameters. "
                "Keeping all contigs."
            )
        data_by_chrom_by_sample = filtered_data_by_chrom_by_sample

        # Filter contigs based on coverage cutoff. Contigs only will be excluded if they
        # fail the cutoff in _all_ samples.
        filtered_data_by_chrom_by_sample = defaultdict(dict)
        if (cutoff := self.cfg.get("perchrom_fraction_cutoff", 0)) > 0:
            passing_contigs = set()
            for s_name, data_by_chrom in data_by_chrom_by_sample.items():
                for contig, data in data_by_chrom.items():
                    pct = data.get("DATA(%)")
                    if pct is not None and isinstance(pct, float) and pct > cutoff:
                        passing_contigs.add(contig)

            rejected_contigs = set()
            for s_name, perchrom in data_by_chrom_by_sample.items():
                for contig, data in perchrom.items():
                    if contig not in passing_contigs:
                        rejected_contigs.add(contig)
                    else:
                        filtered_data_by_chrom_by_sample[s_name][contig] = data

            if rejected_contigs:
                if self.cfg.get("show_excluded_debug_logs") is True:
                    log.debug(
                        f"Skipping {len(rejected_contigs)} contigs not passing the "
                        f"cutoff of {self.cfg['perchrom_fraction_cutoff']}%:. "
                        f"Skipped contigs: {''.join(rejected_contigs)}"
                    )
            if not filtered_data_by_chrom_by_sample:
                filtered_data_by_chrom_by_sample = data_by_chrom_by_sample
                log.warning(
                    f"All contigs would be filtered out by the cutoff of "
                    f"{self.cfg['perchrom_fraction_cutoff']}%. Consider changing the cut-off. "
                    f"Keeping all contigs."
                )
            data_by_chrom_by_sample = filtered_data_by_chrom_by_sample

        return data_by_chrom_by_sample

    def _build_per_chrom_plot(self, data_by_chrom_by_sample):
        def _chrom_key(name):
            k = "".join(c for c in name if c.isdigit())
            return int(k) if k else math.inf

        # Sort chromosomes by name. Sort smartly: chr2 should go before chr10
        for s_name, data_by_chrom in data_by_chrom_by_sample.items():
            data_by_chrom_by_sample[s_name] = dict(sorted(data_by_chrom.items(), key=lambda x: _chrom_key(x[0])))

        # Filter contigs
        filt_data_by_chrom_by_sample = self._filter_contigs(data_by_chrom_by_sample)
        contigs = set.union(*[set(d.keys()) for d in data_by_chrom_by_sample.values()])
        filt_contigs = set.union(*[set(d.keys()) for d in filt_data_by_chrom_by_sample.values()])
        if contigs == filt_contigs or len(filt_contigs) == 0:
            datasets = [data_by_chrom_by_sample]
            data_labels = None
        else:
            datasets = [filt_data_by_chrom_by_sample, data_by_chrom_by_sample]
            data_labels = ["Main contigs", "All contigs"]

        depth_datasets = []
        cov_datasets = []
        for dataset in datasets:
            # Extracting two metrics to plot: average depth and coverage percentage
            depth_datasets.append(defaultdict(dict))
            cov_datasets.append(defaultdict(dict))
            for s_name, data_by_chrom in dataset.items():
                for chrom, d in data_by_chrom.items():
                    depth_datasets[-1][s_name][chrom] = d["Avg depth"]
                    cov_datasets[-1][s_name][chrom] = d["Coverage%"]

        num_chroms = max([len(by_chrom.keys()) for by_chrom in depth_datasets[0].values()])
        if num_chroms > 1:
            perchrom_depth_plot = linegraph.plot(
                depth_datasets,
                pconfig={
                    "id": "bamdst-depth-per-contig-plot",
                    "title": "Bamdst: average depth per contig",
                    "xlab": "Region",
                    "ylab": "Average depth",
                    "categories": True,
                    "tt_decimals": 1,
                    "tt_suffix": "x",
                    "smooth_points": 500,
                    "logswitch": True,
                    "hide_zero_cats": False,
                    "ymin": 0,
                    "data_labels": data_labels,
                },
            )
            perchrom_cov_plot = linegraph.plot(
                cov_datasets,
                pconfig={
                    "id": "bamdst-cov-per-contig-plot",
                    "title": "Bamdst: coverage percentage of each contig",
                    "xlab": "Region",
                    "ylab": "Coverage %",
                    "categories": True,
                    "tt_decimals": 1,
                    "tt_suffix": "%",
                    "smooth_points": 500,
                    "logswitch": True,
                    "hide_zero_cats": False,
                    "ymax": 100,
                    "ymin": 0,
                    "data_labels": data_labels,
                },
            )
        else:
            perchrom_depth_plot = bargraph.plot(
                depth_datasets,
                pconfig={
                    "id": "bamdst-depth-per-contig-plot",
                    "title": "Bamdst: average depth",
                    "xlab": "Sample",
                    "ylab": "Average depth",
                    "tt_suffix": "x",
                    "hide_zero_cats": False,
                    "ymin": 0,
                    "data_labels": data_labels,
                },
            )
            perchrom_cov_plot = bargraph.plot(
                cov_datasets,
                pconfig={
                    "id": "bamdst-cov-per-contig-plot",
                    "title": "Bamdst: coverage percentage",
                    "xlab": "Sample",
                    "ylab": "Coverage %",
                    "tt_suffix": "%",
                    "hide_zero_cats": False,
                    "ymax": 100,
                    "ymin": 0,
                    "data_labels": data_labels,
                },
            )

        self.add_section(
            name="Average depth per chromosome",
            anchor="bamdst-depth-per-contig",
            plot=perchrom_depth_plot,
        )
        self.add_section(
            name="Coverage % of each contig",
            anchor="bamdst-coverage-per-contig",
            plot=perchrom_cov_plot,
        )
