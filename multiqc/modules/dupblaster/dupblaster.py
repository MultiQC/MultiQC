import logging
from typing import Dict, List, Optional, Tuple, Union

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table
from multiqc.plots.table_object import InputRow
from multiqc.types import LoadedFileDict, SampleGroup, SampleName

log = logging.getLogger(__name__)

Number = Union[int, float]
Metrics = Dict[str, Optional[Number]]

TEXT_COLUMNS = ("sample", "library", "dupblaster_version", "category", "sequencing_unit")

# Duplicate metrics that add up over a sample's libraries.
SUMMED_COLUMNS = (
    "total_templates",
    "duplicate_templates",
    "mapped_pairs",
    "unmapped_pairs",
    "duplicate_pairs",
    "raw_sequencing_duplicate_pairs",
    "corrected_sequencing_duplicate_pairs",
    "library_duplicate_pairs",
    "estimated_library_size",
    "mapped_orphans",
    "duplicate_orphans",
    "unmapped_orphans",
    "unmated_templates",
)

# Occurrence-count bins (label, smallest count) for the duplication spectrum, as in FastQC's duplication levels.
SPECTRUM_BINS: List[Tuple[str, int]] = [
    *((f"{k}x", k) for k in range(1, 10)),
    ("10-49x", 10),
    ("50-99x", 50),
    ("100-499x", 100),
    ("500-999x", 500),
    ("1,000-4,999x", 1_000),
    ("5,000-9,999x", 5_000),
    ("10,000x+", 10_000),
]
SPECTRUM_COLORS = [
    "#313695",
    "#384b9f",
    "#3e60aa",
    "#4575b4",
    "#5588be",
    "#649ac7",
    "#74add1",
    "#a2bebb",
    "#d0cfa6",
    "#fee090",
    "#fbba76",
    "#f7935d",
    "#f46d43",
    "#da4939",
    "#bf2430",
    "#a50026",
]


class MultiqcModule(BaseMultiqcModule):
    """
    dupblaster marks or removes PCR duplicates in query-grouped SAM/BAM, as a streaming successor to samblaster.
    It also splits duplicate pairs into _sequencing_ duplicates, made on the flowcell (optical or ExAmp
    duplicates), and _library_ duplicates, which come from independent molecules (PCR copies, or distinct
    molecules that share a locus). It tells the two apart by imaging-tile identity and corrects for tiles that
    collide by chance, so it needs no pixel distance.

    Every run writes its metrics under `--metrics-prefix <PREFIX>`, and the module reads each file by its header:

    - `<PREFIX>.duplicate-metrics.tsv`: one row per library, giving the duplicate categories plot, the duplicate
      metrics table and the General Statistics columns. The `--stats` file of dupblaster 0.1 and 0.2 has the same
      header and is read too; it has no sequencing duplicate columns.
    - `<PREFIX>.sequencing-units.tsv`: sequencing duplicates per flowcell lane.
    - `<PREFIX>.duplication-sampled.tsv`: the duplicate rate as templates accumulate.
    - `<PREFIX>.duplication-spectrum.tsv`: how many molecules were seen once, twice, and so on. dupblaster writes it
      only with `--duplication-spectrum on`.

    A compatible run:

    ```bash
    samtools sort -n sample.bam \\
      | dupblaster --metrics-prefix sample.dupblaster --duplication-spectrum on -o sample.dups.bam
    ```

    #### Sample names

    Each row is named by its `sample` column (dupblaster's `--sample`, or the read groups' `SM` values), or by the
    file name when that column is empty. A file holding several libraries gives each library its own row, named
    `<sample> (<library>)`, and the duplicate metrics table nests those rows under a row for the sample. To name rows
    by file instead, use:

    ```yaml
    use_filename_as_sample_name:
      - dupblaster
    ```

    The `.dupblaster` stem of the recommended prefix and the four file suffixes are removed from file names.
    """

    def __init__(self):
        super().__init__(
            name="dupblaster",
            anchor="dupblaster",
            href="https://github.com/fulcrumgenomics/dupblaster",
            info="Marks PCR duplicates in query-grouped SAM/BAM and splits them into sequencing and library duplicates.",
            doi="10.5281/zenodo.21445780",
            license="MIT License",
            license_url="https://github.com/fulcrumgenomics/dupblaster/blob/main/LICENSE",
        )

        self.sample_of: Dict[str, str] = {}
        metrics = self.parse_duplicate_metrics()
        units = self.parse_sequencing_units()
        sampled = self.parse_duplication_sampled()
        spectrum = self.parse_duplication_spectrum()
        if not (metrics or units or sampled or spectrum):
            raise ModuleNoSamplesFound
        log.info(
            f"Found {len(metrics)} duplicate metrics, {len(units)} sequencing unit, {len(sampled)} duplication "
            f"sampled and {len(spectrum)} duplication spectrum reports"
        )

        if metrics:
            self.add_general_stats(metrics)
            self.add_categories_section(metrics)
            self.add_metrics_section(metrics)
        if units:
            self.add_sequencing_units_section(units)
        if sampled:
            self.add_duplication_sampled_section(sampled)
        if spectrum:
            self.add_duplication_spectrum_section(spectrum)

        if metrics:
            self.write_data_file(metrics, "multiqc_dupblaster")
        if units:
            self.write_data_file(
                {f"{key} {unit}": d for key, by_unit in units.items() for unit, d in by_unit.items()},
                "multiqc_dupblaster_sequencing_units",
            )
        if sampled:
            self.write_data_file(
                {
                    key: {f"{label}_{total}": pct for label, line in d.items() for total, pct in line.items()}
                    for key, d in sampled.items()
                },
                "multiqc_dupblaster_duplication_sampled",
            )
        if spectrum:
            self.write_data_file(
                {key: {str(k): n for k, n in d.items()} for key, d in spectrum.items()},
                "multiqc_dupblaster_duplication_spectrum",
            )

    def row_names(self, f: LoadedFileDict[str], rows: List[Dict[str, str]]) -> Dict[Tuple[str, str], str]:
        """The report name of each (sample, library) in a file: the sample, plus the library if the file has several."""
        pairs = list(dict.fromkeys((row["sample"], row["library"]) for row in rows))
        names: Dict[Tuple[str, str], str] = {}
        for sample, library in pairs:
            s_name = self.clean_s_name(sample, f) if sample else f["s_name"]
            name = f"{s_name} ({library})" if len(pairs) > 1 else s_name
            names[(sample, library)] = name
            self.sample_of[name] = s_name
        return names

    def parse_duplicate_metrics(self) -> Dict[str, Metrics]:
        data: Dict[str, Metrics] = {}
        for f in self.find_log_files("dupblaster/duplicate_metrics"):
            rows = read_rows(f)
            names = self.row_names(f, rows)
            for row in rows:
                s_name = names[(row["sample"], row["library"])]
                if s_name in data:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                d = {col: to_number(val) for col, val in row.items() if col not in TEXT_COLUMNS}
                # dupblaster writes a rate of 0, not an empty cell, when the rate has no denominator.
                if required(d, "mapped_pairs") == 0:
                    d["frac_duplicate_pairs"] = None
                if required(d, "mapped_pairs") + required(d, "mapped_orphans") == 0:
                    d["frac_duplicates"] = None
                # Absent before dupblaster 0.3.0, and empty when the sequencing duplicate split was not computed.
                sequencing = d.get("corrected_sequencing_duplicate_pairs")
                duplicate_pairs = required(d, "duplicate_pairs")
                if sequencing is not None and duplicate_pairs > 0:
                    d["frac_duplicate_pairs_sequencing"] = sequencing / duplicate_pairs
                data[s_name] = d
                self.add_data_source(f, s_name, section="duplicate_metrics")
                self.add_software_version(row["dupblaster_version"], s_name)
        return self.ignore_samples(data)

    def parse_sequencing_units(self) -> Dict[str, Dict[str, Metrics]]:
        data: Dict[str, Dict[str, Metrics]] = {}
        for f in self.find_log_files("dupblaster/sequencing_units"):
            rows = read_rows(f)
            names = self.row_names(f, rows)
            for row in rows:
                # A library with no pairs gets one placeholder row with no sequencing unit.
                if row["sequencing_unit"] == "":
                    continue
                s_name = names[(row["sample"], row["library"])]
                data.setdefault(s_name, {})[row["sequencing_unit"]] = {
                    col: to_number(val) for col, val in row.items() if col not in TEXT_COLUMNS
                }
                self.add_data_source(f, s_name, section="sequencing_units")
        return self.ignore_samples(data)

    def parse_duplication_sampled(self) -> Dict[str, Dict[str, Dict[int, float]]]:
        data: Dict[str, Dict[str, Dict[int, float]]] = {}
        for f in self.find_log_files("dupblaster/duplication_sampled"):
            rows = read_rows(f)
            names = self.row_names(f, rows)
            for row in rows:
                s_name = names[(row["sample"], row["library"])]
                d = data.setdefault(s_name, {"window": {}, "cumulative": {}})
                total = int(row["total"])
                d["window"][total] = float(row["window_frac_duplicates"]) * 100.0
                d["cumulative"][total] = float(row["frac_duplicates"]) * 100.0
                self.add_data_source(f, s_name, section="duplication_sampled")
        return self.ignore_samples(data)

    def parse_duplication_spectrum(self) -> Dict[str, Dict[int, int]]:
        data: Dict[str, Dict[int, int]] = {}
        for f in self.find_log_files("dupblaster/duplication_spectrum"):
            rows = read_rows(f)
            names = self.row_names(f, rows)
            for row in rows:
                s_name = names[(row["sample"], row["library"])]
                data.setdefault(s_name, {})[int(row["n_observations"])] = int(row["n_molecules"])
                self.add_data_source(f, s_name, section="duplication_spectrum")
        return self.ignore_samples(data)

    def add_general_stats(self, metrics: Dict[str, Metrics]) -> None:
        headers = {
            "frac_duplicates": {
                "title": "% Dups",
                "description": "Duplicate reads as a percentage of mapped reads (Picard's PERCENT_DUPLICATION)",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "modify": lambda x: x * 100.0,
            },
            "frac_sequencing_duplicate_pairs": {
                "title": "% Seq Dups",
                "description": "Mapped pairs that are sequencing (flowcell) duplicates, corrected for tiles that "
                "collide by chance",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "modify": lambda x: x * 100.0,
                "hidden": True,
            },
            "estimated_library_size": {
                "title": "Library Size",
                "description": "Estimated distinct molecules in the library, with sequencing duplicates excluded",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "YlGn",
                "hidden": True,
            },
        }
        general_stats_headers = self.get_general_stats_headers(all_headers=headers)
        if general_stats_headers:
            self.general_stats_addcols(
                {s_name: {k: v for k, v in d.items() if v is not None} for s_name, d in metrics.items()},
                general_stats_headers,
            )

    def add_categories_section(self, metrics: Dict[str, Metrics]) -> None:
        plot_data: Dict[str, Dict[str, Number]] = {}
        for s_name, d in metrics.items():
            duplicate_pairs = required(d, "duplicate_pairs")
            duplicate_orphans = required(d, "duplicate_orphans")
            categories = {
                "unique_pairs": required(d, "mapped_pairs") - duplicate_pairs,
                "unique_unpaired": required(d, "mapped_orphans") - duplicate_orphans,
                "duplicate_unpaired": duplicate_orphans,
                "unmapped": required(d, "unmapped_pairs") + required(d, "unmapped_orphans"),
            }
            sequencing = d.get("corrected_sequencing_duplicate_pairs")
            library = d.get("library_duplicate_pairs")
            if sequencing is not None and library is not None:
                categories["library_duplicate_pairs"] = library
                categories["sequencing_duplicate_pairs"] = sequencing
            else:
                categories["unsplit_duplicate_pairs"] = duplicate_pairs
            plot_data[s_name] = categories

        cats = {
            "unique_pairs": {"name": "Unique pairs", "color": "#437bb1"},
            "unique_unpaired": {"name": "Unique unpaired", "color": "#7cb5ec"},
            "library_duplicate_pairs": {"name": "Library duplicate pairs", "color": "#b1084c"},
            "sequencing_duplicate_pairs": {"name": "Sequencing duplicate pairs", "color": "#f7a35c"},
            "unsplit_duplicate_pairs": {"name": "Duplicate pairs (not split)", "color": "#7f0000"},
            "duplicate_unpaired": {"name": "Duplicate unpaired", "color": "#e4a2c0"},
            "unmapped": {"name": "Unmapped", "color": "#999999"},
        }
        self.add_section(
            name="Duplicate Categories",
            anchor=f"{self.anchor}-categories",
            description="Templates of each library by duplicate state. Each template counts once, so a library's "
            "bars add up to the templates dupblaster examined.",
            helptext="""
            A template is one read pair, or one read for unpaired data.

            * _Unique pairs_ and _Unique unpaired_: templates kept as the representative of their duplicate set.
            * _Library duplicate pairs_: copies of a molecule that existed before sequencing, from PCR, or distinct
              molecules that share a locus. Many of these point to an over-amplified or low-complexity library.
            * _Sequencing duplicate pairs_: copies made on the flowcell (optical or ExAmp duplicates), corrected for
              tiles that collide by chance. Many of these point to an overloaded flowcell, not to the library.
            * _Duplicate pairs (not split)_: duplicate pairs of a library whose split was not computed: the run used
              `--sequencing-duplicate-detection off`, all its pairs sit on one tile, or dupblaster predates 0.3.0.
            * _Unpaired_: templates with one mapped read, either a single-end read or a pair whose mate is unmapped
              or missing.
            * _Unmapped_: templates with no mapped read. They are never checked for duplicates.
            """,
            plot=bargraph.plot(
                plot_data,
                cats,
                {
                    "id": f"{self.anchor}_categories",
                    "title": "dupblaster: Duplicate Categories",
                    "ylab": "Templates",
                    "cpswitch_counts_label": "Number of Templates",
                    "cpswitch_c_active": False,
                },
            ),
        )

    def add_metrics_section(self, metrics: Dict[str, Metrics]) -> None:
        headers = {
            "total_templates": {
                "title": f"{config.read_count_prefix} Templates",
                "description": f"Templates examined ({config.read_count_desc})",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "mapped_pairs": {
                "title": f"{config.read_count_prefix} Mapped Pairs",
                "description": f"Templates with both reads mapped ({config.read_count_desc})",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
            "frac_duplicates": {
                "title": "% Dups",
                "description": "Duplicate reads as a percentage of mapped reads (Picard's PERCENT_DUPLICATION)",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "modify": lambda x: x * 100.0,
            },
            "frac_duplicate_pairs": {
                "title": "% Dup Pairs",
                "description": "Duplicate pairs as a percentage of mapped pairs",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "modify": lambda x: x * 100.0,
            },
            "frac_sequencing_duplicate_pairs": {
                "title": "% Seq Dup Pairs",
                "description": "Sequencing (flowcell) duplicate pairs as a percentage of mapped pairs, corrected for "
                "tiles that collide by chance; part of % Dup Pairs",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "modify": lambda x: x * 100.0,
            },
            "frac_duplicate_pairs_sequencing": {
                "title": "Seq Share",
                "description": "Share of duplicate pairs that are sequencing (flowcell) duplicates; the rest are "
                "library duplicates",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "PuOr",
                "modify": lambda x: x * 100.0,
            },
            "estimated_library_size": {
                "title": "Library Size",
                "description": "Lander-Waterman estimate of distinct molecules, with sequencing duplicates excluded",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "YlGn",
            },
            "mapped_orphans": {
                "title": "Unpaired",
                "description": "Templates with exactly one mapped read",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Blues",
                "hidden": True,
            },
            "unmated_templates": {
                "title": "Unmated",
                "description": "Paired templates whose mate was missing (only counted under --ignore-unmated)",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Reds",
                "hidden": True,
            },
        }
        rows_by_group: Dict[Union[str, SampleGroup], List[InputRow]] = {}
        for s_name, d in metrics.items():
            rows_by_group.setdefault(SampleGroup(self.sample_of[s_name]), []).append(
                InputRow(sample=SampleName(s_name), data=d)
            )
        for sample, rows in rows_by_group.items():
            if len(rows) > 1:
                summary = sum_libraries([metrics[row.sample] for row in rows])
                rows.insert(0, InputRow(sample=SampleName(sample), data=summary))

        self.add_section(
            name="Duplicate Metrics",
            anchor=f"{self.anchor}-metrics",
            description="Duplicate rates and library size for each library.",
            helptext="""
            * _% Dups_ counts reads, like Picard's `PERCENT_DUPLICATION`: each duplicate pair adds two reads.
            * _% Dup Pairs_ and _% Seq Dup Pairs_ share the mapped-pair denominator, so the sequencing rate is the
              part of the pair rate that the flowcell made.
            * _Seq Share_ is the sequencing duplicates' share of all duplicate pairs. A high share says to load the
              flowcell less densely; a low share with a high duplicate rate says the library is over-amplified.
            * _Library Size_ is the Lander-Waterman estimate of distinct molecules. Sequencing duplicates are left
              out of the observed total, as Picard does, because they are no evidence that the library is exhausted.

            A sample with several libraries has a summary row; click it to show each library. Its counts are the sums
            over its libraries and its rates are recomputed from those sums. Its library size is the sum of the
            libraries' estimates, since each library is its own pool of molecules.
            """,
            plot=table.plot(
                rows_by_group,
                headers,
                {
                    "id": f"{self.anchor}_metrics_table",
                    "title": "dupblaster: Duplicate Metrics",
                    "namespace": "dupblaster",
                },
            ),
        )

    def add_sequencing_units_section(self, units: Dict[str, Dict[str, Metrics]]) -> None:
        rows_by_group: Dict[Union[str, SampleGroup], List[InputRow]] = {}
        for s_name, by_unit in units.items():
            unit_rows = [InputRow(sample=SampleName(f"{s_name} {unit}"), data=d) for unit, d in by_unit.items()]
            if len(unit_rows) == 1:
                rows_by_group[SampleGroup(s_name)] = unit_rows
                continue
            templates = sum(required(d, "templates") for d in by_unit.values())
            duplicates = sum(required(d, "sequencing_duplicate_pairs") for d in by_unit.values())
            total = {
                "units": len(by_unit),
                "templates": templates,
                "tiles": sum(required(d, "tiles") for d in by_unit.values()),
                "sequencing_duplicate_pairs": duplicates,
                "frac_sequencing_duplicate_pairs": duplicates / templates if templates > 0 else None,
            }
            rows_by_group[SampleGroup(s_name)] = [InputRow(sample=SampleName(s_name), data=total), *unit_rows]

        headers = {
            "units": {
                "title": "Lanes",
                "description": "Flowcell lanes the library was sequenced on",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Greys",
            },
            "templates": {
                "title": f"{config.read_count_prefix} Mapped Pairs",
                "description": f"Mapped pairs from this flowcell lane ({config.read_count_desc})",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "tiles": {
                "title": "Tiles",
                "description": "Distinct imaging tiles seen on this flowcell lane",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Greys",
            },
            "sequencing_duplicate_pairs": {
                "title": "Seq Dup Pairs",
                "description": "Sequencing duplicate pairs on this lane's tiles, before the chance-collision correction",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "OrRd",
                "hidden": True,
            },
            "frac_sequencing_duplicate_pairs": {
                "title": "% Seq Dup Pairs",
                "description": "Sequencing duplicate pairs as a percentage of this lane's mapped pairs, before the "
                "chance-collision correction",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "modify": lambda x: x * 100.0,
            },
        }
        self.add_section(
            name="Sequencing Duplicates by Lane",
            anchor=f"{self.anchor}-sequencing-units",
            description="Sequencing (flowcell) duplicates of each library on each flowcell lane.",
            helptext="""
            The sequencing duplicate rate can differ a lot between the flowcells and lanes of one library, which a
            per-library average hides. Each row is named by the library and its flowcell and lane, as dupblaster read
            them from the read names. A library sequenced on several lanes has a summary row; click it to show each
            lane.

            These counts are not corrected for tiles that collide by chance, so they are a little higher than the
            corrected figures in the duplicate metrics table. They add up to the library's uncorrected
            sequencing duplicate count.
            """,
            plot=table.plot(
                rows_by_group,
                headers,
                {
                    "id": f"{self.anchor}_sequencing_units_table",
                    "title": "dupblaster: Sequencing Duplicates by Lane",
                    "namespace": "dupblaster",
                },
            ),
        )

    def add_duplication_sampled_section(self, sampled: Dict[str, Dict[str, Dict[int, float]]]) -> None:
        self.add_section(
            name="Duplication Rate by Depth",
            anchor=f"{self.anchor}-duplication-sampled",
            description="Duplicate rate of each library as templates accumulate, in the order dupblaster read them.",
            helptext="""
            dupblaster takes a snapshot every `--sampling-interval` templates (1,000,000 by default) and one at
            the end. _Per window_ is the duplicate rate of the templates added since the previous snapshot; a curve
            that keeps rising says deeper sequencing would mostly add duplicates. _Cumulative_ is the duplicate rate
            of all templates so far. A library with any mapped pairs is measured on its pairs, and an unpaired
            library on its reads.

            The curve depends on input order. On unsorted aligner output from one lane it reads as a saturation
            curve. On sorted or multi-flowcell input, flowcell duplicates land next to each other, so the shape is
            a fingerprint of the flowcells rather than of library complexity; use the duplication spectrum for that.
            """,
            plot=linegraph.plot(
                [
                    {s_name: d["window"] for s_name, d in sampled.items()},
                    {s_name: d["cumulative"] for s_name, d in sampled.items()},
                ],
                {
                    "id": f"{self.anchor}_duplication_sampled",
                    "title": "dupblaster: Duplication Rate by Depth",
                    "xlab": "Templates",
                    "ylab": "% Duplicates",
                    "ysuffix": "%",
                    "ymin": 0,
                    "data_labels": [
                        {"name": "Per window", "ylab": "% Duplicates in window"},
                        {"name": "Cumulative", "ylab": "% Duplicates"},
                    ],
                },
            ),
        )

    def add_duplication_spectrum_section(self, spectrum: Dict[str, Dict[int, int]]) -> None:
        molecules: Dict[str, Dict[str, int]] = {}
        templates: Dict[str, Dict[str, int]] = {}
        for s_name, histogram in spectrum.items():
            molecules[s_name] = {label: 0 for label, _ in SPECTRUM_BINS}
            templates[s_name] = {label: 0 for label, _ in SPECTRUM_BINS}
            for k, n in histogram.items():
                label = next(label for label, smallest in reversed(SPECTRUM_BINS) if k >= smallest)
                molecules[s_name][label] += n
                templates[s_name][label] += k * n

        cats = {label: {"name": label, "color": color} for (label, _), color in zip(SPECTRUM_BINS, SPECTRUM_COLORS)}
        self.add_section(
            name="Duplication Spectrum",
            anchor=f"{self.anchor}-duplication-spectrum",
            description="How many times each distinct molecule was seen, as a share of molecules and of templates.",
            helptext="""
            dupblaster counts how many times it saw each distinct molecule, and the bars bin those counts as
            FastQC does for duplicated sequences. _Molecules_ shows the share of distinct molecules seen once,
            twice, and so on. _Templates_ weights each molecule by its count, so it shows where the sequenced
            templates went. Templates riding well above molecules in the high bins point to a few heavily copied
            molecules, from PCR or flowcell jackpots, carrying a large share of the data. Unlike the duplication rate
            by depth, the spectrum does not depend on input order.

            dupblaster writes the spectrum only with `--duplication-spectrum on`. Counts above 65,535 are reported
            as 65,535.
            """,
            plot=bargraph.plot(
                [molecules, templates],
                [cats, cats],
                {
                    "id": f"{self.anchor}_duplication_spectrum",
                    "title": "dupblaster: Duplication Spectrum",
                    "ylab": "Molecules",
                    "cpswitch_c_active": False,
                    "data_labels": [
                        {"name": "Molecules", "ylab": "Molecules"},
                        {"name": "Templates", "ylab": "Templates"},
                    ],
                },
            ),
        )


def read_rows(f: LoadedFileDict[str]) -> List[Dict[str, str]]:
    """The rows of a dupblaster TSV as `{column: value}`; a row with the wrong number of fields is a hard error."""
    lines = f["f"].splitlines()
    header = lines[0].split("\t")
    rows = []
    for line_num, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"dupblaster: {f['fn']} line {line_num} has {len(fields)} fields, expected {len(header)}")
        rows.append(dict(zip(header, fields)))
    return rows


def to_number(value: str) -> Optional[Number]:
    """A metric cell as a number, or None for the empty cell dupblaster writes when it cannot compute a metric."""
    if value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return float(value)


def sum_libraries(libraries: List[Metrics]) -> Metrics:
    """A sample's duplicate metrics from its libraries: counts summed, rates recomputed as dupblaster defines them."""
    d: Metrics = {}
    for column in SUMMED_COLUMNS:
        values = [library.get(column) for library in libraries]
        d[column] = None if any(v is None for v in values) else sum(v for v in values if v is not None)
    mapped_pairs = required(d, "mapped_pairs")
    duplicate_pairs = required(d, "duplicate_pairs")
    mapped_reads = required(d, "mapped_orphans") + 2 * mapped_pairs
    duplicate_reads = required(d, "duplicate_orphans") + 2 * duplicate_pairs
    sequencing = d["corrected_sequencing_duplicate_pairs"]
    d["frac_duplicates"] = duplicate_reads / mapped_reads if mapped_reads > 0 else None
    d["frac_duplicate_pairs"] = duplicate_pairs / mapped_pairs if mapped_pairs > 0 else None
    d["frac_sequencing_duplicate_pairs"] = (
        sequencing / mapped_pairs if sequencing is not None and mapped_pairs > 0 else None
    )
    d["frac_duplicate_pairs_sequencing"] = (
        sequencing / duplicate_pairs if sequencing is not None and duplicate_pairs > 0 else None
    )
    return d


def required(d: Metrics, column: str) -> Number:
    """A metric dupblaster always writes; an empty cell means the file is not what the module expects."""
    value = d[column]
    if value is None:
        raise ValueError(f"dupblaster: required column '{column}' is empty")
    return value
