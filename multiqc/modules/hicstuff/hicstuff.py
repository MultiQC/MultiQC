import ast
import logging
import math
import re
from typing import Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)

# Colour palette (tab10) used to assign one colour per sample in the distance-law plot
_SAMPLE_COLOURS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses two file types from the
    [hicstuff](https://github.com/koszullab/hicstuff) Hi-C pipeline:

    - **Pipeline log files** (`*.log`, `*.txt`), identified by the
      `## hicstuff:` header line. The end-of-run summary stats dictionary
      feeds the General Statistics table and the Read Fate stacked bar plot.
    - **Distance law tables** (default `distance_law.txt`, also commonly
      seen with `.tsv` extensions), identified by a `## distance_law`
      header. These drive the P(s) contact-probability line graph and its
      log-log slope plot.

    hicstuff writes a summary stats dictionary to the log at the end of each run:

    ```
    ## hicstuff: v3.2.2 log file
    ## date: 2024-02-16 14:00:23
    ## enzyme: DpnII,HinfI
    ## input1: ../tinyMapper/tests/testHiC_R1.fq.gz
    ## input2: ../tinyMapper/tests/testHiC_R2.fq.gz
    ## ref: /home/rsg/genomes/S288c/S288c.fa
    ---
    ...
    2024-02-16,14:00:43 :: INFO :: 77% reads (single ends) mapped with Q >= 30 (154272/200000)
    2024-02-16,14:00:44 :: INFO :: 66943 pairs successfully mapped (66.94%)
    2024-02-16,14:00:46 :: INFO :: Fetching mapping and pairing stats
    2024-02-16,14:00:46 :: INFO :: {'Sample': 'testHiC^CGNT57', 'Total read pairs': 100000, 'Mapped reads': 154272, 'Unmapped reads': 45728, 'Recovered contacts': 66943, 'Final contacts': 66943, 'Removed contacts': 0, 'Filtered out': 0, 'Loops': 0, 'Uncuts': 0, 'Weirds': 0, 'PCR duplicates': 0}
    2024-02-16,14:00:46 :: INFO :: Contact map generated after 0h 0m 23s
    ```
    """

    def __init__(self):
        super().__init__(
            name="hicstuff",
            anchor="hicstuff",
            href="https://github.com/koszullab/hicstuff",
            info="Hi-C pipeline that generates contact maps from sequencing reads.",
            doi="10.5281/zenodo.4066363",
        )

        self.hicstuff_data = {}

        for f in self.find_log_files("hicstuff/pipeline_stats"):
            parsed = self.parse_hicstuff_log(f)
            if parsed is not None:
                s_name = self.clean_s_name(parsed["Sample"], f)
                self.hicstuff_data[s_name] = parsed
                self.add_data_source(f, s_name)
                if "version" in parsed:
                    self.add_software_version(str(parsed["version"]), sample=s_name)

        self.add_software_version(None)

        self.hicstuff_data = self.ignore_samples(self.hicstuff_data)

        # Distance-law files
        # Structure: {s_name: {chrom: {start_bp: p_s, ...}}}
        self.hicstuff_distancelaw: dict[str, dict[str, dict[int, float]]] = {}
        for f in self.find_log_files("hicstuff/distancelaw"):
            parsed = self.parse_distancelaw(f)
            if parsed is not None:
                s_name = self.clean_s_name(f["s_name"], f)
                self.hicstuff_distancelaw[s_name] = parsed
                self.add_data_source(f, s_name)

        self.hicstuff_distancelaw = self.ignore_samples(self.hicstuff_distancelaw)

        if not self.hicstuff_data and not self.hicstuff_distancelaw:
            raise ModuleNoSamplesFound

        if self.hicstuff_data:
            log.info(f"Found {len(self.hicstuff_data)} pipeline reports")
        if self.hicstuff_distancelaw:
            log.info(f"Found {len(self.hicstuff_distancelaw)} distance-law files")

        if self.hicstuff_data:
            self.hicstuff_stats_table()
            self.hicstuff_read_fate_plot()

        if self.hicstuff_distancelaw:
            self.hicstuff_distance_law_plot()

        if self.hicstuff_data:
            self.write_data_file(self.hicstuff_data, "multiqc_hicstuff")

    def parse_hicstuff_log(self, f):
        """Parse a hicstuff log file and return the stats dict, or None if not found."""
        content = f["f"]
        version_match = re.search(r"^## hicstuff:\s+(v[\d.]+)", content, re.MULTILINE)
        stats_match = re.search(r"(\{'Sample':.*\})\s*$", content, re.MULTILINE)
        if stats_match is None:
            return None
        try:
            stats = ast.literal_eval(stats_match.group(1))
        except (ValueError, SyntaxError) as exc:
            log.error(f"Could not parse hicstuff stats dict: {exc}")
            return None
        if version_match:
            stats["version"] = version_match.group(1)
        return stats

    def parse_distancelaw(self, f) -> Optional[dict[str, dict[int, float]]]:
        """Parse a hicstuff distance-law TSV file.

        Returns a dict mapping chromosome name -> {start_bp: p_s}, or None if no
        data rows were parsed.
        """
        content = f["f"]
        chroms: dict[str, dict[int, float]] = {}
        for line in content.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                start_bp = int(float(parts[0]))
                p_s = float(parts[1])
            except ValueError:
                continue
            chrom = parts[2].strip() if len(parts) >= 3 else "genome"
            chroms.setdefault(chrom, {})[start_bp] = p_s

        return chroms if chroms else None

    def hicstuff_distance_law_plot(self):
        """Generate P(s) distance-law line graph and slope analysis."""
        # Build flat series dicts and colours dict
        plot_data: dict[str, dict[int, float]] = {}
        slope_data: dict[str, dict[int, float]] = {}
        colours: dict[str, str] = {}

        sample_names = sorted(self.hicstuff_distancelaw.keys())
        for idx, s_name in enumerate(sample_names):
            colour = _SAMPLE_COLOURS[idx % len(_SAMPLE_COLOURS)]
            chrom_data = self.hicstuff_distancelaw[s_name]
            chroms = sorted(chrom_data.keys())
            for chrom in chroms:
                series_key = f"{s_name} ({chrom})" if len(chroms) > 1 else s_name
                plot_data[series_key] = chrom_data[chrom]

                # Calculate slopes (derivative in log-log space)
                sorted_bp = sorted(chrom_data[chrom].keys())
                slopes = {}
                for i in range(1, len(sorted_bp)):
                    x1, x2 = sorted_bp[i - 1], sorted_bp[i]
                    y1, y2 = chrom_data[chrom][x1], chrom_data[chrom][x2]
                    # Avoid log of zero or negative values
                    if y1 > 0 and y2 > 0 and x1 > 0 and x2 > 0:
                        # Slope in log-log space: dlog(y)/dlog(x)
                        slope = (math.log10(y2) - math.log10(y1)) / (math.log10(x2) - math.log10(x1))
                        slopes[x2] = slope
                slope_data[series_key] = slopes
                colours[series_key] = colour

        self.add_section(
            name="Distance law P(s)",
            anchor="hicstuff-distance-law-ps",
            description=(
                "Contact probability `P(s)` as a function of genomic distance `s`. "
                "Each curve represents one sample (or chromosome within a sample). "
                "Chromosomes from the same sample share the same colour."
            ),
            helptext="""
                The distance law (also called `P(s)` curve) describes how the contact
                probability between two genomic loci decreases as the distance between
                them increases. For a typical mammalian genome, `P(s)` follows a power law
                with a slope of approximately -1 at intermediate distances.

                When a file contains data for multiple chromosomes, each chromosome is
                shown as a separate line but all lines for the same sample share the same
                colour so they can be visually grouped.
            """,
            plot=linegraph.plot(
                plot_data,
                pconfig={
                    "id": "hicstuff_distance_law_ps",
                    "title": "hicstuff: Distance law P(s)",
                    "xlab": "Genomic distance (bp)",
                    "ylab": "P(s)",
                    "colors": colours,
                    "xlog": True,
                    "ylog": True,
                    "xmin": 1000,
                    "smooth_points": 500,
                    "tt_label": "{point.x} bp: {point.y:.2e}",
                },
            ),
        )

        self.add_section(
            name="Distance law Slope",
            anchor="hicstuff-distance-law-slope",
            description=(
                "Power-law slope (derivative in log-log space) of the distance law. "
                "A slope near `-1` indicates a typical power law."
            ),
            helptext="""
                The slope is calculated as `dlog(P)/dlog(s)` in log-log space.
                A slope of approximately `-1` is typical for mammalian genomes at
                intermediate genomic distances. Steeper slopes (more negative) indicate
                faster decay of contact probability with distance.
            """,
            plot=linegraph.plot(
                slope_data,
                pconfig={
                    "id": "hicstuff_distance_law_slope",
                    "title": "hicstuff: Distance law Slope",
                    "xlab": "Genomic distance (bp)",
                    "ylab": "dlog(P)/dlog(s)",
                    "colors": colours,
                    "xlog": True,
                    "ylog": False,
                    "xmin": 1000,
                    "ymin": -5,
                    "ymax": 2,
                    "smooth_points": 500,
                    "tt_label": "{point.x} bp: {point.y:.2f}",
                },
            ),
        )

    def hicstuff_stats_table(self):
        """Add hicstuff stats to the general stats table."""
        gstats = {}
        for s_name, data in self.hicstuff_data.items():
            n = data["Total read pairs"]
            if not n:
                continue
            recovered = data["Recovered contacts"]
            gstats[s_name] = {
                "total_read_pairs": n,
                "pct_mapped": data["Mapped reads"] / (2.0 * n) * 100.0,
                "pct_pcr_dups": (data["PCR duplicates"] / recovered * 100.0) if recovered else 0.0,
                "pct_final": data["Final contacts"] / n * 100.0,
            }

        headers = {
            "total_read_pairs": {
                "title": "Read pairs",
                "description": "Total number of read pairs sequenced",
                "scale": "Blues",
                "shared_key": "read_count",
                "hidden": True,
            },
            "pct_mapped": {
                "title": "% Mapped",
                "description": "Percentage of single-end reads mapped (Q ≥ 30)",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.1f}",
            },
            "pct_pcr_dups": {
                "title": "% PCR dups",
                "description": "Percentage of recovered contacts identified as PCR duplicates",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "scale": "OrRd",
                "format": "{:,.1f}",
            },
            "pct_final": {
                "title": "% Final",
                "description": "Percentage of read pairs retained as final Hi-C contacts",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.1f}",
            },
        }

        self.general_stats_addcols(gstats, headers, namespace="hicstuff")

    def hicstuff_read_fate_plot(self):
        """Generate the read fate stacked bar chart."""
        keys = {
            "Final contacts": {"color": "#3a7d3d", "name": "Final contacts"},
            "Weirds": {"color": "#e8b372", "name": "Weirds"},
            "Uncuts": {"color": "#f58a50", "name": "Uncuts"},
            "Loops": {"color": "#f95422", "name": "Loops"},
            "PCR duplicates": {"color": "#8b3a39", "name": "PCR duplicates"},
            "Unpaired reads": {"color": "#545454", "name": "Unpaired reads"},
            "Unmapped reads": {"color": "#B0B0B0", "name": "Unmapped reads"},
        }

        plot_data = {}
        for s_name, data in self.hicstuff_data.items():
            if not data["Total read pairs"]:
                continue
            mapped = data["Mapped reads"]
            recovered = data["Recovered contacts"]
            plot_data[s_name] = {
                "Final contacts": data["Final contacts"],
                "Weirds": data["Weirds"],
                "Uncuts": data["Uncuts"],
                "Loops": data["Loops"],
                "PCR duplicates": data["PCR duplicates"],
                "Unpaired reads": max(0.0, mapped / 2.0 - recovered),
                "Unmapped reads": data["Unmapped reads"] / 2.0,
            }

        self.add_section(
            name="Read fate",
            anchor="hicstuff-read-fate",
            description=(
                "Breakdown of read pairs by fate in the Hi-C pipeline. "
                "Values are expressed as read-pair equivalents and sum to the total read pairs."
            ),
            helptext="""
                - **Unmapped reads** - single-end reads that failed to align to the reference genome
                (divided by 2 to give read-pair equivalents).
                - **Unpaired reads** - mapped reads that could not be paired into a valid contact.
                - **PCR duplicates** - contact pairs removed as PCR duplicates.
                - **Loops** - pairs where both ends map to the same restriction fragment (self-ligation artefacts).
                - **Uncuts** - pairs from adjacent fragments that were not digested (dangling ends).
                - **Weirds** - pairs with an unexpected read orientation.
                - **Final contacts** - pairs retained for the Hi-C contact map.
            """,
            plot=bargraph.plot(
                plot_data,
                keys,
                {
                    "id": "hicstuff_read_fate",
                    "title": "hicstuff: Read fate",
                    "ylab": "Read pairs",
                    "cpswitch": True,
                    "cpswitch_c_active": True,
                    "tt_decimals": 0,
                },
            ),
        )
