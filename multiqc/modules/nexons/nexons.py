import json
import logging
from typing import Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table
from multiqc.types import SectionAlert

log = logging.getLogger(__name__)

OUTCOMES = {
    "Total_Reads": "Total Reads",
    "No_Alignment": "No Alignment",
    "Primary_Alignment": "Primary Alignment",
    "Secondary_Alignment": "Secondary Alignment",
    "No_Gene": "No Gene Match",
    "No_Hit": "No Transcript Hit",
    "Gene": "Unique Gene Match",
    "Multi_Gene": "Multi Gene Match",
    "Partial": "Partial Unique Transcript Match",
    "Unique": "Full Transcript Match",
    "Same_Strand_Hit": "Same Strand Matches",
    "Opposing_Strand_Hit": "Opposing Strand Matches",
}


class MultiqcModule(BaseMultiqcModule):
    """
    Nexons quantifies long RNA sequencing reads against annotated genes and transcripts.
    This module reads the JSON in `*nexons_stats.txt` files and reproduces the metrics
    and distributions from the nexons combined QC report: read fate, gene and transcript
    matching, alignment directionality, read lengths, transcript coverage and exon flexibility.

    Sample names come from the input BAM path stored in the statistics file, subject to
    MultiQC sample name cleaning and `use_filename_as_sample_name` configuration.
    Match counts overlap: partial and full transcript matches are subsets of gene matches.
    All outcome percentages use total reads as the denominator, including directionality.
    Distribution counts and bins are preserved, including the collapsed read length tail.
    Nexons statistics files do not contain a software version or run parameters.
    """

    def __init__(self):
        super().__init__(
            name="Nexons",
            anchor="nexons",
            href="https://github.com/s-andrews/nexons",
            info="Quantifies long RNA sequencing reads against annotated genes and transcripts.",
            doi=None,
            license="GNU General Public License v3.0",
            license_url="https://github.com/s-andrews/nexons/blob/master/LICENSE",
        )
        self.nexons_data: dict[str, dict] = {}
        for f in self.find_log_files("nexons"):
            try:
                data = self.parse_stats(f["f"])
            except (ValueError, TypeError, KeyError) as exc:
                raise ValueError(f"Invalid nexons statistics in {f['root']}/{f['fn']}: {exc}") from exc
            s_name = self.clean_s_name(data["file"], f)
            if s_name in self.nexons_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.nexons_data[s_name] = data
            self.add_data_source(f, s_name=s_name)
            self.add_software_version(None, s_name)

        self.nexons_data = self.ignore_samples(self.nexons_data)
        if not self.nexons_data:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.nexons_data)} reports")
        self.outcomes = {s: d["outcomes"] for s, d in self.nexons_data.items()}
        self.percentages = {
            s: {k: 100 * v / d["Total_Reads"] if d["Total_Reads"] else 0 for k, v in d.items()}
            for s, d in self.outcomes.items()
        }
        self.add_general_stats()
        self.add_section(
            name="Summary",
            anchor="nexons_summary",
            description="All reported outcome counts. Gene and transcript match categories overlap.",
            plot=table.plot(
                self.outcomes,
                {k: {"title": title, "format": "{:,.0f}"} for k, title in OUTCOMES.items()},
                {"id": "nexons_summary_table", "title": "Nexons: Summary"},
            ),
        )
        self.add_outcome_plot("read_fate", "Read Fate", ["No_Alignment", "Primary_Alignment", "Secondary_Alignment"], "Breakdown of types of alignment seen in the input BAM file")
        for key in ("Gene", "Partial", "Unique", "Multi_Gene", "No_Gene"):
            self.add_outcome_plot(key.lower(), OUTCOMES[key], [key], "Reads assigned to features with different degrees of specificity. Individual reads can be in multiple classes, so all transcipt matching reads will also be gene matching.")
        self.add_outcome_plot("directionality", "Alignment Directionality", ["Same_Strand_Hit", "Opposing_Strand_Hit"],"Direction of matches relative to annotated features")
        self.add_distributions()
        self.write_data_file(self.outcomes, "multiqc_nexons")
        self.write_data_file(self.nexons_data, "multiqc_nexons_distributions", data_format="json")

    @staticmethod
    def parse_stats(contents: str) -> dict:
        data = json.loads(contents)
        if not isinstance(data, dict) or not isinstance(data["file"], str) or not data["file"]:
            raise ValueError("Expected an object with an input filename")

        def check_count(value):
            if type(value) is not int or value < 0:
                raise ValueError(f"Expected a non-negative integer count, found {value!r}")

        for key in OUTCOMES:
            check_count(data["outcomes"][key])
        for key in ("read_lengths", "coverage"):
            if not isinstance(data[key], list):
                raise TypeError(f"Expected a list for {key}")
        for pair in data["read_lengths"]:
            if not isinstance(pair, list) or len(pair) != 2:
                raise ValueError("Expected read length/count pairs")
            check_count(pair[0])
            check_count(pair[1])
        for value in data["coverage"]:
            check_count(value)
        for key in ("end_flex", "inner_flex"):
            if not isinstance(data[key], dict):
                raise TypeError(f"Expected an object for {key}")
            for position, value in data[key].items():
                int(position)
                check_count(value)
        return data

    def add_general_stats(self):
        data = {
            s: {"Total_Reads": d["Total_Reads"], **{k: self.percentages[s][k] for k in ("Gene", "Unique", "Partial")}}
            for s, d in self.outcomes.items()
        }
        headers: dict[str, dict[str, Any]] = {
            "Total_Reads": {"title": "Reads", "description": "Total reads", "shared_key": "read_count"},
            **{
                k: {
                    "title": title,
                    "description": f"{OUTCOMES[k]} as a percentage of all reads",
                    "suffix": "%",
                    "min": 0,
                    "max": 100,
                    "scale": "Blues",
                    "hidden": k == "Partial",
                }
                for k, title in (
                    ("Gene", "% Gene"),
                    ("Unique", "% Full Transcript"),
                    ("Partial", "% Partial Transcript"),
                )
            },
        }
        self.general_stats_addcols(data, headers)

    def add_outcome_plot(self, anchor: str, title: str, keys: list, description: str):
        # Explicit datasets avoid MultiQC normalising overlapping counts to their sum.
        self.add_section(
            name=title,
            anchor=f"nexons_{anchor}",
            description=description,
            plot=bargraph.plot(
                [
                    {s: {k: d[k] for k in keys} for s, d in dataset.items()}
                    for dataset in (self.percentages, self.outcomes)
                ],
                {k: {"name": OUTCOMES[k]} for k in keys},
                {
                    "id": f"nexons_{anchor}_plot",
                    "title": f"Nexons: {title}",
                    "cpswitch": False,
                    "hide_zero_cats": False,
                    "data_labels": [
                        {"name": "Percentage of all reads", "ylab": "Percentage of all reads (%)", "ymax": 100},
                        {"name": "Counts", "ylab": "Count", "ymax": None},
                    ],
                },
            ),
        )

    def add_distributions(self):
        for key, title, xlab, ylab, description in (
            ("read_lengths", "Read Lengths", "Read length (bp, reported bins)", "Read count","Lengths of all reads processed (200bp bin size)"),
            ("coverage", "Transcript Coverage", "Average percentile coverage over transcripts (5′ to 3′)", "Read count","Average relative coverage of transcript area from 5' to 3' - useful for detecting coverage bias in your libraries."),
            ("inner_flex", "Inner Exon Flex", "Observed distances from annotated exon junctions (bp)", "Junction count","Distribution of distances of observed exon ends compared to the annotation in the GTF file. Limited to whatever value was set for 'flex' in the original analysis."),
            ("end_flex", "Transcript End Flex", "Observed distances from annotated transcript ends (bp)", "Junction count","Distribution of distances of observed transcript ends compared to the annotation in the GTF file.  Limited to whatever values was set for 'endflex' in the original analysis"),
        ):
            data = {}
            for sample, stats in self.nexons_data.items():
                values = stats[key]
                if key == "coverage":
                    points = dict(enumerate(values))
                elif key == "read_lengths":
                    points = dict(values)
                else:
                    points = {int(x): y for x, y in values.items()}
                if points:
                    data[sample] = dict(sorted(points.items()))
            empty_samples = [s for s in self.nexons_data if s not in data]
            self.add_section(
                name=title,
                anchor=f"nexons_{key}",
                description=description,
                alerts=[SectionAlert(level="info", message="No observations reported.", affected_samples=empty_samples)]
                if empty_samples
                else [],
                plot=linegraph.plot(
                    data,
                    {
                        "id": f"nexons_{key}_plot",
                        "title": f"Nexons: {title}",
                        "xlab": xlab,
                        "ylab": ylab,
                        "ymin": 0,
                        "hide_empty": False,
                        "smooth_points": max(len(points) for points in data.values()),
                    },
                )
                if data
                else None,
            )
