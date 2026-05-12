import json
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.samtools.stats import parse_samtools_stats_lines
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Parse HiFi-Trimmer JSON summaries and optionally merge sample totals from
    matching samtools stats reports.
    """

    def __init__(self):
        super().__init__(
            name="HiFi-Trimmer",
            anchor="hifi_trimmer",
            href="https://github.com/sanger-tol/hifi-trimmer",
            info="Filters and trims adapter sequences from HiFi reads using BLAST.",
            # doi="",  # No DOI available
        )

        self.hifi_trimmer_data = {}

        # Samtools totals are optional, but when present they are the best denominator
        # for kept percentages and for separating processed from unprocessed signal.
        samtools_data = self._load_samtools_stats()

        for f in self.find_log_files("hifi_trimmer", filehandles=True):
            parsed = self.parse_hifi_trimmer_json(f)
            if parsed is None:
                continue

            s_name, data = parsed
            if s_name in self.hifi_trimmer_data:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

            # Merge totals first so kept percentages prefer full-sample counts over
            # HiFi-Trimmer's processed-only counts.
            self._merge_samtools_totals(s_name, data, samtools_data)
            self._add_kept_percentage(
                s_name,
                data,
                total_key="sample_total_reads",
                processed_key="total_reads_processed",
                removed_key="total_reads_discarded",
                pct_key="pct_reads_kept",
            )
            self._add_kept_percentage(
                s_name,
                data,
                total_key="sample_total_bases",
                processed_key="total_bases_processed",
                removed_key="total_bases_removed",
                pct_key="pct_bases_kept",
            )

            self.hifi_trimmer_data[s_name] = data
            self.add_data_source(f, s_name=s_name)

        self.hifi_trimmer_data = self.ignore_samples(self.hifi_trimmer_data)

        if not self.hifi_trimmer_data:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.hifi_trimmer_data)} reports")

        self.add_software_version(None)

        self.hifi_trimmer_general_stats_table()

        self.add_section(
            name="Read Statistics",
            anchor="hifi_trimmer_reads",
            description="Summary of read processing statistics from HiFi-Trimmer.",
            helptext=(
                "This plot shows discarded, trimmed, and unchanged reads derived from "
                "`total_reads_processed`. If matching samtools stats are found, reads "
                "outside the HiFi-Trimmer input are shown as unprocessed."
            ),
            plot=self.hifi_trimmer_reads_barplot(),
        )

        self.add_section(
            name="Base Statistics",
            anchor="hifi_trimmer_bases",
            description="Summary of base processing statistics from HiFi-Trimmer.",
            helptext=(
                "This plot shows removed and unchanged bases derived from "
                "`total_bases_processed`. If matching samtools stats are found, bases "
                "outside the HiFi-Trimmer input are shown as unprocessed."
            ),
            plot=self.hifi_trimmer_bases_barplot(),
        )

        self.write_data_file(self.hifi_trimmer_data, "multiqc_hifi_trimmer")

    def _load_samtools_stats(self):
        samtools_data = {}
        for f in self.find_log_files("samtools/stats"):
            parsed_data, _, _ = parse_samtools_stats_lines(f["f"])
            # Ignore partial samtools outputs that do not provide a usable total.
            if "raw_total_sequences" not in parsed_data and "total_length" not in parsed_data:
                continue

            samtools_data[f["s_name"]] = parsed_data
            log.debug(f"Found samtools stats for {f['s_name']}")

        return samtools_data

    def _merge_samtools_totals(self, s_name, data, samtools_data):
        stats = samtools_data.get(s_name)
        if stats is None:
            return

        if "raw_total_sequences" in stats:
            data["sample_total_reads"] = stats["raw_total_sequences"]
        if "total_length" in stats:
            data["sample_total_bases"] = stats["total_length"]

        log.debug(f"Merged samtools stats data for {s_name}")

    @staticmethod
    def _add_kept_percentage(s_name, data, total_key, processed_key, removed_key, pct_key):
        # Prefer the full-sample total from samtools, otherwise fall back to the
        # number of reads or bases HiFi-Trimmer says it processed.
        total = data.get(total_key) or data.get(processed_key)
        if not total:
            return

        removed = data[removed_key]
        if removed > total:
            # Likely indicates a wrong samtools file matched, or a real format bug.
            # Skip the pct rather than show a deceptively-clean value.
            log.warning(
                f"HiFi-Trimmer reports {removed_key}={removed} > total ({total}) for '{s_name}'; "
                "skipping kept percentage. Check that the matched samtools stats file is correct."
            )
            return

        data[pct_key] = (total - removed) / total * 100

    @staticmethod
    def _get_unprocessed_total(data, total_key, processed_key):
        total = data.get(total_key)
        if not total:
            return 0

        return max(total - data[processed_key], 0)

    def parse_hifi_trimmer_json(self, f):
        """Parse HiFi-Trimmer JSON output and return summary statistics."""
        required_keys = (
            "total_reads_discarded",
            "total_reads_trimmed",
            "total_bases_removed",
            "total_reads_processed",
            "total_bases_processed",
        )
        try:
            summary = json.load(f["f"]).get("summary")
            data = {key: summary[key] for key in required_keys}
        except (json.JSONDecodeError, AttributeError, TypeError, KeyError) as e:
            log.debug(f"Could not parse HiFi-Trimmer JSON: '{f['fn']}': {e}")
            return None

        s_name = self.clean_s_name(f["fn"], f)

        data["total_reads_unchanged"] = (
            data["total_reads_processed"] - data["total_reads_discarded"] - data["total_reads_trimmed"]
        )
        data["total_bases_unchanged"] = data["total_bases_processed"] - data["total_bases_removed"]

        return s_name, data

    def hifi_trimmer_general_stats_table(self):
        """Add key statistics to the general stats table."""
        headers = {}

        headers["total_reads_processed"] = {
            "title": "Reads processed",
            "description": "Total number of reads processed by HiFi-Trimmer",
            "min": 0,
            "scale": "Blues",
            "shared_key": "read_count",
        }
        headers["pct_reads_kept"] = {
            "title": "% Reads kept",
            "description": "Percentage of reads kept unchanged or trimmed",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.1f}",
        }

        headers["total_bases_processed"] = {
            "title": "Bases processed",
            "description": "Total number of bases processed by HiFi-Trimmer",
            "min": 0,
            "scale": "Purples",
            "shared_key": "base_count",
            "hidden": True,
        }
        headers["pct_bases_kept"] = {
            "title": "% Bases kept",
            "description": "Percentage of bases kept (not removed) by HiFi-Trimmer",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.1f}",
            "hidden": True,
        }

        self.general_stats_addcols(self.hifi_trimmer_data, headers)

    def hifi_trimmer_reads_barplot(self):
        """Generate a bar plot showing read statistics."""
        plot_data = {}
        for s_name, data in self.hifi_trimmer_data.items():
            entry = {
                "total_reads_unchanged": data["total_reads_unchanged"],
                "total_reads_trimmed": data["total_reads_trimmed"],
                "total_reads_discarded": data["total_reads_discarded"],
            }

            unprocessed = self._get_unprocessed_total(
                data,
                total_key="sample_total_reads",
                processed_key="total_reads_processed",
            )
            if unprocessed > 0:
                entry["unprocessed_reads"] = unprocessed

            plot_data[s_name] = entry

        cats = {
            "unprocessed_reads": {"name": "Reads unprocessed", "color": "#95a5a6"},
            "total_reads_unchanged": {"name": "Reads unchanged", "color": "#2ecc71"},
            "total_reads_trimmed": {"name": "Reads trimmed", "color": "#f39c12"},
            "total_reads_discarded": {"name": "Reads discarded", "color": "#e74c3c"},
        }

        config = {
            "id": "hifi_trimmer_reads_plot",
            "title": "HiFi-Trimmer: Read Statistics",
            "ylab": "Number of Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_zero_cats": True,
        }

        return bargraph.plot(plot_data, cats, config)

    def hifi_trimmer_bases_barplot(self):
        """Generate a bar plot showing base statistics."""
        plot_data = {}
        for s_name, data in self.hifi_trimmer_data.items():
            entry = {
                "total_bases_unchanged": data["total_bases_unchanged"],
                "total_bases_removed": data["total_bases_removed"],
            }

            unprocessed = self._get_unprocessed_total(
                data,
                total_key="sample_total_bases",
                processed_key="total_bases_processed",
            )
            if unprocessed > 0:
                entry["unprocessed_bases"] = unprocessed

            plot_data[s_name] = entry

        cats = {
            "unprocessed_bases": {"name": "Bases unprocessed", "color": "#95a5a6"},
            "total_bases_unchanged": {"name": "Bases unchanged", "color": "#3498db"},
            "total_bases_removed": {"name": "Bases removed", "color": "#e74c3c"},
        }

        config = {
            "id": "hifi_trimmer_bases_plot",
            "title": "HiFi-Trimmer: Base Statistics",
            "ylab": "Number of Bases",
            "cpswitch_counts_label": "Number of Bases",
            "hide_zero_cats": True,
        }

        return bargraph.plot(plot_data, cats, config)
