import logging
from collections.abc import Mapping
from typing import Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import ColumnDict

from .parser import (
    ParsedReport,
    SummaryValidationError,
    count_mapping,
    length_label,
    mapping,
    non_negative_integer,
    parse_summary,
)

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    `vcf-sv-stats` calculates standards-aware descriptive statistics for structural-variant
    and copy-number VCF or BCF callsets.

    The module discovers `*.vcf-sv-stats.json` summaries written by `vcf-sv-stats stats`
    or `vcf-sv-stats run`. It validates the schema, content signature, and RFC 8785 payload
    digest before displaying any value. The digest is an integrity check, not proof of
    authorship.

    Report identifiers remain separate from analysis-unit context and VCF sample columns.
    Counts are descriptive callset statistics. They are not precision, recall, concordance,
    or truth-set accuracy measurements.
    """

    def __init__(self):
        super().__init__(
            name="vcf-sv-stats",
            anchor="vcf-sv-stats",
            href="https://pypi.org/project/vcf-sv-stats/",
            info="Summarizes structural-variant and copy-number VCF or BCF callsets.",
            doi="",
        )

        records: dict[str, ParsedReport] = {}
        self.duplicate_sources: list[str] = []
        for f in self.find_log_files("vcf_sv_stats"):
            for parsed in parse_summary(f["f"], f["fn"]):
                previous = records.get(parsed.report_id)
                if previous is not None:
                    if previous.payload_sha256 != parsed.payload_sha256:
                        raise SummaryValidationError(f"Conflicting vcf-sv-stats report identifier: {parsed.report_id}")
                    self.duplicate_sources.append(f["fn"])
                    continue
                records[parsed.report_id] = parsed
                self.add_data_source(f, s_name=parsed.report_id)
                self.add_software_version(parsed.producer_version, sample=parsed.report_id)

        self.vcf_sv_stats_data: dict[str, dict[str, Any]] = {
            report_id: {
                "callset": parsed.callset,
                "validation": parsed.validation,
                "report": parsed.report,
                "report_payload_sha256": parsed.payload_sha256,
            }
            for report_id, parsed in records.items()
        }
        self.vcf_sv_stats_data = self.ignore_samples(self.vcf_sv_stats_data)
        if not self.vcf_sv_stats_data:
            raise ModuleNoSamplesFound
        records = {key: records[key] for key in self.vcf_sv_stats_data}
        log.info(f"Found {len(records)} vcf-sv-stats reports")

        self._add_general_stats(records)
        self._add_overview(records)
        self._add_sv_types(records)
        self._add_length_section(records)
        self._add_flat_count_section(records, "filters", "Record filter states")
        self._add_breakend_section(records)
        self._add_merged_support_section(records)
        self._add_flat_count_section(records, "genotypes", "Genotype states")
        self._add_flat_count_section(records, "copy_number", "Declared copy number")
        self._add_validation_section(records)
        self._add_normalization_section(records)
        self._add_provenance_section(records)

        self.write_data_file(self.vcf_sv_stats_data, "multiqc_vcf_sv_stats")

    @staticmethod
    def _statistics(parsed: ParsedReport) -> dict[str, Any]:
        return mapping(parsed.report["statistics"], "report.statistics")

    def _add_general_stats(self, records: Mapping[str, ParsedReport]) -> None:
        data: dict[str, dict[str, Any]] = {}
        for report_id, parsed in records.items():
            statistics = self._statistics(parsed)
            breakends = mapping(statistics["breakends"], "report.statistics.breakends")
            diagnostic_counts = count_mapping(parsed.validation["diagnostic_counts"], "validation.diagnostic_counts")
            bnd_total = non_negative_integer(breakends["total"], "breakends.total")
            unresolved = non_negative_integer(
                breakends["without_declared_mate"], "breakends.without_declared_mate"
            ) + non_negative_integer(
                breakends["unresolved_mate_references"],
                "breakends.unresolved_mate_references",
            )
            data[report_id] = {
                "events": non_negative_integer(statistics["events"]["resolved"], "events.resolved"),
                "alleles": non_negative_integer(statistics["alleles"]["total"], "alleles.total"),
                "records": non_negative_integer(statistics["source_records"]["total"], "source_records.total"),
                "errors": diagnostic_counts.get("error", 0) + diagnostic_counts.get("fatal", 0),
                "unresolved_bnd_pct": (0 if bnd_total == 0 else unresolved * 100 / bnd_total),
            }
        all_headers: dict[str, ColumnDict] = {
            "events": {
                "title": "SV events",
                "description": "Resolved structural-variant events",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "alleles": {
                "title": "SV alleles",
                "description": ("Parsed alternate alleles; this is not a source-record count"),
                "scale": "Purples",
                "format": "{:,.0f}",
            },
            "records": {
                "title": "VCF records",
                "description": "Parsed source VCF or BCF records",
                "scale": "Greens",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "errors": {
                "title": "Validation errors",
                "description": ("Error and fatal diagnostics from vcf-sv-stats validation"),
                "scale": "Reds",
                "format": "{:,.0f}",
            },
            "unresolved_bnd_pct": {
                "title": "Unresolved BND",
                "description": "Breakends without a declared or resolvable mate",
                "scale": "OrRd",
                "suffix": "%",
                "min": 0,
                "max": 100,
                "hidden": True,
            },
        }
        headers = self.get_general_stats_headers(all_headers=all_headers)
        if headers:
            self.general_stats_addcols(data, headers, namespace="vcf-sv-stats")

    def _add_overview(self, records: Mapping[str, ParsedReport]) -> None:
        data: dict[str, dict[str, Any]] = {}
        for report_id, parsed in records.items():
            statistics = self._statistics(parsed)
            data[report_id] = {
                "records": statistics["source_records"]["total"],
                "alleles": statistics["alleles"]["total"],
                "events": statistics["events"]["resolved"],
                "producer": parsed.callset["producer"]["producer"],
                "adapter_status": parsed.callset["producer"]["status"],
                "analysis_unit": parsed.report["analysis_unit"]["status"],
                "vcf_samples": (", ".join(parsed.report["mapped_vcf_sample_ids"]) or "Unmapped"),
            }
        self.add_section(
            name="Callset overview",
            anchor="vcf-sv-stats-overview",
            description=(
                "Records, alternate alleles, and resolved events are distinct statistical "
                "grains. The report identifier is a compatibility key, not a biological "
                "sample identity."
            ),
            plot=table.plot(
                data,
                {
                    "records": {
                        "title": "VCF records",
                        "format": "{:,.0f}",
                        "scale": "Greens",
                    },
                    "alleles": {
                        "title": "ALT alleles",
                        "format": "{:,.0f}",
                        "scale": "Purples",
                    },
                    "events": {
                        "title": "Resolved events",
                        "format": "{:,.0f}",
                        "scale": "Blues",
                    },
                    "producer": {"title": "Producer"},
                    "adapter_status": {"title": "Adapter status"},
                    "analysis_unit": {"title": "Analysis-unit status"},
                    "vcf_samples": {"title": "Mapped VCF samples"},
                },
                pconfig={
                    "id": "vcf-sv-stats-overview-table",
                    "title": "vcf-sv-stats: callset overview",
                },
            ),
        )

    def _add_sv_types(self, records: Mapping[str, ParsedReport]) -> None:
        data: dict[str, dict[str, int]] = {}
        keys: set[str] = set()
        for report_id, parsed in records.items():
            alleles = mapping(self._statistics(parsed)["alleles"], "statistics.alleles")
            counts = count_mapping(alleles["types"], "statistics.alleles.types")
            data[report_id] = counts
            keys.update(counts)
        self._add_bargraph(
            data,
            keys,
            title="SV type spectrum",
            anchor="vcf-sv-stats-types",
            plot_id="vcf-sv-stats-types-plot",
            ylab="Alleles",
        )

    def _add_flat_count_section(
        self,
        records: Mapping[str, ParsedReport],
        key: str,
        title: str,
    ) -> None:
        data: dict[str, dict[str, int]] = {}
        keys: set[str] = set()
        for report_id, parsed in records.items():
            counts = count_mapping(self._statistics(parsed)[key], f"statistics.{key}")
            data[report_id] = counts
            keys.update(counts)
        anchor_key = key.replace("_", "-")
        self._add_bargraph(
            data,
            keys,
            title=title,
            anchor=f"vcf-sv-stats-{anchor_key}",
            plot_id=f"vcf-sv-stats-{anchor_key}-plot",
            ylab="Count",
        )

    def _add_bargraph(
        self,
        data: dict[str, dict[str, int]],
        keys: set[str],
        *,
        title: str,
        anchor: str,
        plot_id: str,
        ylab: str,
    ) -> None:
        categories = {key: {"name": str(key).replace("_", " ")} for key in sorted(keys)}
        self.add_section(
            name=title,
            anchor=anchor,
            description=(
                "Counts retain the metric scope and denominator declared in the producer "
                "summary. They are descriptive statistics, not truth-set accuracy measures."
            ),
            plot=bargraph.plot(
                data,
                categories,
                pconfig={
                    "id": plot_id,
                    "title": f"vcf-sv-stats: {title}",
                    "ylab": ylab,
                    "cpswitch": False,
                },
            ),
        )

    def _add_length_section(self, records: Mapping[str, ParsedReport]) -> None:
        data: dict[str, dict[str, int]] = {}
        expected_labels: tuple[str, ...] | None = None
        for report_id, parsed in records.items():
            length = mapping(self._statistics(parsed)["length_bp"], "statistics.length_bp")
            labels = tuple(length_label(boundary) for boundary in length["boundaries"])
            if expected_labels is not None and labels != expected_labels:
                raise SummaryValidationError("Conflicting length histogram boundaries")
            expected_labels = labels
            data[report_id] = dict(zip(labels, length["counts"]))
        categories = {label: {"name": label} for label in expected_labels or ()}
        self.add_section(
            name="Structural-variant length",
            anchor="vcf-sv-stats-length",
            description=(
                "Alleles with an applicable length, grouped using the producer's fixed "
                "histogram policy. Missing, invalid, and not-applicable values are not zero."
            ),
            plot=bargraph.plot(
                data,
                categories,
                pconfig={
                    "id": "vcf-sv-stats-length-plot",
                    "title": "vcf-sv-stats: structural-variant length",
                    "ylab": "Alleles",
                    "cpswitch": False,
                },
            ),
        )

    def _add_breakend_section(self, records: Mapping[str, ParsedReport]) -> None:
        data = {report_id: self._statistics(parsed)["breakends"] for report_id, parsed in records.items()}
        self.add_section(
            name="Breakend relationships",
            anchor="vcf-sv-stats-breakends",
            description=("Pair resolution follows explicit record and relationship evidence only."),
            plot=table.plot(
                data,
                {
                    "total": {
                        "title": "Breakends",
                        "format": "{:,.0f}",
                        "scale": "Blues",
                    },
                    "reciprocal_pairs": {
                        "title": "Reciprocal pairs",
                        "format": "{:,.0f}",
                    },
                    "without_declared_mate": {
                        "title": "No declared mate",
                        "format": "{:,.0f}",
                        "scale": "OrRd",
                    },
                    "unresolved_mate_references": {
                        "title": "Unresolved mate",
                        "format": "{:,.0f}",
                        "scale": "Reds",
                    },
                },
                pconfig={
                    "id": "vcf-sv-stats-breakends-table",
                    "title": "vcf-sv-stats: breakends",
                },
            ),
        )

    def _add_merged_support_section(self, records: Mapping[str, ParsedReport]) -> None:
        data: dict[str, dict[str, Any]] = {}
        for report_id, parsed in records.items():
            support = mapping(self._statistics(parsed)["merged_support"], "statistics.merged_support")
            source_counts = count_mapping(
                support["source_count_states"], "statistics.merged_support.source_count_states"
            )
            consistency = count_mapping(
                support["vector_count_consistency"],
                "statistics.merged_support.vector_count_consistency",
            )
            data[report_id] = {
                "producer_kind": parsed.callset["producer_kind"],
                "status": support["status"],
                "records_with_support": support["records_with_support"],
                "records_without_support": support["records_without_support"],
                "declared_source_counts": ", ".join(sorted(source_counts)) or "Not reported",
                "inconsistent_vectors": consistency.get("inconsistent", 0),
            }
        self.add_section(
            name="Merged support provenance",
            anchor="vcf-sv-stats-merged-support",
            description=(
                "Merger support fields describe source provenance. They are not biological "
                "replication, precision, recall, concordance, or truth evidence."
            ),
            plot=table.plot(
                data,
                {
                    "producer_kind": {"title": "Producer kind"},
                    "status": {"title": "Support status"},
                    "records_with_support": {
                        "title": "Records with support",
                        "format": "{:,.0f}",
                        "scale": "Blues",
                    },
                    "records_without_support": {
                        "title": "Records without support",
                        "format": "{:,.0f}",
                        "scale": "Oranges",
                    },
                    "declared_source_counts": {"title": "Declared source-count states"},
                    "inconsistent_vectors": {
                        "title": "Count/vector conflicts",
                        "format": "{:,.0f}",
                        "scale": "Reds",
                    },
                },
                pconfig={
                    "id": "vcf-sv-stats-merged-support-table",
                    "title": "vcf-sv-stats: merger support provenance",
                },
            ),
        )

    def _add_validation_section(self, records: Mapping[str, ParsedReport]) -> None:
        data: dict[str, dict[str, Any]] = {}
        for report_id, parsed in records.items():
            states = parsed.validation["states"]
            counts = count_mapping(parsed.validation["diagnostic_counts"], "validation.diagnostic_counts")
            data[report_id] = {
                "container": states["container_state"],
                "parse": states["parse_state"],
                "conformance": states["vcf_conformance_state"],
                "sv_semantics": states["sv_semantic_state"],
                "statistics": states["statistics_state"],
                "warnings": counts.get("warning", 0),
                "errors": counts.get("error", 0) + counts.get("fatal", 0),
            }
        self.add_section(
            name="Validation states",
            anchor="vcf-sv-stats-validation",
            description=(
                "Container, parse, conformance, structural-variant semantics, and statistics states remain independent."
            ),
            plot=table.plot(
                data,
                {
                    "container": {"title": "Container"},
                    "parse": {"title": "Parse"},
                    "conformance": {"title": "VCF conformance"},
                    "sv_semantics": {"title": "SV semantics"},
                    "statistics": {"title": "Statistics"},
                    "warnings": {
                        "title": "Warnings",
                        "format": "{:,.0f}",
                        "scale": "Oranges",
                    },
                    "errors": {
                        "title": "Errors",
                        "format": "{:,.0f}",
                        "scale": "Reds",
                    },
                },
                pconfig={
                    "id": "vcf-sv-stats-validation-table",
                    "title": "vcf-sv-stats: validation states",
                },
            ),
        )

    def _add_normalization_section(self, records: Mapping[str, ParsedReport]) -> None:
        data = {
            report_id: {
                "operation_safety": parsed.validation["states"]["operation_safety_state"],
                "adapter_status": parsed.callset["producer"]["status"],
            }
            for report_id, parsed in records.items()
        }
        self.add_section(
            name="Normalization safety",
            anchor="vcf-sv-stats-normalization",
            description=(
                "Operation safety reports whether rewriting is supported by the observed "
                "evidence. It does not imply that normalization occurred."
            ),
            plot=table.plot(
                data,
                {
                    "operation_safety": {"title": "Operation safety"},
                    "adapter_status": {"title": "Adapter status"},
                },
                pconfig={
                    "id": "vcf-sv-stats-normalization-table",
                    "title": "vcf-sv-stats: normalization safety",
                },
            ),
        )

    def _add_provenance_section(self, records: Mapping[str, ParsedReport]) -> None:
        data: dict[str, dict[str, Any]] = {}
        for report_id, parsed in records.items():
            producer = parsed.callset["producer"]
            analysis = parsed.report["analysis_unit"]
            data[report_id] = {
                "producer": producer["producer"],
                "version": producer.get("version") or "Unknown",
                "adapter": producer["adapter_id"],
                "adapter_status": producer["status"],
                "analysis_unit_id": analysis.get("analysis_unit_id") or "Unresolved",
                "display_id": analysis.get("display_id") or "Unresolved",
                "algorithm_id": analysis.get("algorithm_id") or "Unresolved",
                "vcf_samples": (", ".join(parsed.report["mapped_vcf_sample_ids"]) or "Unmapped"),
            }
        self.add_section(
            name="Provenance and analysis context",
            anchor="vcf-sv-stats-provenance",
            description=(
                "Producer support is provenance, not accuracy evidence. Analysis-unit "
                "identifiers and mapped VCF sample columns remain separate."
            ),
            plot=table.plot(
                data,
                {
                    "producer": {"title": "Producer"},
                    "version": {"title": "Producer version"},
                    "adapter": {"title": "Adapter"},
                    "adapter_status": {"title": "Adapter status"},
                    "analysis_unit_id": {"title": "Analysis unit"},
                    "display_id": {"title": "Display ID"},
                    "algorithm_id": {"title": "Algorithm"},
                    "vcf_samples": {"title": "Mapped VCF samples"},
                },
                pconfig={
                    "id": "vcf-sv-stats-provenance-table",
                    "title": "vcf-sv-stats: provenance",
                },
            ),
        )
