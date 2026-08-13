import copy
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest
import rfc8785

from multiqc import config, report, validation
from multiqc.modules.vcf_sv_stats import MultiqcModule
from multiqc.modules.vcf_sv_stats.parser import SummaryValidationError, parse_summary


@pytest.fixture(autouse=True)
def reset_multiqc_state():
    yield
    report.reset()
    config.reset()
    validation.reset()


def _metric_contract(scope: str, unit: str, comparability: str) -> dict[str, str]:
    return {
        "scope": scope,
        "denominator": f"all_{scope}",
        "unit": unit,
        "comparability": comparability,
    }


def _statistics() -> dict[str, Any]:
    return {
        "source_records": {"total": 12},
        "alleles": {"total": 14, "types": {"BND": 2, "DEL": 7, "DUP": 5}},
        "events": {"resolved": 13},
        "breakends": {
            "total": 2,
            "reciprocal_pairs": 1,
            "without_declared_mate": 0,
            "unresolved_mate_references": 0,
        },
        "filters": {"PASS": 10, "filtered_any": 2, "LowQuality": 2},
        "genotypes": {
            "heterozygous_or_other_alt": 8,
            "homozygous_alt": 2,
            "no_call": 2,
        },
        "copy_number": {"1": 7, "2": 3, "missing": 2},
        "merged_support": {
            "status": "not_declared",
            "declared_fields": [],
            "records_with_support": 0,
            "records_without_support": 0,
            "supporting_sources": {},
            "source_count_states": {},
            "vector_count_consistency": {},
            "interpretation": "merger_provenance_only",
        },
        "length_bp": {
            "boundaries": [[0, 50], [50, 100], [100, None]],
            "counts": [2, 5, 5],
            "n": 12,
            "missing": 0,
            "invalid": 0,
            "not_applicable": 2,
            "policy_id": "vss-bins/1",
        },
        "metric_contracts": {
            "source_records": _metric_contract("source_records", "records", "canonical-observation/1"),
            "alleles": _metric_contract("alternate_alleles", "alleles", "canonical-observation/1"),
            "events": _metric_contract("resolved_events", "events", "event-resolution/1"),
            "breakends": _metric_contract("breakend_observations", "breakends", "event-resolution/1"),
            "filters": _metric_contract("source_records", "records", "vcf-filter-state/1"),
            "genotypes": _metric_contract("sample_calls", "sample_calls", "vcf-genotype-state/1"),
            "copy_number": _metric_contract("sample_calls", "declared_copy_number", "declared-cn/1"),
            "length_bp": _metric_contract("alternate_alleles", "base_pairs", "vss-bins/1"),
            "merged_support": _metric_contract("source_records", "records", "merger-support-provenance/1"),
        },
    }


def make_summary(report_count: int = 1) -> dict[str, Any]:
    statistics = _statistics()
    value: dict[str, Any] = {
        "schema_name": "vcf-sv-stats.summary",
        "schema_version": "1.0.0",
        "content_signature": "vcf-sv-stats:summary:1",
        "producer": {"name": "vcf-sv-stats", "version": "0.8.0"},
        "input": {
            "display_name": "neutral.hg002.vcf.gz",
            "complete": True,
            "sha256": "a" * 64,
        },
        "callset": {
            "record_count": 12,
            "allele_count": 14,
            "producer_kind": "caller",
            "producer": {
                "adapter_id": "urn:vcf-sv-stats:adapter:generic:1",
                "producer": "unknown",
                "version": None,
                "status": "supported",
            },
        },
        "validation": {
            "diagnostic_counts": {},
            "states": {
                "container_state": "accepted",
                "parse_state": "accepted",
                "vcf_conformance_state": "conformant",
                "sv_semantic_state": "consistent",
                "operation_safety_state": "safe",
                "statistics_state": "complete",
            },
        },
        "statistics": statistics,
        "reports": [
            {
                "report_id": f"vss1-test-{index:04d}",
                "analysis_unit": {"status": "unresolved"},
                "mapped_vcf_sample_ids": [],
                "statistics": copy.deepcopy(statistics),
            }
            for index in range(report_count)
        ],
    }
    value["payload_sha256"] = hashlib.sha256(rfc8785.dumps(value)).hexdigest()
    return value


def render(value: dict[str, Any]) -> str:
    return json.dumps(value, sort_keys=True)


def resign(value: dict[str, Any]) -> dict[str, Any]:
    value.pop("payload_sha256", None)
    value["payload_sha256"] = hashlib.sha256(rfc8785.dumps(value)).hexdigest()
    return value


def test_parser_validates_digest_schema_and_nested_contract() -> None:
    value = make_summary()
    parsed = parse_summary(render(value), "valid.vcf-sv-stats.json")
    assert len(parsed) == 1
    assert parsed[0].report_id == "vss1-test-0000"
    assert parsed[0].report["statistics"]["events"]["resolved"] == 13

    tampered = copy.deepcopy(value)
    tampered["reports"][0]["statistics"]["events"]["resolved"] = 99
    with pytest.raises(SummaryValidationError, match="digest does not match"):
        parse_summary(render(tampered), "tampered.vcf-sv-stats.json")

    future = copy.deepcopy(value)
    future["schema_version"] = "2.0.0"
    with pytest.raises(SummaryValidationError, match="schema major"):
        parse_summary(render(future), "future.vcf-sv-stats.json")

    malformed = copy.deepcopy(value)
    malformed["reports"][0]["statistics"]["length_bp"]["n"] = 11
    resign(malformed)
    with pytest.raises(SummaryValidationError, match="counts do not reconcile"):
        parse_summary(render(malformed), "malformed.vcf-sv-stats.json")

    missing = copy.deepcopy(value)
    del missing["reports"][0]["statistics"]["events"]
    resign(missing)
    with pytest.raises(SummaryValidationError, match="required nested field"):
        parse_summary(render(missing), "missing.vcf-sv-stats.json")


@pytest.mark.parametrize(
    "case",
    [
        "unknown-producer",
        "cnv",
        "merged",
        "multi-unit",
        "unresolved",
        "invalid",
        "normalized-name",
        "missing-optional-context",
    ],
)
def test_parser_golden_context_matrix(case: str) -> None:
    value = make_summary(report_count=2 if case == "multi-unit" else 1)
    if case == "unknown-producer":
        value["callset"]["producer"]["status"] = "provisional"
    elif case == "cnv":
        value["reports"][0]["statistics"]["alleles"]["types"] = {"CNV": 14}
    elif case == "merged":
        value["callset"]["producer_kind"] = "merger"
        value["reports"][0]["statistics"]["merged_support"] = {
            "status": "present",
            "declared_fields": ["SUPP", "SUPP_VEC"],
            "records_with_support": 12,
            "records_without_support": 0,
            "supporting_sources": {"1": 7, "2": 5},
            "source_count_states": {"6": 12},
            "vector_count_consistency": {"consistent": 12},
            "interpretation": "merger_provenance_only",
        }
    elif case == "invalid":
        value["validation"]["diagnostic_counts"] = {"error": 2, "warning": 1}
        value["validation"]["states"]["vcf_conformance_state"] = "nonconformant"
    elif case == "normalized-name":
        value["input"]["display_name"] = "neutral.normalized.vcf.gz"
    elif case == "missing-optional-context":
        assert value["reports"][0]["analysis_unit"] == {"status": "unresolved"}
    resign(value)
    parsed = parse_summary(render(value), f"{case}.vcf-sv-stats.json")
    assert len(parsed) == (2 if case == "multi-unit" else 1)


def test_parser_handles_one_thousand_reports_deterministically() -> None:
    value = make_summary(report_count=1000)
    first = parse_summary(render(value), "scale.vcf-sv-stats.json")
    second = parse_summary(render(value), "scale.vcf-sv-stats.json")
    assert first == second
    assert len(first) == 1000
    assert first[0].report_id == "vss1-test-0000"
    assert first[-1].report_id == "vss1-test-0999"


def _write_summary(path: Path, value: dict[str, Any]) -> Path:
    path.write_text(render(value), encoding="utf-8")
    return path


def test_module_renders_all_sections_and_preserves_report_identity(tmp_path: Path) -> None:
    summary = make_summary()
    summary["reports"][0]["analysis_unit"] = {
        "status": "resolved",
        "analysis_unit_id": "analysis-001",
        "display_id": "Demonstration",
        "algorithm_id": "generic-sv-1",
    }
    summary["reports"][0]["mapped_vcf_sample_ids"] = ["HG002"]
    resign(summary)
    source = _write_summary(tmp_path / "valid.vcf-sv-stats.json", summary)

    report.analysis_files = [str(source)]
    report.search_files(["vcf_sv_stats"])
    module = MultiqcModule()

    assert list(module.vcf_sv_stats_data) == ["vss1-test-0000"]
    assert len(module.sections) == 11
    assert {section.id for section in module.sections} == {
        "vcf-sv-stats-overview",
        "vcf-sv-stats-types",
        "vcf-sv-stats-length",
        "vcf-sv-stats-filters",
        "vcf-sv-stats-breakends",
        "vcf-sv-stats-merged-support",
        "vcf-sv-stats-genotypes",
        "vcf-sv-stats-copy-number",
        "vcf-sv-stats-validation",
        "vcf-sv-stats-normalization",
        "vcf-sv-stats-provenance",
    }
    assert module.vcf_sv_stats_data["vss1-test-0000"]["report"]["mapped_vcf_sample_ids"] == ["HG002"]


def test_module_deduplicates_identical_and_rejects_conflicting_reports(
    tmp_path: Path,
) -> None:
    first = make_summary()
    _write_summary(tmp_path / "first.vcf-sv-stats.json", first)
    _write_summary(tmp_path / "second.vcf-sv-stats.json", first)
    report.analysis_files = [str(tmp_path)]
    report.search_files(["vcf_sv_stats"])
    module = MultiqcModule()
    assert len(module.vcf_sv_stats_data) == 1
    assert module.duplicate_sources == ["second.vcf-sv-stats.json"]

    conflict = make_summary()
    conflict["reports"][0]["statistics"]["events"]["resolved"] = 12
    resign(conflict)
    _write_summary(tmp_path / "third.vcf-sv-stats.json", conflict)
    report.reset()
    report.analysis_files = [str(tmp_path)]
    report.search_files(["vcf_sv_stats"])
    with pytest.raises(SummaryValidationError, match="Conflicting"):
        MultiqcModule()
