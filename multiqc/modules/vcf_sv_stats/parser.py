import hashlib
import json
from dataclasses import dataclass
from typing import Any

import jsonschema
import rfc8785

SUMMARY_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "urn:vcf-sv-stats:schema:summary:1.0.0",
    "type": "object",
    "additionalProperties": False,
    "required": [
        "schema_name",
        "schema_version",
        "content_signature",
        "producer",
        "input",
        "callset",
        "validation",
        "statistics",
        "reports",
        "payload_sha256",
    ],
    "properties": {
        "schema_name": {"const": "vcf-sv-stats.summary"},
        "schema_version": {"const": "1.0.0"},
        "content_signature": {"const": "vcf-sv-stats:summary:1"},
        "producer": {
            "type": "object",
            "additionalProperties": False,
            "required": ["name", "version"],
            "properties": {
                "name": {"const": "vcf-sv-stats"},
                "version": {"type": "string", "minLength": 1},
            },
        },
        "input": {"type": "object"},
        "callset": {"type": "object"},
        "validation": {"type": "object"},
        "statistics": {"type": "object"},
        "reports": {"type": "array", "items": {"type": "object"}},
        "payload_sha256": {"type": "string", "pattern": "^[0-9a-f]{64}$"},
        "execution": {"type": ["object", "null"]},
    },
}


class SummaryValidationError(ValueError):
    pass


@dataclass(frozen=True)
class ParsedReport:
    report_id: str
    payload_sha256: str
    producer_version: str
    callset: dict[str, Any]
    validation: dict[str, Any]
    report: dict[str, Any]


def mapping(value: Any, field: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise SummaryValidationError(f"{field} must be an object")
    return value


def non_negative_integer(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise SummaryValidationError(f"{field} must be a non-negative integer")
    return value


def non_empty_string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise SummaryValidationError(f"{field} must be a non-empty string")
    return value


def count_mapping(value: Any, field: str) -> dict[str, int]:
    raw = mapping(value, field)
    return {str(key): non_negative_integer(count, f"{field}.{key}") for key, count in raw.items()}


def _validate_statistics(statistics: dict[str, Any], prefix: str) -> None:
    source_records = mapping(statistics["source_records"], f"{prefix}.source_records")
    non_negative_integer(source_records["total"], f"{prefix}.source_records.total")
    alleles = mapping(statistics["alleles"], f"{prefix}.alleles")
    non_negative_integer(alleles["total"], f"{prefix}.alleles.total")
    count_mapping(alleles["types"], f"{prefix}.alleles.types")
    events = mapping(statistics["events"], f"{prefix}.events")
    non_negative_integer(events["resolved"], f"{prefix}.events.resolved")
    breakends = mapping(statistics["breakends"], f"{prefix}.breakends")
    for key in (
        "total",
        "reciprocal_pairs",
        "without_declared_mate",
        "unresolved_mate_references",
    ):
        non_negative_integer(breakends[key], f"{prefix}.breakends.{key}")
    count_mapping(statistics["filters"], f"{prefix}.filters")
    count_mapping(statistics["genotypes"], f"{prefix}.genotypes")
    count_mapping(statistics["copy_number"], f"{prefix}.copy_number")

    support = mapping(statistics["merged_support"], f"{prefix}.merged_support")
    non_empty_string(support["status"], f"{prefix}.merged_support.status")
    declared_fields = support["declared_fields"]
    if not isinstance(declared_fields, list) or not all(isinstance(field, str) and field for field in declared_fields):
        raise SummaryValidationError(f"{prefix}.merged_support.declared_fields must be a string array")
    records_with_support = non_negative_integer(
        support["records_with_support"], f"{prefix}.merged_support.records_with_support"
    )
    non_negative_integer(support["records_without_support"], f"{prefix}.merged_support.records_without_support")
    supporting_sources = count_mapping(support["supporting_sources"], f"{prefix}.merged_support.supporting_sources")
    if sum(supporting_sources.values()) != records_with_support:
        raise SummaryValidationError(f"{prefix}.merged_support counts do not reconcile")
    count_mapping(support["source_count_states"], f"{prefix}.merged_support.source_count_states")
    consistency = count_mapping(
        support["vector_count_consistency"],
        f"{prefix}.merged_support.vector_count_consistency",
    )
    if sum(consistency.values()) != records_with_support:
        raise SummaryValidationError(f"{prefix}.merged_support consistency does not reconcile")
    if support["interpretation"] != "merger_provenance_only":
        raise SummaryValidationError(f"{prefix}.merged_support interpretation is unsupported")

    length = mapping(statistics["length_bp"], f"{prefix}.length_bp")
    boundaries = length["boundaries"]
    counts = length["counts"]
    if not isinstance(boundaries, list) or not isinstance(counts, list):
        raise SummaryValidationError(f"{prefix}.length_bp histogram must use arrays")
    if len(boundaries) != len(counts):
        raise SummaryValidationError(f"{prefix}.length_bp histogram cardinality differs")
    total = sum(non_negative_integer(count, f"{prefix}.length_bp.counts") for count in counts)
    if total != non_negative_integer(length["n"], f"{prefix}.length_bp.n"):
        raise SummaryValidationError(f"{prefix}.length_bp counts do not reconcile")

    contracts = mapping(statistics["metric_contracts"], f"{prefix}.metric_contracts")
    for metric in (
        "alleles",
        "breakends",
        "copy_number",
        "events",
        "filters",
        "genotypes",
        "length_bp",
        "merged_support",
        "source_records",
    ):
        contract = mapping(contracts[metric], f"{prefix}.metric_contracts.{metric}")
        for field in ("scope", "denominator", "unit", "comparability"):
            non_empty_string(contract[field], f"{prefix}.metric_contracts.{metric}.{field}")


def _parse_summary(content: str, source_name: str) -> tuple[ParsedReport, ...]:
    try:
        payload = json.loads(content)
    except json.JSONDecodeError as exc:
        raise SummaryValidationError(f"Invalid JSON in {source_name}") from exc
    if not isinstance(payload, dict):
        raise SummaryValidationError(f"Summary root must be an object in {source_name}")

    version = str(payload.get("schema_version", ""))
    try:
        major = int(version.split(".", 1)[0])
    except ValueError as exc:
        raise SummaryValidationError(f"Invalid schema version in {source_name}") from exc
    if major != 1:
        raise SummaryValidationError(f"Unsupported summary schema major in {source_name}")
    if version != "1.0.0":
        raise SummaryValidationError(f"Unsupported summary schema version in {source_name}")
    try:
        jsonschema.Draft202012Validator(SUMMARY_SCHEMA).validate(payload)
    except jsonschema.ValidationError as exc:
        raise SummaryValidationError(f"Summary schema validation failed in {source_name}") from exc

    expected_digest = payload["payload_sha256"]
    digest_payload = dict(payload)
    digest_payload.pop("payload_sha256")
    digest_payload.pop("execution", None)
    observed_digest = hashlib.sha256(rfc8785.dumps(digest_payload)).hexdigest()
    if expected_digest != observed_digest:
        raise SummaryValidationError(f"Summary payload digest does not match in {source_name}")

    producer = mapping(payload["producer"], "producer")
    producer_version = non_empty_string(producer["version"], "producer.version")
    callset = mapping(payload["callset"], "callset")
    callset_producer = mapping(callset["producer"], "callset.producer")
    for field in ("adapter_id", "producer", "status"):
        non_empty_string(callset_producer[field], f"callset.producer.{field}")
    non_negative_integer(callset["record_count"], "callset.record_count")
    non_negative_integer(callset["allele_count"], "callset.allele_count")

    validation = mapping(payload["validation"], "validation")
    count_mapping(validation["diagnostic_counts"], "validation.diagnostic_counts")
    states = mapping(validation["states"], "validation.states")
    for field in (
        "container_state",
        "parse_state",
        "vcf_conformance_state",
        "sv_semantic_state",
        "operation_safety_state",
        "statistics_state",
    ):
        non_empty_string(states[field], f"validation.states.{field}")

    reports = payload["reports"]
    if not reports:
        raise SummaryValidationError(f"Summary contains no reports in {source_name}")
    parsed: list[ParsedReport] = []
    seen_ids: set[str] = set()
    for index, raw_report in enumerate(reports):
        report = mapping(raw_report, f"reports[{index}]")
        report_id = non_empty_string(report["report_id"], f"reports[{index}].report_id")
        if report_id in seen_ids:
            raise SummaryValidationError(f"Duplicate report identifier in {source_name}")
        seen_ids.add(report_id)
        analysis_unit = mapping(report["analysis_unit"], f"reports[{index}].analysis_unit")
        non_empty_string(analysis_unit["status"], f"reports[{index}].analysis_unit.status")
        mapped_samples = report["mapped_vcf_sample_ids"]
        if not isinstance(mapped_samples, list) or not all(
            isinstance(sample, str) and sample for sample in mapped_samples
        ):
            raise SummaryValidationError(f"reports[{index}].mapped_vcf_sample_ids must be a string array")
        statistics = mapping(report["statistics"], f"reports[{index}].statistics")
        _validate_statistics(statistics, f"reports[{index}].statistics")
        report_digest = hashlib.sha256(rfc8785.dumps(report)).hexdigest()
        parsed.append(
            ParsedReport(
                report_id=report_id,
                payload_sha256=report_digest,
                producer_version=producer_version,
                callset=callset,
                validation=validation,
                report=report,
            )
        )
    return tuple(parsed)


def parse_summary(content: str, source_name: str = "summary") -> tuple[ParsedReport, ...]:
    try:
        return _parse_summary(content, source_name)
    except SummaryValidationError:
        raise
    except (KeyError, TypeError) as exc:
        raise SummaryValidationError(f"Summary is missing a required nested field in {source_name}") from exc


def length_label(boundary: Any) -> str:
    if not isinstance(boundary, list) or len(boundary) != 2:
        raise SummaryValidationError("Length boundary must contain lower and upper values")
    lower, upper = boundary
    if not isinstance(lower, int) or (upper is not None and not isinstance(upper, int)):
        raise SummaryValidationError("Length boundary values must be integers or null")
    return f"{lower:,}+ bp" if upper is None else f"{lower:,}-{upper - 1:,} bp"
