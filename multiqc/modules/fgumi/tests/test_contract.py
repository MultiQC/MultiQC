"""Pins every fgumi schema to fgumi's column contract (``crates/fgumi-metrics/metric_columns.json``, vendored from
fulcrumgenomics/fgumi#979).

Update ``metric_columns.json`` here whenever fgumi's contract changes; a schema column that fgumi no longer
emits (or emits in a different order) then fails this test instead of silently parsing nothing.
"""

import json
from pathlib import Path

import pytest

from multiqc.modules.fgumi import schemas

MANIFEST = json.loads((Path(__file__).parent / "metric_columns.json").read_text())

SCHEMA_BY_MANIFEST_KEY = {
    "clip.metrics": schemas.ClippingMetric,
    "consensus.stats": schemas.ConsensusStatMetric,
    "copy_umi.metrics": schemas.CopyUmiMetric,
    "correct.metrics": schemas.UmiCorrectionMetric,
    "dedup.duplication_ladder": schemas.DuplicationLadderMetric,
    "dedup.family_sizes": schemas.FamilySizeMetric,
    "dedup.metrics": schemas.DeduplicationMetric,
    "downsample.histogram_kept": schemas.DownsampleHistogramMetric,
    "downsample.histogram_rejected": schemas.DownsampleHistogramMetric,
    "duplex.duplex_family_sizes": schemas.DuplexFamilySizeMetric,
    "duplex.duplex_umi_counts": schemas.DuplexUmiMetric,
    "duplex.duplex_yield_metrics": schemas.DuplexYieldMetric,
    "duplex.family_sizes": schemas.DuplexStrandFamilySizeMetric,
    "duplex.umi_counts": schemas.UmiMetric,
    "filter.stats": schemas.FilterStatsMetric,
    "group.family_sizes": schemas.FamilySizeMetric,
    "group.grouping_metrics": schemas.UmiGroupingMetric,
    "group.position_group_sizes": schemas.PositionGroupSizeMetric,
    "retag.metrics": schemas.RetagMetric,
    "review.details": schemas.ReviewDetailMetric,
    "simplex.family_sizes": schemas.SimplexFamilySizeMetric,
    "simplex.simplex_yield_metrics": schemas.SimplexYieldMetric,
    "simplex.umi_counts": schemas.UmiMetric,
}


def test_every_manifest_entry_has_a_schema():
    assert set(SCHEMA_BY_MANIFEST_KEY) == set(MANIFEST)


@pytest.mark.parametrize("key", sorted(SCHEMA_BY_MANIFEST_KEY))
def test_schema_columns_are_an_ordered_subset_of_the_manifest(key):
    manifest_columns = MANIFEST[key]
    schema_columns = SCHEMA_BY_MANIFEST_KEY[key].columns()
    missing = [c for c in schema_columns if c not in manifest_columns]
    assert not missing, f"{key}: schema reads column(s) fgumi does not emit: {missing}"
    positions = [manifest_columns.index(c) for c in schema_columns]
    assert positions == sorted(positions), f"{key}: schema column order differs from fgumi's"


@pytest.mark.parametrize("key", sorted(SCHEMA_BY_MANIFEST_KEY))
def test_schema_reads_a_manifest_header(key):
    # A header-only file with exactly fgumi's columns must parse (to zero rows) with the mapped schema.
    assert SCHEMA_BY_MANIFEST_KEY[key].read("\t".join(MANIFEST[key]) + "\n", key) == []


def _fgumi_search_patterns():
    import yaml

    patterns = yaml.safe_load((Path(__file__).parents[3] / "search_patterns.yaml").read_text())
    return {key: value for key, value in patterns.items() if key.startswith("fgumi/")}


def _matches(pattern, header_line):
    import re

    specs = pattern if isinstance(pattern, list) else [pattern]
    for spec in specs:
        if "contents" in spec and spec["contents"] in header_line:
            return True
        if "contents_re" in spec and re.search(spec["contents_re"], header_line):
            return True
    return False


@pytest.mark.parametrize("key", sorted(MANIFEST))
def test_each_manifest_header_matches_exactly_one_search_pattern(key):
    # Ties search_patterns.yaml to fgumi's contract: a renamed column used in a pattern would otherwise stop
    # detection silently while the schema tests stay green.
    header = "\t".join(MANIFEST[key])
    matching = [name for name, pattern in _fgumi_search_patterns().items() if _matches(pattern, header)]
    assert len(matching) == 1, f"{key}: header matched {matching}"


def test_module_is_in_module_order():
    from multiqc import config

    names = [next(iter(entry)) if isinstance(entry, dict) else entry for entry in config.module_order]
    assert "fgumi" in names


def test_docstring_names_the_fgbio_outputs_the_module_also_reads():
    from multiqc.modules.fgumi import MultiqcModule

    for tool in ("GroupReadsByUmi", "CollectDuplexSeqMetrics", "CorrectUmis", "ClipBam"):
        assert tool in MultiqcModule.__doc__, tool


def test_every_section_has_a_description_and_help_text(run_fgumi):
    # MultiQC's new-module checklist asks for both on every section.
    data_dir = Path(__file__).parents[5] / "MultiQC-test-data" / "data" / "modules" / "fgumi"
    if not data_dir.is_dir():
        pytest.skip("MultiQC test-data checkout not found next to the MultiQC repo")
    module = run_fgumi({p.name: p.read_text() for p in data_dir.iterdir() if p.suffix == ".txt"})
    missing = [s.name for s in module.sections if not s.description or not s.helptext]
    assert not missing, missing
