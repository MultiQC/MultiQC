"""Pins every fgumi schema to fgumi's published column contract (``crates/fgumi-metrics/metric_columns.json``).

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
