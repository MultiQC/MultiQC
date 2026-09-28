"""Typed row schemas for fgumi metric files.

``FgumiMetric`` mirrors the API of fgmetric (https://github.com/fg-labs/fgmetric): one pydantic model per
row of a delimited metric file. fgmetric is not a dependency because it requires Python >= 3.12 and
pydantic >= 2.11, while MultiQC supports Python >= 3.9 and pydantic >= 2.7. Keeping the same shape makes a
later switch mechanical.

Column names and order follow fgumi's published contract, ``crates/fgumi-metrics/metric_columns.json``.
"""

from typing import (
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Type,
    TypeVar,
    Union,
    get_args,
    get_origin,
)

from pydantic import BaseModel, ConfigDict, ValidationError, model_validator

M = TypeVar("M", bound="FgumiMetric")


class MetricFormatError(ValueError):
    """A metric file whose header or rows do not match the expected schema."""


def _allows_none(annotation: object) -> bool:
    return get_origin(annotation) is Union and type(None) in get_args(annotation)


class FgumiMetric(BaseModel):
    """One row of an fgbio-style metric TSV: a header row of field names, then one row per record.

    Empty cells become ``None`` for ``Optional`` fields. Float cells accept fgbio's non-finite tokens
    (``NaN``, ``Infinity``, ``-Infinity``). Columns not declared on the schema are ignored.
    """

    model_config = ConfigDict(extra="ignore")

    @model_validator(mode="before")
    @classmethod
    def _coerce_cells(cls, data: object) -> object:
        if not isinstance(data, dict):
            return data
        data = dict(data)
        for name, field in cls.model_fields.items():
            key = field.alias or name
            value = data.get(key)
            if not isinstance(value, str):
                continue
            if value == "" and _allows_none(field.annotation):
                data[key] = None
        # pydantic itself parses fgbio's "NaN" / "Infinity" / "-Infinity" cells into float fields.
        return data

    @classmethod
    def columns(cls) -> List[str]:
        """The column names this schema reads, in declaration order."""
        return [field.alias or name for name, field in cls.model_fields.items()]

    @classmethod
    def matches(cls, header: Iterable[str]) -> bool:
        """Whether a file with column names ``header`` has every column this schema reads."""
        return set(cls.columns()) <= set(header)

    @classmethod
    def validate_row(cls: Type[M], row: Dict[str, str], source: str, number: int) -> M:
        """Data row ``number`` of ``source`` validated as this schema; ``MetricFormatError`` if it does not fit."""
        try:
            return cls.model_validate(row)
        except ValidationError as error:
            raise MetricFormatError(f"{source}: data row {number}: {error}") from error

    @classmethod
    def iter_dicts(cls, lines: Iterable[str], source: str) -> Iterator[Dict[str, str]]:
        """Yields each data row as ``{column: cell}`` after checking the header has every schema column.

        Raises ``MetricFormatError`` naming ``source`` for an empty file, a missing column, or a row whose
        cell count differs from the header's.
        """
        header: Optional[List[str]] = None
        for number, raw in enumerate(lines, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            cells = line.split("\t")
            if header is None:
                header = cells
                missing = [column for column in cls.columns() if column not in header]
                if missing:
                    raise MetricFormatError(f"{source}: header is missing column(s) {missing}")
                continue
            if len(cells) != len(header):
                raise MetricFormatError(f"{source}: line {number} has {len(cells)} cells, header has {len(header)}")
            yield dict(zip(header, cells))
        if header is None:
            raise MetricFormatError(f"{source}: empty file, expected a header row")

    @classmethod
    def iter_rows(cls: Type[M], lines: Iterable[str], source: str) -> Iterator[M]:
        """Yields each data row validated as this schema."""
        for number, row in enumerate(cls.iter_dicts(lines, source), start=1):
            yield cls.validate_row(row, source, number)

    @classmethod
    def read(cls: Type[M], text: str, source: str) -> List[M]:
        """All rows of ``text`` (a whole metric file) validated as this schema."""
        return list(cls.iter_rows(text.splitlines(), source))


# --- group / dedup ------------------------------------------------------------------------------------------


class FamilySizeMetric(FgumiMetric):
    """`group.family_sizes`, `dedup.family_sizes` (same columns as fgbio GroupReadsByUmi's histogram)."""

    family_size: int
    count: int
    fraction: float
    fraction_gt_or_eq_family_size: float


class PositionGroupSizeMetric(FgumiMetric):
    """`group.position_group_sizes`."""

    position_group_size: int
    count: int
    fraction: float
    fraction_gt_or_eq_position_group_size: float


class UmiGroupingMetric(FgumiMetric):
    """`group.grouping_metrics` (fgbio UmiGroupingMetric columns, including the `umis_to_short` spelling)."""

    accepted_sam_records: int
    discarded_non_pf: int
    discarded_poor_alignment: int
    discarded_ns_in_umi: int
    discarded_umis_to_short: int


class DeduplicationMetric(FgumiMetric):
    """`dedup.metrics`: one row per library plus a final `All Reads` row."""

    sample: str
    library: str
    filtered_templates: int
    total_templates: int
    unique_templates: int
    duplicate_templates: int
    duplicate_rate: float
    total_reads: int
    percent_duplication: float
    estimated_library_size: Optional[int]


class DuplicationLadderMetric(FgumiMetric):
    """`dedup.duplication_ladder`."""

    library: str
    templates_seen: int
    duplicate_fraction: float
    window_templates: int
    window_duplicate_fraction: float


# --- simplex / duplex metrics ---------------------------------------------------------------------------------


class SimplexFamilySizeMetric(FgumiMetric):
    """`simplex.family_sizes`."""

    family_size: int
    cs_count: int
    cs_fraction: float
    cs_fraction_gt_or_eq_size: float
    ss_count: int
    ss_fraction: float
    ss_fraction_gt_or_eq_size: float


class DuplexStrandFamilySizeMetric(SimplexFamilySizeMetric):
    """`duplex.family_sizes` (the simplex columns plus double-strand families)."""

    ds_count: int
    ds_fraction: float
    ds_fraction_gt_or_eq_size: float


class DuplexFamilySizeMetric(FgumiMetric):
    """`duplex.duplex_family_sizes`: counts by (AB strand size, BA strand size)."""

    ab_size: int
    ba_size: int
    count: int
    fraction: float
    fraction_gt_or_eq_size: float


class SimplexYieldMetric(FgumiMetric):
    """`simplex.simplex_yield_metrics`."""

    fraction: float
    read_pairs: int
    cs_families: int
    ss_families: int
    mean_ss_family_size: float
    ss_singletons: int
    ss_singleton_fraction: float
    ss_consensus_families: int


class DuplexYieldMetric(FgumiMetric):
    """`duplex.duplex_yield_metrics`."""

    fraction: float
    read_pairs: int
    cs_families: int
    ss_families: int
    ds_families: int
    ds_duplexes: int
    ds_fraction_duplexes: float
    ds_fraction_duplexes_ideal: float


class UmiMetric(FgumiMetric):
    """`simplex.umi_counts`, `duplex.umi_counts`: one row per UMI sequence."""

    umi: str
    raw_observations: int
    raw_observations_with_errors: int
    unique_observations: int
    fraction_raw_observations: float
    fraction_unique_observations: float


class DuplexUmiMetric(UmiMetric):
    """`duplex.duplex_umi_counts`: one row per duplex UMI pair."""

    fraction_unique_observations_expected: float


class ConsensusStatMetric(FgumiMetric):
    """`consensus.stats`: the `--stats` key/value/description rows of `simplex`, `duplex` and `codec`."""

    key: str
    value: str
    description: str


# --- other commands -------------------------------------------------------------------------------------------


class UmiCorrectionMetric(FgumiMetric):
    """`correct.metrics`: one row per expected UMI, plus an all-`N` row for unmatched UMIs."""

    umi: str
    total_matches: int
    perfect_matches: int
    one_mismatch_matches: int
    two_mismatch_matches: int
    other_matches: int
    fraction_of_matches: float
    representation: float


class ClippingMetric(FgumiMetric):
    """`clip.metrics`: rows for Fragment, ReadOne, ReadTwo, Pair and All."""

    read_type: str
    reads: int
    reads_clipped_five_prime: int
    reads_clipped_three_prime: int
    reads_clipped_overlapping: int
    reads_clipped_extending: int
    bases: int
    bases_clipped_five_prime: int
    bases_clipped_three_prime: int
    bases_clipped_overlapping: int
    bases_clipped_extending: int


class FilterStatsMetric(FgumiMetric):
    """`filter.stats`."""

    total_reads: int
    passed_reads: int
    failed_reads: int
    pass_rate: float


class CopyUmiMetric(FgumiMetric):
    """`copy_umi.metrics`."""

    total_records: int
    rx_written: int
    rx_overwritten: int
    names_trimmed: int


class RetagMetric(FgumiMetric):
    """`retag.metrics`: one row per tag operation."""

    operation: str
    kind: str
    records_applied: int
    dst_overwritten: int
    src_missing: int


class DownsampleHistogramMetric(FgumiMetric):
    """`downsample.histogram_kept`, `downsample.histogram_rejected`."""

    family_size: int
    count: int


class ReviewDetailMetric(FgumiMetric):
    """`review.details`: one row per (variant site, consensus read) observation."""

    chrom: str
    pos: int
    ref: str
    consensus_read: str
