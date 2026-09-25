"""Helpers shared by the fgumi parsers."""

import logging
import math
from pathlib import PurePath
from typing import Any, Callable, Dict, Iterable, Iterator, List, Mapping, Optional, Tuple, Type, TypeVar

from multiqc import config, report
from multiqc.base_module import BaseMultiqcModule
from multiqc.types import ModuleId

from .schemas import FgumiMetric, M, MetricFormatError

log = logging.getLogger(__name__)

# fgumi's fixed output suffixes. MultiQC's default name cleaning would otherwise leave e.g.
# "S1.duplex_family_sizes" as the sample name. Checked in order, so a suffix that ends with another
# (".duplex_umi_counts.txt" / ".umi_counts.txt") must come first.
FGUMI_SUFFIXES = (
    ".simplex_yield_metrics.txt",
    ".duplex_yield_metrics.txt",
    ".duplex_family_sizes.txt",
    ".position_group_sizes.txt",
    ".duplex_umi_counts.txt",
    ".grouping_metrics.txt",
    ".family_sizes.txt",
    ".umi_counts.txt",
    # `fgumi runall --all-metrics <prefix>` names its non-prefix outputs `<prefix>.<stage>.<suffix>`.
    ".correct.metrics.txt",
    ".filter.stats.txt",
)
# The stage `fgumi runall` inserts before a `--metrics <prefix>` suffix, e.g. `S1.duplex.umi_counts.txt`.
RUNALL_STAGES = (".group", ".simplex", ".duplex", ".codec")


def sample_name(module: BaseMultiqcModule, f: Any) -> str:
    """The sample name for fgumi file ``f``: its fgumi suffix (and any `runall` stage before it) stripped, then
    MultiQC's usual cleaning. Nothing is stripped when sample name cleaning is off (``--fullnames``)."""
    fn = f["fn"]
    if config.fn_clean_sample_names:
        for suffix in FGUMI_SUFFIXES:
            if fn.endswith(suffix):
                fn = fn[: -len(suffix)]
                stage = next((stage for stage in RUNALL_STAGES if fn.endswith(stage)), None)
                fn = fn[: -len(stage)] if stage else fn
                break
    return module.clean_s_name(fn, f)


def register(module: BaseMultiqcModule, f: Any, s_name: Optional[str] = None, section: Optional[str] = None) -> str:
    """Records ``f`` as a source for its sample (``s_name``, else derived from the file name) and returns the name.

    The data source is keyed by ``section`` (default: the file's search pattern), so a sample's several fgumi files
    are all listed in ``multiqc_sources.txt`` rather than each replacing the last. A second file for the same sample
    and section replaces the first in the report, which is logged.
    """
    s_name = sample_name(module, f) if s_name is None else s_name
    section = f["sp_key"].split("/", 1)[-1] if section is None else section
    if s_name in report.data_sources[module.name][section]:
        log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name} ({section})")
    module.add_data_source(f, s_name, section=section)
    module.add_software_version(None, s_name)
    return s_name


def flatten(data: Mapping[str, Mapping[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """``data`` as the two-level ``{sample: {column: value}}`` table MultiQC's data files expect.

    A nested ``{tab: {x: value}}`` cell becomes ``{"<tab>_<x>": value}`` columns, and ``None`` values are left
    out so they are written as empty cells rather than the text "None".
    """
    table: Dict[str, Dict[str, Any]] = {}
    for s_name, row in data.items():
        flat: Dict[str, Any] = {}
        for column, value in row.items():
            if isinstance(value, Mapping):
                flat.update({f"{column}_{x}": v for x, v in value.items() if v is not None})
            elif value is not None:
                flat[column] = value
        table[s_name] = flat
    return table


# Search patterns only fgumi can match: their files show the run has fgumi output. The consensus --stats pattern
# is left out because its generic key/value/description header can match other tools' files.
FGUMI_ONLY_PATTERNS = (
    "fgumi/position_group_sizes",
    "fgumi/dedup",
    "fgumi/dedup_ladder",
    "fgumi/filter_stats",
    "fgumi/copy_umi",
    "fgumi/retag",
    "fgumi/downsample_histogram",
)
# Search patterns only fgbio can match.
FGBIO_ONLY_PATTERNS = ("fgbio/errorratebyreadposition",)
FAMILY_SIZE_MODULES = ("fgumi", "fgbio")


def _directory_and_ancestors(root: str) -> List[str]:
    path = PurePath(root)
    return [str(path)] + [str(parent) for parent in path.parents]


def family_sizes_evidence() -> Dict[str, Tuple[bool, bool]]:
    """For every directory holding (at any depth) fgumi-only or fgbio-only files: whether it has each kind."""
    evidence: Dict[str, Tuple[bool, bool]] = {}
    for keys, is_fgumi in [(FGUMI_ONLY_PATTERNS, True), (FGBIO_ONLY_PATTERNS, False)]:
        for key in keys:
            for g in report.files.get(ModuleId(key), []):
                for directory in _directory_and_ancestors(g["root"]):
                    has_fgumi, has_fgbio = evidence.get(directory, (False, False))
                    evidence[directory] = (has_fgumi or is_fgumi, has_fgbio or not is_fgumi)
    return evidence


def family_sizes_module(f: Any, evidence: Dict[str, Tuple[bool, bool]]) -> str:
    """The module ("fgumi" or "fgbio") that reports family-size histogram ``f`` when both modules find it.

    fgumi and fgbio GroupReadsByUmi write byte-identical histograms, so the choice comes from the files around it.
    The nearest directory (the file's own, then each parent) holding fgumi-only or fgbio-only files decides: fgumi
    if it has fgumi-only files and no fgbio-only files, fgbio otherwise. With no such directory it is fgbio, which
    reported these files before the fgumi module existed. ``fgumi_config: {family_sizes_module: fgumi|fgbio}``
    overrides this for every file. ``evidence`` comes from ``family_sizes_evidence()``.
    """
    choice = (getattr(config, "fgumi_config", None) or {}).get("family_sizes_module")
    if choice in FAMILY_SIZE_MODULES:
        return choice
    if choice is not None:
        log.warning(f"Ignoring fgumi_config.family_sizes_module {choice!r}: expected one of {FAMILY_SIZE_MODULES}")
    for directory in _directory_and_ancestors(f["root"]):
        if directory in evidence:
            has_fgumi, has_fgbio = evidence[directory]
            return "fgumi" if has_fgumi and not has_fgbio else "fgbio"
    return "fgbio"


def found_by(f: Any, sp_key: str) -> bool:
    """Whether MultiQC file ``f`` was also found by search pattern ``sp_key`` (so that module will read it)."""
    return any(g["root"] == f["root"] and g["fn"] == f["fn"] for g in report.files.get(ModuleId(sp_key), []))


def finite(value: Optional[float]) -> Optional[float]:
    """``value``, or ``None`` when it is missing or non-finite (MultiQC plots cannot show NaN/Infinity)."""
    return value if value is not None and math.isfinite(value) else None


def pct(fraction: Optional[float]) -> Optional[float]:
    """A fraction as a finite percentage, or ``None``."""
    value = finite(fraction)
    return None if value is None else value * 100.0


K = TypeVar("K")
V = TypeVar("V")


def drop_none(mapping: Mapping[K, Optional[V]]) -> Dict[K, V]:
    """``mapping`` without its ``None`` values: a missing point or cell, rather than a ``None`` MultiQC would
    have to interpret."""
    return {key: value for key, value in mapping.items() if value is not None}


def header_columns(f: Any) -> List[str]:
    """The column names on the first line of MultiQC file ``f`` (loaded as text)."""
    return f["f"].split("\n", 1)[0].rstrip("\r").split("\t")


def load_rows(f: Any, schema: Type[M]) -> Optional[List[M]]:
    """All rows of MultiQC file ``f`` validated as ``schema``; ``None`` (with a warning) if malformed."""
    if f["f"] is None:
        return None
    try:
        return schema.read(f["f"], f["fn"])
    except MetricFormatError as error:
        log.warning(f"Skipping {error}")
        return None


def iter_samples(
    module: BaseMultiqcModule, sp_key: str, schema: Type[M], skip: Optional[Callable[[Any], bool]] = None
) -> Iterator[Tuple[str, List[M]]]:
    """``(sample name, rows)`` for each well-formed file found by ``sp_key`` whose sample is not ignored.

    Each file yielded is registered as a data source. Files for which ``skip`` returns true are passed over
    without being registered. A malformed file is skipped with a warning; a file with a header but no rows is
    yielded with no rows.
    """
    for f in module.find_log_files(sp_key):
        if skip is not None and skip(f):
            continue
        rows = load_rows(f, schema)
        if rows is None:
            continue
        s_name = sample_name(module, f)
        if module.is_ignore_sample(s_name):
            continue
        yield register(module, f, s_name), rows


# Rows of a streamed file validated against its schema; later rows are only split into cells.
STREAM_VALIDATE_ROWS = 100


def stream_dicts(lines: Iterable[str], source: str, schema: Type[FgumiMetric]) -> Iterator[Dict[str, str]]:
    """Raw ``{column: cell}`` rows of a large file (found with ``filehandles=True``) named ``source``.

    Only the header and the first ``STREAM_VALIDATE_ROWS`` rows are validated as ``schema``, so a file with a
    million per-UMI rows is not turned into a million pydantic objects. Raises ``MetricFormatError`` if malformed.
    """
    for number, row in enumerate(schema.iter_dicts(lines, source), start=1):
        if number <= STREAM_VALIDATE_ROWS:
            schema.validate_row(row, source, number)
        yield row


def skip_unreadable(f: Any, error: Exception) -> None:
    """Warns that streamed file ``f`` is skipped because of ``error`` (a format, value or decoding error)."""
    # MetricFormatError messages already start with the file name.
    log.warning(f"Skipping {error}" if isinstance(error, MetricFormatError) else f"Skipping {f['fn']}: {error}")
