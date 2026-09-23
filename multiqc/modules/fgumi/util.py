"""Helpers shared by the fgumi parsers."""

import logging
import math
import re
from typing import Any, Dict, Iterator, List, Mapping, Optional, Type, TypeVar

from pydantic import ValidationError

from multiqc.base_module import BaseMultiqcModule

from .schemas import FgumiMetric, M, MetricFormatError

log = logging.getLogger(__name__)

# fgumi's fixed output suffixes, longest first so e.g. ".duplex_family_sizes.txt" wins over ".family_sizes.txt".
# MultiQC's default name cleaning would otherwise leave "S1.duplex_family_sizes" as the sample name.
FGUMI_SUFFIXES = (
    ".simplex_yield_metrics.txt",
    ".duplex_yield_metrics.txt",
    ".duplex_family_sizes.txt",
    ".position_group_sizes.txt",
    ".duplex_umi_counts.txt",
    ".grouping_metrics.txt",
    ".family_sizes.txt",
    ".umi_counts.txt",
)


def sample_name(module: BaseMultiqcModule, f: Any) -> str:
    """The sample name for fgumi file ``f``: its fgumi suffix stripped, then MultiQC's usual cleaning."""
    fn = f["fn"][:-3] if f["fn"].endswith(".gz") else f["fn"]
    for suffix in FGUMI_SUFFIXES:
        if fn.endswith(suffix):
            return module.clean_s_name(fn[: -len(suffix)], f)
    return module.clean_s_name(fn, f)


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


def safe_id(text: str) -> str:
    """``text`` reduced to characters that are safe in a plot or section id."""
    return re.sub(r"[^A-Za-z0-9_-]+", "_", text)


def load_rows(f: Any, schema: Type[M]) -> Optional[List[M]]:
    """All rows of MultiQC file ``f`` validated as ``schema``; ``None`` (with a warning) if malformed."""
    if f["f"] is None:
        return None
    try:
        return schema.read(f["f"], f["fn"])
    except MetricFormatError as error:
        log.warning(f"Skipping {error}")
        return None


def stream_dicts(f: Any, schema: Type[FgumiMetric], validate_first: int = 100) -> Iterator[Dict[str, str]]:
    """Raw ``{column: cell}`` rows of a large file ``f`` (found with ``filehandles=True``).

    Only the header and the first ``validate_first`` rows are validated as ``schema``, so a file with a million
    per-UMI rows is not turned into a million pydantic objects. Raises ``MetricFormatError`` if malformed.
    """
    for number, row in enumerate(schema.iter_dicts(f["f"], f["fn"]), start=1):
        if number <= validate_first:
            try:
                schema.model_validate(row)
            except ValidationError as error:
                raise MetricFormatError(f"{f['fn']}: data row {number}: {error}") from error
        yield row
