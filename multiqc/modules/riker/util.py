"""Shared helpers for the riker submodules."""

from typing import Dict, Iterator, List, TextIO


def read_tsv(handle: TextIO) -> Iterator[Dict[str, str]]:
    """
    Yield rows of a riker TSV output as ``{column: value}`` dicts.

    Riker outputs are plain TSVs with a header on the first line, no comment
    or metadata lines, and ``sample`` as the first column. Empty / blank rows
    and rows whose column count does not match the header are skipped.
    """
    header_line = handle.readline()
    if not header_line:
        return
    header: List[str] = header_line.rstrip("\n").split("\t")

    for line in handle:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            continue
        yield dict(zip(header, fields))


def to_float(value: str) -> float:
    """Parse a numeric string from riker output to float, returning ``nan`` on failure."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def to_int(value: str) -> int:
    """Parse an integer column. Falls back to int(float(...)) for values written as floats."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(float(value))
