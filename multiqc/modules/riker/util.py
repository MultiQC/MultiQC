"""Shared helpers for the riker submodules."""

import logging
from typing import Dict, Iterator, List, TextIO

log = logging.getLogger(__name__)


def read_tsv(handle: TextIO, source: str = "<unknown>") -> Iterator[Dict[str, str]]:
    """
    Yield rows of a riker TSV output as ``{column: value}`` dicts.

    Riker outputs are plain TSVs with a header on the first line, no comment
    or metadata lines, and ``sample`` as the first column. Empty/blank rows
    are skipped. A row whose column count does not match the header is a
    sign that the file format has drifted or the file is corrupt; we log a
    warning naming ``source`` and the line number and skip the row, rather
    than silently dropping it.
    """
    header_line = handle.readline()
    if not header_line:
        return
    header: List[str] = header_line.rstrip("\n").split("\t")

    for line_num, line in enumerate(handle, start=2):
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            log.warning(
                f"riker: skipping malformed row in {source} line {line_num}: "
                f"got {len(fields)} fields, expected {len(header)}"
            )
            continue
        yield dict(zip(header, fields))


def to_int(value: str) -> int:
    """Parse an integer column. Falls back to int(float(...)) for values written as floats."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(float(value))
