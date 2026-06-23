"""Parse BibTeX (.bib) citations into the shared Citation model.

BibTeX is an optional alternative input to CSL-JSON. It requires the
`bibtexparser` package, which is not a core MultiQC dependency; install it with
`pip install multiqc[citations]`. If a `.bib` file is found without
`bibtexparser` installed, the file is skipped with a warning (CSL-JSON still
works).

BibTeX has no standard field for the tool name or the runtime version, so:

- tool name comes from a `tool` field if present, else the entry citation key;
- runtime version comes from a non-standard `version` field if present.
"""

import logging
import re
from typing import List, Optional, Tuple

from .citation import Authors, Citation, clean_doi

log = logging.getLogger(__name__)

_AND_SPLIT = re.compile(r"\s+and\s+", flags=re.IGNORECASE)


def bibtexparser_available() -> bool:
    """Whether the optional `bibtexparser` dependency can be imported."""
    try:
        import bibtexparser  # noqa: F401

        return True
    except ImportError:
        return False


def _debrace(value: Optional[str]) -> Optional[str]:
    """Strip BibTeX grouping braces and collapse whitespace."""
    if value is None:
        return None
    return re.sub(r"\s+", " ", value.replace("{", "").replace("}", "")).strip() or None


def _split_bibtex_name(name: str) -> Tuple[str, str]:
    """Split a single BibTeX name into (surname, given), per BibTeX semantics.

    "Andrews, S" -> ("Andrews", "S"); "Heng Li" -> ("Li", "Heng").
    """
    if "," in name:
        last, _, given = name.partition(",")
        return last.strip(), given.strip()
    parts = name.split()
    if len(parts) <= 1:
        return name.strip(), ""
    return parts[-1], " ".join(parts[:-1])


def _normalize_bibtex_authors(field: Optional[str]) -> Authors:
    """Normalise a BibTeX `author` value to an `Authors`."""
    field = _debrace(field)
    if not field:
        return Authors.from_names([], False)

    raw_names = [n.strip() for n in _AND_SPLIT.split(field) if n.strip()]
    has_etal = any(n.lower() == "others" for n in raw_names)  # BibTeX "and others"
    raw_names = [n for n in raw_names if n.lower() != "others"]

    names: List[Tuple[str, str]] = []
    for name in raw_names:
        last, given = _split_bibtex_name(name)
        names.append((last, f"{last} {given}".strip() if given else last))

    return Authors.from_names(names, has_etal)


def _entry_to_citation(entry: dict) -> Citation:
    tool = _debrace(entry.get("tool")) or entry.get("ID")
    if not tool:
        raise ValueError(f"BibTeX entry is missing both a `tool` field and a citation key: {entry!r}")

    year: Optional[int] = None
    raw_year = _debrace(entry.get("year"))
    if raw_year:
        match = re.search(r"\d{4}", raw_year)
        if match:
            year = int(match.group())

    return Citation(
        tool=str(tool),
        version=_debrace(entry.get("version")),
        title=_debrace(entry.get("title")),
        authors=_normalize_bibtex_authors(entry.get("author")),
        year=year,
        container_title=_debrace(entry.get("journal")) or _debrace(entry.get("booktitle")),
        doi=clean_doi(_debrace(entry.get("doi"))),
        url=_debrace(entry.get("url")) or _debrace(entry.get("howpublished")),
        csl_type=entry.get("ENTRYTYPE"),
    )


def parse_bibtex(content: str, fn: str = "<citations>") -> List[Citation]:
    """Parse BibTeX file content into Citations.

    Returns an empty list (with a warning) if `bibtexparser` is not installed.
    Raises ValueError with the file path on a parse failure.
    """
    if not bibtexparser_available():
        log.warning(
            "Found BibTeX citations file '%s' but the 'bibtexparser' package is not installed; "
            "skipping. Install it with: pip install multiqc[citations]",
            fn,
        )
        return []

    import bibtexparser
    from bibtexparser.bparser import BibTexParser

    try:
        parser = BibTexParser(common_strings=True, ignore_nonstandard_types=False)
        database = bibtexparser.loads(content, parser=parser)
    except Exception as exc:  # bibtexparser raises a variety of exception types
        raise ValueError(f"Could not parse BibTeX citations file '{fn}': {exc}") from exc

    return [_entry_to_citation(entry) for entry in database.entries]
