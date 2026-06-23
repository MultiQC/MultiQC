"""The format-agnostic Citation model and Harvard short-form rendering.

`Citation` is the shared value object that both the CSL-JSON adapter (`csl.py`)
and the BibTeX adapter (`bibtex.py`) produce: each adapter reads its own input
format and emits `Citation` objects. This module also holds the rendering of
those citations into report HTML, a methods sentence and a bibliography, plus
the shared helpers both adapters reuse (`clean_doi`, author normalisation).

The Harvard short form is "Surname (year)" for one author and
"Surname et al. (year)" for several (italic "et al.").
"""

import html
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


def _esc(value: Any) -> str:
    """Escape text destined for HTML element content."""
    return html.escape(str(value), quote=False)


def _esc_attr(value: Any) -> str:
    """Escape text destined for an HTML attribute value (e.g. an href)."""
    return html.escape(str(value), quote=True)


def clean_doi(raw: Optional[str]) -> Optional[str]:
    """Normalise a DOI to its bare form (strip `doi:` and resolver prefixes)."""
    if not raw:
        return None
    doi = str(raw).strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix):
            doi = doi[len(prefix) :].strip()
    return doi or None


def _surname_from_name(name: str) -> str:
    """Extract a surname from a free-text personal name.

    Uses the "surname-first" convention common in bioinformatics citations:
    the first whitespace token before the first comma. "Andrews S" -> "Andrews";
    "Li H, Durbin R" -> "Li". This is a heuristic for CSL `literal` and
    plain-string authors; structured authors use their `family` field directly.
    """
    head = name.split(",", 1)[0].strip()
    tokens = head.split()
    return tokens[0] if tokens else head


@dataclass
class Authors:
    """The shared author shape: surnames, a display string, and is_multiple.

    Each format adapter parses its own raw author field into format-specific
    `(surname, display_piece)` pairs plus an "et al." flag; `from_names` owns the
    rest (the display join, the "et al." suffix, and the is_multiple rule).
    """

    surnames: List[str]
    display: str
    is_multiple: bool

    @classmethod
    def from_names(cls, names: List[Tuple[str, str]], has_etal: bool) -> "Authors":
        # names: list of (surname, display_piece) in document order
        surnames = [s for s, _ in names]
        parts = [d for _, d in names]
        if has_etal:
            parts.append("et al.")
        return cls(surnames=surnames, display=", ".join(parts), is_multiple=has_etal or len(names) > 1)


def _normalize_authors(author_field: Any) -> Authors:
    """Normalise a CSL `author` value to an `Authors`.

    Accepts the CSL-standard list of `{family, given}` / `{literal}` dicts, and
    also a plain author string (e.g. "Ewels P, Magnusson M, Lundin S, et al.").
    """
    if not author_field:
        return Authors.from_names([], False)

    if isinstance(author_field, str):
        return _authors_from_string(author_field)

    if isinstance(author_field, list):
        # A single literal element may itself hold several comma-joined names.
        if len(author_field) == 1 and "literal" in author_field[0]:
            return _authors_from_string(str(author_field[0]["literal"]))

        names: List[Tuple[str, str]] = []
        for person in author_field:
            family = person.get("family")
            if family:
                given = person.get("given")
                names.append((str(family), f"{family} {given}".strip() if given else str(family)))
            elif person.get("literal"):
                literal = str(person["literal"])
                names.append((_surname_from_name(literal), literal))
        return Authors.from_names(names, False)

    return Authors.from_names([], False)


def _authors_from_string(text: str) -> Authors:
    chunks = [c.strip() for c in text.strip().split(",") if c.strip()]
    has_etal = any(c.lower().rstrip(".") == "et al" for c in chunks)
    name_chunks = [c for c in chunks if c.lower().rstrip(".") != "et al"]
    names = [(_surname_from_name(c), c) for c in name_chunks]
    return Authors.from_names(names, has_etal)


@dataclass
class Citation:
    """A single tool citation, parsed from one CSL-JSON item or BibTeX entry."""

    tool: str
    version: Optional[str] = None
    title: Optional[str] = None
    authors: Authors = field(default_factory=lambda: Authors([], "", False))
    year: Optional[int] = None
    container_title: Optional[str] = None
    doi: Optional[str] = None
    url: Optional[str] = None
    csl_type: Optional[str] = None

    @property
    def surnames(self) -> List[str]:
        return self.authors.surnames

    @property
    def author_display(self) -> str:
        return self.authors.display

    @property
    def is_multiple_authors(self) -> bool:
        return self.authors.is_multiple

    @classmethod
    def from_csl(cls, item: Dict[str, Any]) -> "Citation":
        custom = item.get("custom") or {}
        tool = custom.get("tool") or item.get("id")
        if not tool:
            raise ValueError(f"CSL item is missing both `custom.tool` and `id`: {item!r}")

        year: Optional[int] = None
        issued = item.get("issued") or {}
        date_parts = issued.get("date-parts") if isinstance(issued, dict) else None
        if date_parts and date_parts[0]:
            year = int(date_parts[0][0])

        return cls(
            tool=str(tool),
            version=custom.get("version"),
            title=item.get("title"),
            authors=_normalize_authors(item.get("author")),
            year=year,
            container_title=item.get("container-title"),
            doi=clean_doi(item.get("DOI")),
            url=item.get("URL"),
            csl_type=item.get("type"),
        )

    @property
    def has_publication(self) -> bool:
        """True when there is enough to form a Harvard short citation."""
        return bool(self.surnames) and self.year is not None

    @property
    def doi_url(self) -> Optional[str]:
        return f"https://doi.org/{self.doi}" if self.doi else None

    @property
    def link(self) -> Optional[str]:
        """Short-citation link target: the DOI if present, else the homepage."""
        return self.doi_url or self.url

    def short_citation_html(self) -> Optional[str]:
        """Harvard short form ("Surname (year)"), or None without a publication."""
        if not self.has_publication:
            return None
        if self.is_multiple_authors:
            name = f"{_esc(self.surnames[0])} <i>et al.</i>"
        else:
            name = _esc(self.surnames[0])
        return f"{name} ({self.year})"

    def inline_html(self) -> str:
        """One entry in the inline methods sentence."""
        short = self.short_citation_html()
        link = self.link
        if short:
            linked = f'<a href="{_esc_attr(link)}">{short}</a>' if link else short
            if self.version:
                return f"{_esc(self.tool)} ({_esc(self.version)}, {linked})"
            return f"{_esc(self.tool)} ({linked})"

        # No publication: identify by (linked) tool name rather than anonymously.
        tool = f'<a href="{_esc_attr(link)}">{_esc(self.tool)}</a>' if link else _esc(self.tool)
        return f"{tool} ({_esc(self.version)})" if self.version else tool

    def bibliography_html(self) -> str:
        """A full <li> body: one bibliography reference for this tool."""
        if self.has_publication:
            segments = []
            if self.author_display:
                # Avoid a doubled period when the author string ends in "et al."
                author = self.author_display if self.author_display.endswith(".") else f"{self.author_display}."
                segments.append(_esc(author))
            segments.append(f"({self.year}).")
            if self.title:
                segments.append(f"{_esc(self.title)}.")
            if self.container_title:
                segments.append(f"{_esc(self.container_title)}.")
            text = " ".join(segments)
            if self.doi_url:
                text += f' <a href="{_esc_attr(self.doi_url)}">{_esc(self.doi_url)}</a>'
            elif self.url:
                text += f' <a href="{_esc_attr(self.url)}">{_esc(self.url)}</a>'
            return text

        # No publication: lead with the tool name so it is not anonymous.
        if self.doi_url:
            return f'{_esc(self.tool)}. doi: <a href="{_esc_attr(self.doi_url)}">{_esc(self.doi_url)}</a>'
        if self.url:
            return f'{_esc(self.tool)}. Available at: <a href="{_esc_attr(self.url)}">{_esc(self.url)}</a>'
        return f"{_esc(self.tool)}."

    def to_dict(self) -> Dict[str, Any]:
        """Flat record for the data file."""
        return {
            "tool": self.tool,
            "version": self.version,
            "title": self.title,
            "author": self.author_display,
            "year": self.year,
            "source": self.container_title,
            "doi": self.doi,
            "url": self.url,
        }


def render_inline(citations: List[Citation]) -> str:
    """Render the inline methods sentence as HTML."""
    parts = [c.inline_html() for c in citations]
    return "Tools used in the workflow included: " + ", ".join(parts) + "."


def render_bibliography(citations: List[Citation]) -> str:
    """Render the bibliography as an ordered list."""
    items = "".join(f"<li>{c.bibliography_html()}</li>" for c in citations)
    return f'<ol class="citations-bibliography">{items}</ol>'
