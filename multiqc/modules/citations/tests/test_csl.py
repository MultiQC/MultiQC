"""Tests for CSL-JSON parsing and Harvard rendering in the citations module."""

import pytest

from multiqc.modules.citations.csl import (
    Citation,
    parse_csl,
    render_bibliography,
    render_inline,
)

SINGLE_AUTHOR = """
[
  {
    "id": "fastqc",
    "type": "software",
    "title": "FastQC: A Quality Control Tool for High Throughput Sequence Data",
    "author": [{"family": "Andrews", "given": "S"}],
    "issued": {"date-parts": [[2010]]},
    "container-title": "Bioinformatics",
    "DOI": "10.1093/bioinformatics/btv033",
    "URL": "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
    "custom": {"tool": "fastqc", "version": "0.12.1"}
  }
]
"""

MULTI_AUTHOR = """
[
  {
    "id": "samtools",
    "type": "article-journal",
    "title": "Twelve years of SAMtools and BCFtools",
    "author": [
      {"family": "Danecek", "given": "P"},
      {"family": "Bonfield", "given": "J K"},
      {"family": "Liddle", "given": "J"}
    ],
    "issued": {"date-parts": [[2021]]},
    "container-title": "GigaScience",
    "DOI": "10.1093/gigascience/giab008",
    "custom": {"tool": "samtools", "version": "1.21"}
  }
]
"""

DOI_ONLY = """
[{"id": "snpsift", "type": "software", "DOI": "10.5281/zenodo.1234567", "custom": {"tool": "snpsift", "version": "5.2"}}]
"""

HOMEPAGE_ONLY = """
[{"id": "customtool", "type": "software", "URL": "https://example.org/customtool", "custom": {"tool": "customtool", "version": "1.0"}}]
"""

PARTIAL = """
[{"id": "partialtool", "type": "article-journal", "title": "A partial record without a source or DOI",
  "author": [{"literal": "Smith J"}], "issued": {"date-parts": [[2019]]}, "custom": {"tool": "partialtool", "version": "3.1"}}]
"""


def test_single_author_inline():
    citations = parse_csl(SINGLE_AUTHOR)
    assert render_inline(citations) == (
        "Tools used in the workflow included: fastqc (0.12.1, "
        '<a href="https://doi.org/10.1093/bioinformatics/btv033">Andrews (2010)</a>).'
    )


def test_single_author_bibliography():
    citations = parse_csl(SINGLE_AUTHOR)
    assert render_bibliography(citations) == (
        '<ol class="citations-bibliography"><li>'
        "Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. Bioinformatics. "
        '<a href="https://doi.org/10.1093/bioinformatics/btv033">https://doi.org/10.1093/bioinformatics/btv033</a>'
        "</li></ol>"
    )


def test_multi_author_uses_et_al():
    c = parse_csl(MULTI_AUTHOR)[0]
    assert c.is_multiple_authors is True
    assert c.inline_html() == (
        'samtools (1.21, <a href="https://doi.org/10.1093/gigascience/giab008">Danecek <i>et al.</i> (2021)</a>)'
    )
    assert c.bibliography_html().startswith("Danecek P, Bonfield J K, Liddle J. (2021).")


def test_doi_only_identified_by_tool_name():
    c = parse_csl(DOI_ONLY)[0]
    assert c.has_publication is False
    assert c.inline_html() == '<a href="https://doi.org/10.5281/zenodo.1234567">snpsift</a> (5.2)'
    assert c.bibliography_html() == (
        'snpsift. doi: <a href="https://doi.org/10.5281/zenodo.1234567">https://doi.org/10.5281/zenodo.1234567</a>'
    )


def test_homepage_only_fallback():
    c = parse_csl(HOMEPAGE_ONLY)[0]
    assert c.inline_html() == '<a href="https://example.org/customtool">customtool</a> (1.0)'
    assert c.bibliography_html() == (
        'customtool. Available at: <a href="https://example.org/customtool">https://example.org/customtool</a>'
    )


def test_partial_publication_renders_present_fields():
    c = parse_csl(PARTIAL)[0]
    assert c.inline_html() == "partialtool (3.1, Smith (2019))"
    assert c.bibliography_html() == "Smith J. (2019). A partial record without a source or DOI."


def test_plain_string_author_etal_no_double_period():
    item = {
        "id": "multiqc",
        "title": "MultiQC: summarize analysis results for multiple tools and samples in a single report",
        "author": "Ewels P, Magnusson M, Lundin S, et al.",
        "issued": {"date-parts": [[2016]]},
        "container-title": "Bioinformatics",
        "DOI": "10.1093/bioinformatics/btw354",
        "custom": {"tool": "multiqc", "version": "1.21"},
    }
    c = Citation.from_csl(item)
    assert c.surnames[0] == "Ewels"
    assert c.is_multiple_authors is True
    bib = c.bibliography_html()
    assert "et al.. " not in bib  # regression: no doubled period
    assert bib.startswith("Ewels P, Magnusson M, Lundin S, et al. (2016).")


def test_doi_is_normalised():
    assert Citation.from_csl({"id": "x", "DOI": "https://doi.org/10.1/abc", "custom": {"tool": "x"}}).doi == "10.1/abc"
    assert Citation.from_csl({"id": "y", "DOI": "doi: 10.2/def", "custom": {"tool": "y"}}).doi == "10.2/def"


def test_link_prefers_doi_then_homepage():
    both = Citation.from_csl({"id": "t", "DOI": "10.1/a", "URL": "https://h", "custom": {"tool": "t"}})
    assert both.link == "https://doi.org/10.1/a"
    assert Citation.from_csl({"id": "t", "URL": "https://h", "custom": {"tool": "t"}}).link == "https://h"
    assert Citation.from_csl({"id": "t", "custom": {"tool": "t"}}).link is None


def test_accepts_object_keyed_by_id_and_single_item():
    keyed = '{"fastqc": {"id": "fastqc", "custom": {"tool": "fastqc", "version": "1.0"}}}'
    assert [c.tool for c in parse_csl(keyed)] == ["fastqc"]
    single = '{"id": "fastqc", "custom": {"tool": "fastqc", "version": "1.0"}}'
    assert [c.tool for c in parse_csl(single)] == ["fastqc"]


def test_malformed_json_raises_with_filename():
    with pytest.raises(ValueError, match="bad.csl.json"):
        parse_csl("{not json", "bad.csl.json")


def test_html_is_escaped():
    item = {
        "id": "evil",
        "title": 'A <script> & "quotes"',
        "author": [{"family": "O'Brien", "given": "X"}],
        "issued": {"date-parts": [[2020]]},
        "custom": {"tool": "evil"},
    }
    bib = Citation.from_csl(item).bibliography_html()
    assert "<script>" not in bib
    assert "&lt;script&gt;" in bib
