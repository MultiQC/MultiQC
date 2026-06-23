"""Tests for BibTeX parsing in the citations module.

Skipped when the optional `bibtexparser` dependency is not installed.
"""

import pytest

from multiqc.modules.citations.bibtex import (
    _split_bibtex_name,
    bibtexparser_available,
    parse_bibtex,
)
from multiqc.modules.citations.csl import render_bibliography, render_inline

pytestmark = pytest.mark.skipif(not bibtexparser_available(), reason="bibtexparser not installed")

SAMPLE = """
@software{fastqc,
  title    = {{FastQC}: A Quality Control Tool for High Throughput Sequence Data},
  author   = {Andrews, S},
  year     = {2010},
  journal  = {Bioinformatics},
  doi      = {10.1093/bioinformatics/btv033},
  url      = {https://www.bioinformatics.babraham.ac.uk/projects/fastqc/},
  version  = {0.12.1}
}

@article{samtools,
  title    = {Twelve years of {SAMtools} and {BCFtools}},
  author   = {Danecek, Petr and Bonfield, James K and others},
  year     = {2021},
  journal  = {GigaScience},
  doi      = {10.1093/gigascience/giab008},
  version  = {1.21}
}

@software{snpsift,
  doi      = {10.5281/zenodo.1234567},
  version  = {5.2}
}
"""


@pytest.fixture
def citations():
    return {c.tool: c for c in parse_bibtex(SAMPLE)}


def test_single_author_entry(citations):
    fastqc = citations["fastqc"]
    assert fastqc.version == "0.12.1"
    assert fastqc.title == "FastQC: A Quality Control Tool for High Throughput Sequence Data"
    assert fastqc.surnames == ["Andrews"]
    assert fastqc.container_title == "Bioinformatics"
    assert fastqc.inline_html() == (
        'fastqc (0.12.1, <a href="https://doi.org/10.1093/bioinformatics/btv033">Andrews (2010)</a>)'
    )


def test_multi_author_and_others(citations):
    samtools = citations["samtools"]
    assert samtools.surnames == ["Danecek", "Bonfield"]
    assert samtools.is_multiple_authors is True
    assert samtools.author_display == "Danecek Petr, Bonfield James K, et al."
    bib = samtools.bibliography_html()
    assert "et al.. " not in bib  # regression: no doubled period
    assert bib.startswith("Danecek Petr, Bonfield James K, et al. (2021).")


def test_doi_only_entry(citations):
    snpsift = citations["snpsift"]
    assert snpsift.has_publication is False
    assert snpsift.inline_html() == '<a href="https://doi.org/10.5281/zenodo.1234567">snpsift</a> (5.2)'


def test_renderers_shared_with_csl(citations):
    ordered = [citations[t] for t in sorted(citations, key=str.lower)]
    assert render_inline(ordered).startswith("Tools used in the workflow included: fastqc (0.12.1,")
    assert render_bibliography(ordered).count("<li>") == 3


@pytest.mark.parametrize(
    "name,expected",
    [
        ("Andrews, S", ("Andrews", "S")),
        ("Heng Li", ("Li", "Heng")),
        ("Madonna", ("Madonna", "")),
        ("van der Berg, A", ("van der Berg", "A")),
    ],
)
def test_split_bibtex_name(name, expected):
    assert _split_bibtex_name(name) == expected


def test_tool_field_overrides_citation_key():
    citations = {c.tool: c for c in parse_bibtex("@misc{somekey, tool = {realtool}, version = {2.0}, doi = {10.1/x}}")}
    assert "realtool" in citations
    assert citations["realtool"].version == "2.0"
