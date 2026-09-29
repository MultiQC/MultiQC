import pytest

from multiqc import config, report, reset
from multiqc.modules.cramino.cramino import MultiqcModule

V2_ALIGNED_WITH_KARYOTYPE = """File name\tsample.bam
Number of alignments\t1000
% from total reads\t96.00
Number of reads\t950
Yield [Gb]\t1.50
Mean coverage\t5.00
Yield [Gb] (>25kb)\t0.20
N50\t600
N75\t500
Median length\t550.00
Mean length\t700.00
N50 aligned\t600
N75 aligned\t500
Median length aligned\t550.00
Mean length aligned\t700.00

Median identity\t99.40
Mean identity\t99.00
Modal identity\t100.0



# Normalized read count per chromosome

{contigs}

Path\t/work/sample.bam
Creation time\tNA
"""


@pytest.fixture(autouse=True)
def _reset_report():
    """Isolate from file-search state left behind by other tests running in the same process."""
    reset()


def run_module(tmp_path, contigs, monkeypatch, cramino_config=None):
    contig_lines = "\n".join(f"{name}\t{value}" for name, value in contigs.items())
    (tmp_path / "sample.cramino.txt").write_text(V2_ALIGNED_WITH_KARYOTYPE.format(contigs=contig_lines))
    monkeypatch.setattr(config, "preserve_module_raw_data", True)
    if cramino_config is not None:
        monkeypatch.setattr(config, "cramino_config", cramino_config, raising=False)
    report.analysis_files = [tmp_path]
    report.search_files(["cramino"])
    return MultiqcModule()


def test_v2_report_end_to_end(tmp_path, monkeypatch):
    module = run_module(tmp_path, {"chr1": 0.9, "chrM": 71.0}, monkeypatch)

    assert {section.anchor for section in module.sections} == {"cramino-stats", "cramino-karyotype"}
    assert module.saved_raw_data is not None
    assert module.saved_raw_data["multiqc_cramino"]["sample"]["pct_from_total_alignments"] == 96.0
    # The default contig exclusion only applies to the plot, not to the exported data
    assert module.saved_raw_data["multiqc_cramino_karyotype"]["sample"] == {"chr1": 0.9, "chrM": 71.0}


def test_patterns_removing_all_contigs_keep_all_with_warning(tmp_path, caplog, monkeypatch):
    module = run_module(tmp_path, {"1": 0.9, "2": 1.1}, monkeypatch, {"include_contigs": ["chr*"]})

    karyotype_section = next(section for section in module.sections if section.anchor == "cramino-karyotype")
    assert karyotype_section.plot_anchor is not None
    assert "Keeping all contigs" in caplog.text


def test_unexpected_format_raises_with_file_name(tmp_path, monkeypatch):
    (tmp_path / "sample.cramino.txt").write_text(
        "File name\tsample.bam\nNumber of alignments\t10\nNumber of reads\t10\n"
    )
    report.analysis_files = [tmp_path]
    report.search_files(["cramino"])

    with pytest.raises(ValueError, match="sample.cramino.txt"):
        MultiqcModule()
