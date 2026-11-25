import pytest

from multiqc import config, report
from multiqc.modules.suppa import MultiqcModule
from multiqc.modules.suppa.diffsplice import parse_dpsi
from multiqc.modules.suppa.events import ISOFORM, event_type
from multiqc.modules.suppa.psi import parse_psi

PSI_EVENTS = """S1\tS2\tS3
ENSG01;SE:chr1:100-200:300-400:+\t0.5\t1.0000000000000002\tnan
ENSG01;A5:chr1:100-200:150-200:+\t0.0\tNA\t0.25
ENSG02;MX:chr2:1-2:3-4:5-6:7-8:-\t1.0\t0.75\t-1.0
"""

PSI_ISOFORMS = """S1\tS2
ENSG01;ENST01\t0.9\tnan
ENSG01;ENST02\t0.1\tnan
"""

DPSI_LOCAL = """c1-c2_dPSI\tc1-c2_p-val
ENSG01;SE:chr1:100-200:300-400:+\t0.5\t0.01
ENSG01;A5:chr1:100-200:150-200:+\t-0.3\t0.001
ENSG01;A3:chr1:100-200:150-200:+\t0.04\t0.001
ENSG02;RI:chr2:1-2:3-4:-\t0.6\t0.2
ENSG03;AF:chr3:1-2:3-4:5-6:+\tnan\t1.0
"""

DPSI_MULTI = """c1-c2_dPSI\tc1-c2_p-val\tc2-c3_dPSI\tc2-c3_p-val
ENSG01;ENST01\t0.5\t0.01\tnan\tnan
ENSG01;ENST02\t-0.5\t0.01\t0.2\t0.04
"""


def test_event_type():
    assert event_type("ENSG01;SE:chr1:100-200:300-400:+") == "Skipped exon"
    assert event_type("ENSG01;AL:chr1:1:2-3:4:5-6:+") == "Alternative last exon"
    assert event_type("ENSG01;ENST01") == ISOFORM


def test_parse_psi_events():
    level, data = parse_psi(PSI_EVENTS.splitlines(keepends=True))
    assert level == "event"
    assert data["S1"]["events"] == 3
    assert data["S1"]["quantified"] == 3
    assert data["S2"]["quantified"] == 2
    # nan and the legacy -1 sentinel are not quantified
    assert data["S3"]["quantified"] == 1
    assert data["S3"]["pct_quantified"] == pytest.approx(100 / 3)
    assert data["S1"]["mean_psi"] == pytest.approx(0.5)
    # 1.0000000000000002 lands in the top bin
    assert data["S2"]["histogram"][0.975] == pytest.approx(50.0)
    assert sum(data["S1"]["histogram"].values()) == pytest.approx(100.0)


def test_parse_psi_isoforms_all_missing():
    level, data = parse_psi(PSI_ISOFORMS.splitlines(keepends=True))
    assert level == "isoform"
    assert data["S2"]["quantified"] == 0
    assert "mean_psi" not in data["S2"]
    assert "histogram" not in data["S2"]


@pytest.mark.parametrize(
    "content, message",
    [
        ("", "empty"),
        ("S1\tS2\n", "no events"),
        ("S1\tS1\nENSG01;ENST01\t0.5\t0.5\n", "duplicate"),
        ("S1\tS2\nENSG01;ENST01\t0.5\n", "1 values for 2 samples"),
        ("S1\nENSG01;ENST01\t1.5\n", "outside"),
        ("S1\nENSG01;ENST01\tabc\n", "could not convert"),
        ("S1\nENSG01;ENST01\t0.5\nENSG01;SE:chr1:1-2:3-4:+\t0.5\n", "mixes"),
    ],
)
def test_parse_psi_malformed(content, message):
    with pytest.raises(ValueError, match=message):
        parse_psi(content.splitlines(keepends=True))


def test_parse_dpsi_local():
    data = parse_dpsi(DPSI_LOCAL.splitlines(keepends=True), 0.05, 0.05)
    d = data["c1-c2"]
    assert d["events"] == 5
    assert d["tested"] == 4
    # A3 fails |dPSI| > 0.05, RI fails p < 0.05
    assert d["significant"] == 2
    assert d["significant_up"] == 1
    assert d["significant_down"] == 1
    assert d["significant_by_type"] == {"Skipped exon": 1, "Alternative 5' splice site": 1}
    assert d["tested_by_type"]["Retained intron"] == 1
    assert "Alternative first exon" not in d["tested_by_type"]
    assert d["pct_tested"] == pytest.approx(80.0)
    assert d["pct_significant"] == pytest.approx(50.0)


def test_parse_dpsi_cutoffs():
    data = parse_dpsi(DPSI_LOCAL.splitlines(keepends=True), 0.5, 0.0)
    assert data["c1-c2"]["significant"] == 4


def test_parse_dpsi_multiple_comparisons():
    data = parse_dpsi(DPSI_MULTI.splitlines(keepends=True), 0.05, 0.05)
    assert list(data) == ["c1-c2", "c2-c3"]
    assert data["c1-c2"]["significant_by_type"] == {ISOFORM: 2}
    assert data["c2-c3"]["tested"] == 1
    assert data["c2-c3"]["significant"] == 1


@pytest.mark.parametrize(
    "content, message",
    [
        ("", "empty"),
        ("c1-c2_dPSI\n", "column pairs"),
        ("c1-c2_dPSI\tc2-c3_p-val\n", "unexpected column pair"),
        ("c1-c2_dPSI\tc1-c2_p-val\nENSG01;ENST01\t0.5\n", "expected 2"),
    ],
)
def test_parse_dpsi_malformed(content, message):
    with pytest.raises(ValueError, match=message):
        parse_dpsi(content.splitlines(keepends=True), 0.05, 0.05)


@pytest.fixture
def run_module(tmp_path):
    def _run(files: dict[str, str], module_config: dict | None = None, prepend_dirs: bool = False):
        for name, content in files.items():
            (tmp_path / name).write_text(content)
        report.reset()
        config.reset()
        config.preserve_module_raw_data = True
        config.prepend_dirs = prepend_dirs
        if module_config is not None:
            config.update({"suppa": module_config})
        report.analysis_files = [tmp_path]
        report.search_files(["suppa"])
        return MultiqcModule()

    yield _run
    config.reset()


def test_module_sections(run_module):
    module = run_module({"tpm.psi": PSI_EVENTS, "c1-c2.dpsi": DPSI_LOCAL})
    assert [s.anchor for s in module.sections] == [
        "suppa_psi_quantification",
        "suppa_psi_distribution",
        "suppa_differential_splicing",
        "suppa_significant_events",
    ]
    assert module.saved_raw_data["multiqc_suppa_diffsplice"]["c1-c2"]["significant"] == 2


def test_module_cutoff_config(run_module):
    module = run_module({"c1-c2.dpsi": DPSI_LOCAL}, {"pvalue_cutoff": 0.5, "dpsi_cutoff": 0.0})
    assert module.saved_raw_data["multiqc_suppa_diffsplice"]["c1-c2"]["significant"] == 4


def test_module_prepend_dirs(run_module, tmp_path):
    module = run_module({"c1-c2.dpsi": DPSI_LOCAL, "multi.dpsi": DPSI_MULTI}, prepend_dirs=True)
    names = list(module.saved_raw_data["multiqc_suppa_diffsplice"])
    assert len(names) == 3
    for name in names:
        assert name.count(tmp_path.name) == 1


def test_module_skips_malformed_file(run_module):
    module = run_module({"c1-c2.dpsi": DPSI_LOCAL, "bad.dpsi": "c1-c2_dPSI\tc1-c2_p-val\nENSG01;ENST01\tabc\t0.1\n"})
    assert list(module.saved_raw_data["multiqc_suppa_diffsplice"]) == ["c1-c2"]
