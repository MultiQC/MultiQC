import pytest
from pathlib import Path
from multiqc import report
from multiqc.modules.proteinfold import MultiqcModule

# Will fix when I do a PR against the MultiQC test data repo
# TEST_DATA_DIR = Path(__file__).parent / "test-data"
TEST_DATA_DIR = Path("/srv/scratch/z3374843/test-data/data/modules/proteinfold")


@pytest.fixture
def mod():
    """Fixture to create a MultiqcModule instance with test data."""
    # Reset report state
    report.analysis_files = []
    report.files = {}

    # Set up analysis files from local test data
    report.analysis_files = list(TEST_DATA_DIR.glob("*"))
    report.search_files(["proteinfold"])

    return MultiqcModule()


def test_proteinfold_msa(mod):
    """Test MSA depth is parsed correctly from coverage files."""
    assert len(mod.proteinfold_data) > 0
    assert "T1024" in mod.proteinfold_data
    assert "msa_depth" in mod.proteinfold_data["T1024"]
    assert mod.proteinfold_data["T1024"]["msa_depth"] == 9226


def test_proteinfold_plddt(mod):
    """Test pLDDT scores are parsed correctly from TSV files."""
    assert "T1026" in mod.proteinfold_data
    assert "mean_plddt" in mod.proteinfold_data["T1026"]
    assert "plddt" in mod.proteinfold_data["T1026"]

    plddt_data = mod.proteinfold_data["T1026"]["plddt"]
    assert isinstance(plddt_data["rank_0"], dict)

    # Test top-ranked result mean
    assert mod.proteinfold_data["T1026"]["mean_plddt"] == pytest.approx(75.676, rel=0.001)

    # Parent mean should match rank_0 mean
    assert mod.proteinfold_data["T1026_rank_0"]["mean_plddt"] == pytest.approx(
        mod.proteinfold_data["T1026"]["mean_plddt"], rel=0.00001
    )

    # Test specific position values in rank_0
    assert plddt_data["rank_0"][0] == pytest.approx(34.31, rel=0.01)
    assert plddt_data["rank_0"][50] == pytest.approx(83.67, rel=0.01)

    # Test another rank's position values
    assert plddt_data["rank_2"][0] == pytest.approx(35.51, rel=0.01)
    assert plddt_data["rank_2"][50] == pytest.approx(81.88, rel=0.01)


def test_proteinfold_iptm(mod):
    """Test ipTM (interface pTM) scores are parsed correctly."""
    assert "Jun-Jun" in mod.proteinfold_data
    assert "iptm" in mod.proteinfold_data["Jun-Jun"]

    # Test top-ranked result
    assert mod.proteinfold_data["Jun-Jun"]["iptm"] == pytest.approx(0.618, rel=0.001)

    # Parent should match rank_0
    assert mod.proteinfold_data["Jun-Jun_rank_0"]["iptm"] == pytest.approx(
        mod.proteinfold_data["Jun-Jun"]["iptm"], rel=0.00001
    )

    # Test another rank
    assert mod.proteinfold_data["Jun-Jun_rank_3"]["iptm"] == pytest.approx(0.606, rel=0.001)


def test_proteinfold_ptm(mod):
    """Test pTM (predicted TM-score) values are parsed correctly."""
    assert "Fos-Fos" in mod.proteinfold_data
    assert "ptm" in mod.proteinfold_data["Fos-Fos"]

    # Test top-ranked result
    assert mod.proteinfold_data["Fos-Fos"]["ptm"] == pytest.approx(0.7, rel=0.001)

    # Parent should match rank_0
    assert mod.proteinfold_data["Fos-Fos_rank_0"]["ptm"] == pytest.approx(
        mod.proteinfold_data["Fos-Fos"]["ptm"], rel=0.00001
    )

    # Test another rank
    assert mod.proteinfold_data["Fos-Fos_rank_4"]["ptm"] == pytest.approx(0.696, rel=0.001)
