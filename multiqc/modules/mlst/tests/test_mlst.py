from pathlib import Path

from multiqc import report, reset
from multiqc.modules.mlst.mlst import MultiqcModule


def test_mlst_module(tmp_path: Path) -> None:
    reset()

    sample_log = tmp_path / "Acinetobacter_GCF_008632635.1_mlst.txt"
    sample_log.write_text(
        "FILE\tSCHEME\tST\tSTATUS\tSCORE\tALLELES\n"
        "input_Acinetobacter_GCF_008632635.1.fasta\tabaumannii_2\t46\tPERFECT\t100\t"
        "Pas_cpn60(5);Pas_fusA(12);Pas_gltA(11);Pas_pyrG(2);Pas_recA(14);Pas_rplB(9);Pas_rpoB(14)\n"
    )

    report.analysis_files = [tmp_path]
    report.search_files(["mlst"])
    mod = MultiqcModule()

    sample_key = "input_Acinetobacter_GCF_008632635.1"
    assert sample_key in mod.mlst_data

    parsed = mod.mlst_data[sample_key]
    assert parsed["scheme"] == "abaumannii_2"
    assert parsed["st"] == "46"
    assert parsed["status"] == "PERFECT"
    assert parsed["score"] == 100.0