from multiqc import report
from multiqc.modules.varscan2 import MultiqcModule

LOG_TEMPLATE = """\
Only SNPs will be reported
Min coverage:\t3
Min reads2:\t2
Min var freq:\t0.8
Min avg qual:\t20
P-value thresh:\t0.05
{input_line}12928920 bases in pileup file
204806 variant positions (200602 SNP, 4350 indel)
3348 were failed by the strand-filter
197412 variant positions reported (197412 SNP, 138 indel)
"""


def test_sample_name_without_input_file_name(tmp_path):
    """Input file names are used as sample names, including absolute paths.

    Logs that don't name a regular input file fall back to the log file name.

    This happens when the pileup is piped in, e.g. with process substitution
    (https://github.com/MultiQC/MultiQC/issues/1673) or from stdin.
    """
    (tmp_path / "sample1.log").write_text(LOG_TEMPLATE.format(input_line="Reading input from sample1.mpileup\n"))
    (tmp_path / "sample2.log").write_text(LOG_TEMPLATE.format(input_line="Reading input from /data/sample2.mpileup\n"))
    (tmp_path / "Lu15.log").write_text(LOG_TEMPLATE.format(input_line="Reading input from /dev/fd/63\n"))
    (tmp_path / "Lu16.log").write_text(LOG_TEMPLATE.format(input_line="Reading input from /dev/fd/63\n"))
    (tmp_path / "Lu17.log").write_text(LOG_TEMPLATE.format(input_line=""))

    report.reset()
    report.analysis_files = [tmp_path]
    report.search_files(["varscan2"])
    m = MultiqcModule()

    assert set(m.varscan2_data) == {"sample1", "sample2", "Lu15", "Lu16", "Lu17"}
    assert m.varscan2_data["Lu15"]["variant_reported_snps"] == 197412
    assert "sample_name" not in m.varscan2_data["Lu15"]
