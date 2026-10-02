import multiqc
from multiqc.modules.phantompeakqualtools.phantompeakqualtools import MultiqcModule


def test_comment_header_line_is_skipped(data_dir):
    multiqc.parse_logs(data_dir / "modules/phantompeakqualtools/issue_1295")
    mod = MultiqcModule()
    assert mod.phantompeakqualtools_data == {
        "Example_header": {"Estimated_Fragment_Length_bp": 85, "NSC": 1.076593, "RSC": 1.273581}
    }
