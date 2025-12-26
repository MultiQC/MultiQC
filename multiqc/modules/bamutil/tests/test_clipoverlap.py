from multiqc.modules.bamutil import clipoverlap

bamutil_stderr1 = """
Overlap Statistics:
Number of overlapping pairs: 29
Average # Reference Bases Overlapped: 36.069
Variance of Reference Bases overlapped: 407.638
Number of times orientation causes additional clipping: 0
Number of times the forward strand was clipped: 16
Number of times the reverse strand was clipped: 13
Completed ClipOverlap Successfully.
"""

bamutil_stderr2 = """WARNING - 
Problems encountered parsing command line:

Command line parameter logs/clipOverlap/CZA569_clip.log (#7) ignored

Overlap Statistics:
Number of overlapping pairs: 1543210
Average # Reference Bases Overlapped: 82.341
Variance of Reference Bases overlapped: 983.25
Number of times orientation causes additional clipping: 39875
Number of times the forward strand was clipped: 987634
Number of times the reverse strand was clipped: 555576
Completed ClipOverlap Successfully.
"""

bamutil_stderr3 = """Number of overlapping pairs: 1782320
Average # Reference Bases Overlapped: 79.5654
Variance of Reference Bases overlapped: 1061.87
Number of times orientation causes additional clipping: 47695
Number of times the forward strand was clipped: 1065746
Number of times the reverse strand was clipped: 716574
WARNING: did not find expected overlapping mates for 5145 records.
Completed ClipOverlap Successfully.
"""


def test_parse_logs(monkeypatch):
    class DummyFile:
        def __init__(self, lines):
            self.lines = lines

        def __iter__(self):
            return iter(self.lines)

    def fake_find_log_files(self, pattern, filehandles=False):
        return [
            {
                "s_name": "sample1",
                "f": DummyFile(bamutil_stderr1.splitlines()),
                "root": "/fake/path",
                "fn": "sample1.log",
            },
            {
                "s_name": "sample2",
                "f": DummyFile(bamutil_stderr2.splitlines()),
                "root": "/fake/path",
                "fn": "sample2.log",
            },
            {
                "s_name": "sample3",
                "f": DummyFile(bamutil_stderr3.splitlines()),
                "root": "/fake/path",
                "fn": "sample3.log",
            },
        ]

    monkeypatch.setattr(clipoverlap.MultiqcModule, "find_log_files", fake_find_log_files)
    monkeypatch.setattr(clipoverlap.MultiqcModule, "clean_s_name", lambda self, x: x)

    # Create the module
    mod = clipoverlap.MultiqcModule()

    # Check first sample
    assert "sample1" in mod.bamutil_data
    assert mod.bamutil_data["sample1"]["bamutil_overlapping_pairs"] == 29
    assert mod.bamutil_data["sample1"]["bamutil_avg_ref_bases_overlapped"] == 36.069
    assert mod.bamutil_data["sample1"]["bamutil_variance_ref_bases"] == 407.638
    assert mod.bamutil_data["sample1"]["bamutil_orientation_additional_clipping"] == 0
    assert mod.bamutil_data["sample1"]["bamutil_forward_strand_clipped"] == 16
    assert mod.bamutil_data["sample1"]["bamutil_reverse_strand_clipped"] == 13

    # Check second sample with warning messages
    assert "sample2" in mod.bamutil_data
    assert mod.bamutil_data["sample2"]["bamutil_overlapping_pairs"] == 1543210
    assert mod.bamutil_data["sample2"]["bamutil_avg_ref_bases_overlapped"] == 82.341

    # Check third sample with missing mate warning
    assert "sample3" in mod.bamutil_data
    assert mod.bamutil_data["sample3"]["bamutil_overlapping_pairs"] == 1782320
    assert mod.bamutil_data["sample3"]["bamutil_missing_overlapping_mates"] == 5145
    assert mod.bamutil_data["sample3"]["bamutil_has_warnings"] is True
