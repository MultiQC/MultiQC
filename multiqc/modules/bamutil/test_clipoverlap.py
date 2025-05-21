import sys
sys.path.insert(0, '/workspaces/MultiQC/')
from multiqc.modules.bamutil import clipoverlap

bamutil_stderr = '''
Overlap Statistics:
Number of overlapping pairs: 29
Average # Reference Bases Overlapped: 36.069
Variance of Reference Bases overlapped: 407.638
Number of times orientation causes additional clipping: 0
Number of times the forward strand was clipped: 16
Number of times the reverse strand was clipped: 13
Completed ClipOverlap Successfully.
'''

def test_parse_logs(monkeypatch):
    class DummyFile:
        def __init__(self, lines):
            self.lines = lines
        def __iter__(self):
            return iter(self.lines)
    
    def fake_find_log_files(self, pattern, filehandles=False):
        return [{
            "s_name": "sample1",
            "f": DummyFile(bamutil_stderr.splitlines())
        }]
    
    monkeypatch.setattr(clipoverlap.MultiqcModule, "find_log_files", fake_find_log_files)
    monkeypatch.setattr(clipoverlap.MultiqcModule, "clean_s_name", lambda self, x: x)
    mod = clipoverlap.MultiqcModule()
    assert "sample1" in mod.bamutil_data
    assert mod.bamutil_data["sample1"]["bamutil_overlapping_pairs"] == 29
    assert mod.bamutil_data["sample1"]["bamutil_avg_ref_bases_overlapped"] == 36.069
