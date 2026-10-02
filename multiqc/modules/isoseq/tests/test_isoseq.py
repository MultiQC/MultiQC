import json

import pytest

from multiqc import report
from multiqc.modules.isoseq.isoseq import MultiqcModule, flatten_refine_json
from multiqc.utils import testing

METRICS = {"num_reads_fl", "num_reads_flnc", "num_reads_flnc_polya"}


@pytest.mark.parametrize(
    "path",
    [
        "v3/refine/alz.1perc.subreads.10000.chunk1.filter_summary.json",
        "issue_3673/sample.filter_summary.report.json",
    ],
)
def test_flatten_refine_json(path):
    with open(testing.data_dir() / "modules" / "isoseq" / path) as f:
        data = flatten_refine_json(json.load(f))
    assert METRICS <= set(data)
    assert all(isinstance(data[k], int) for k in METRICS)


def test_refine_pbreports_json_module():
    report.reset()
    report.analysis_files = [str(testing.data_dir() / "modules" / "isoseq" / "issue_3673")]
    report.search_files(["isoseq"])
    MultiqcModule()
