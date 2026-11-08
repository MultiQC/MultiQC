import csv
import io
import json
import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.mark.parametrize("format,delimiter", [("csv", ","), ("tsv", "\t")])
def test_violin_export_units(format, delimiter):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to test browser export code")
    source = Path(__file__).parents[1] / "multiqc/templates/default/src/js/plots/violin.js"
    script = """
const fs = require('node:fs');
const vm = require('node:vm');
const context = vm.createContext({
  Plot: class {}, window: {}, $: () => ({length: 0}),
  applyToolboxSettings: samples => samples.map(name => ({hidden: name === 'hidden'}))
});
vm.runInContext(fs.readFileSync(process.argv[1], 'utf8'), context);
const plot = Object.create(context.window.ViolinPlot.prototype);
plot.tableAnchor = 'test';
plot.activeDatasetIdx = 0;
plot.datasets = [{
  metrics: ['reads', 'percent', 'unitless', 'quoted'],
  header_by_metric: {
    reads: {title: 'Seqs', suffix: ' M'},
    percent: {title: 'GC (%)', suffix: '%'},
    unitless: {title: 'Score', suffix: ''},
    quoted: {title: 'Value "quoted", unit', suffix: 'a,b"c'}
  },
  all_samples: ['sample', 'missing', 'hidden'],
  violin_value_by_sample_by_metric: {
    reads: {sample: 12.345678, hidden: 99},
    percent: {sample: 48, hidden: 50},
    unitless: {sample: 0, hidden: 1},
    quoted: {sample: 3, hidden: 4}
  }
}];
console.log(JSON.stringify(plot.exportData(process.argv[2])));
"""
    result = subprocess.run(
        [node, "-e", script, str(source), format],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = list(csv.reader(io.StringIO(json.loads(result.stdout)), delimiter=delimiter))
    assert rows == [
        ["Sample", "Seqs (M)", "GC (%)", "Score", 'Value "quoted", unit (a,b"c)'],
        ["sample", "12.345678", "48", "0", "3"],
        ["missing", ".", ".", ".", "."],
    ]
