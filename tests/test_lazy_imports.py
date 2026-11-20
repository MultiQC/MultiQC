"""
Check that `__lazy_modules__` (PEP 810) keeps heavy dependencies out of `import multiqc`.

Lazy imports only exist on Python 3.15+; on older versions `__lazy_modules__` is ignored.
"""

import subprocess
import sys

import pytest

# Only needed to build a report, so they should not be loaded by a bare `import multiqc`
DEFERRED_MODULES = [
    "jinja2",
    "jsonschema",
    "markdown",
    "multiqc.base_module",
    "multiqc.interactive",
    "multiqc.report",
    "natsort",
    "numpy",
    "polars",
    "requests",
    "spectra",
]


@pytest.mark.skipif(sys.version_info < (3, 15), reason="lazy imports require Python 3.15+")
def test_import_multiqc_defers_heavy_modules():
    # Run in a subprocess: the test session has already imported everything
    code = "import sys, multiqc; print(' '.join(m for m in sys.argv[1:] if m in sys.modules))"
    result = subprocess.run([sys.executable, "-c", code, *DEFERRED_MODULES], capture_output=True, text=True, check=True)
    assert result.stdout.split() == []
