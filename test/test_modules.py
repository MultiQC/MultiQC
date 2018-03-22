# -*- coding: utf-8 -*-
# Copyright (C) 2017 by Per Unneberg

"""
test_modules
----------------------------------------------------------

Test MultiQC modules parsing

"""
from __future__ import print_function
import re
import os
import pytest
import logging
import importlib
import multiqc.modules
import multiqc.utils.report
import multiqc.utils.util_functions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("multiqc.test.test_modules")

MULTIQC_MODULES_DIR = os.path.dirname(multiqc.modules.__file__)
MULTIQC_MODULES = set([x for x in os.listdir(MULTIQC_MODULES_DIR) if not x.endswith(".py") and not x == "__pycache__"])

has_ngsfixtures = False
fixtures = []
test_fixtures = []

try:
    import pytest_ngsfixtures
    has_ngsfixtures = True
    from pytest_ngsfixtures import factories
    from pytest_ngsfixtures.config import application_fixtures
    fixtures = application_fixtures(end="pe")
    modules = list(set([f[0] for f in fixtures]))
    pytest_modules = [pytest.config.option.module] if pytest.config.option.module else modules
except ImportError as e:
    print("\n\n   pytest-ngsfixtures not installed; install with 'conda install -c percyfal pytest-ngsfixtures'\n\n")

submodule_map = {'feature_counts': 'subread_featureCounts'}

def multiqc_module_command_exists(name, mod, command):
    """Assert whether multiqc module command actually exists.

    There is a mismatch in module content between MultiQC and
    pytest_ngsfixtures. This function checks whether applications
    available in pytest_ngsfixtures also exist in MultiQC.

    Params:
      name (str): pytest_ngsfixtures module name, remapped to MultiQC
        equivalent in case of naming inconsistency
      mod (str): MultiQC module name
      command (str): pytest_ngsfixtures command name, to be matched
        against MultiQC submodules

    Returns:
      bool: whether or not MultiQC module command exists
    """
    loader = importlib.find_loader(mod)
    m = loader.load_module()
    if m is None:
        # For now assume we can do something with data
        return True
    moduledir = os.path.dirname(m.__file__)
    submodules = [os.path.basename(x).replace(".py", "") for x in os.listdir(moduledir)]
    for submod in submodules:
        if submod.startswith("__"):
            continue
        if re.search(submod, command.replace("{}_".format(name), "")):
            return True
        elif submodule_map.get(submod, None) == command:
            return True
    return False

supported_multiqc_modules = []
missing_multiqc_parsers = []

remap = {'bowtie': 'bowtie1', 'fastq-screen': 'fastq_screen', 'subread': 'featureCounts'}

# Loop pytest_ngsfixtures and ascertain which ones are supported by
# MultiQC. Also, keep track of MultiQC modules that lack support by
# pytest_ngsfixtures
for f in fixtures:
    module, command, version, end, fmt = f
    mod = "multiqc.modules.{}".format(remap.get(module, module))
    try:
        importlib.import_module(mod)
        supported_multiqc_modules.append(module)
        if module in pytest_modules and multiqc_module_command_exists(remap.get(module, module), mod, remap.get(command, command)):
            test_fixtures.append(f)
    except ImportError as e:
        # Test data available in pytest-ngsfixtures but no parser in multiqc
        missing_multiqc_parsers.append(module)


missing_modules = MULTIQC_MODULES.difference(set(supported_multiqc_modules))

logger.info("pytest_ngsfixtures test data missing for the following MultiQC modules: {}".format(", ".join(sorted(missing_modules))))
logger.info("MultiQC parsers missing for the following pytest_ngsfixtures test fixtures: {}".format(", ".join(sorted(sorted(set(missing_multiqc_parsers))))))


name_map = {'bowtie1': 'bowtie 1', 'bowtie2': 'bowtie 2', 'fastq_screen': 'fastq screen'}
def evaluate_parser(multiqc_obj, module, command):
    """Evaluate parsed data from multiqc in some way.

    Params:
      multiqc_obj (MultiqcModule): instantiated MultiQC object that performed the parsing
      module (str): name of test module
      command (str): name of submodule/command that was run. Currently not used.
    """
    assert multiqc_obj.name.lower() == name_map.get(module, module).lower()


@pytest.mark.skipif(not has_ngsfixtures, reason="pytest-ngsfixtures not installed")
@pytest.fixture(scope="function", autouse=False, params=test_fixtures,
                ids=["{} {}:{}/{}".format(x[0], x[1], x[2], x[3]) for x in test_fixtures])
def data(request, tmpdir_factory):
    """Data fixture for test_module.

    The request parameter holds information about module name,
    command, version, sequencing configuration (single- or
    paired-end), and a dictionary of formatting strings that specify
    how fixture file names should be formatted. 
    
    pytest_ngsfixtures output files are determined based on the
    formatting strings and linked to a temporary output directory
    fixture in which the test is run.

    Params:
      request (_pytest.fixtures.FixtureRequest): a request for a
        fixture from a test or fixture function.
      tmpdir_factory (_pytest.tmpdir.TempdirFactory): a fixture for
        creating arbitrary temporary directories

    """
    module, command, version, end, fmtdict = request.param
    params = {'version': version, 'end': end}
    # Generate pytest_ngsfixtures application output names relative to
    # applications/module directory
    outputs = [fmt.format(**params) for fmt in fmtdict.values()]
    # Add applications/module prefix
    sources = [os.path.join("applications", module, output) for output in outputs]
    # Extract source basenames
    dests = [os.path.basename(src) for src in sources]
    # Generate a unique test output directory name
    fdir = os.path.join(module, str(version), command, end)
    # Make a temporary directory using unique test directory name
    pdir = factories.safe_mktemp(tmpdir_factory, fdir)
    # Symlink pytest_ngsfixtures files to temporary directory; the
    # safe_symlink function automagically uses pytest_ngsfixtures
    # installation directory to infer location of src
    for src, dst in zip(sources, dests):
        p = factories.safe_symlink(pdir, src, dst)
    # We need to remap some names due to naming inconsistencies with
    # pytest_ngsfixtures
    return remap.get(module, module), remap.get(command, command), pdir



@pytest.mark.skipif(not has_ngsfixtures, reason="pytest-ngsfixtures not installed")
def test_module(data, monkeypatch):
    """Test a MultiQC module on fixture data.

    Params:
      data (list): contains module name, command to run and test
        directory fixture
      monkeypatch (_pytest.monkeypatch.MonkeyPatch): pytest
        monkeypatch for safely setting/deleting attributes and other
        global settings
    """
    module, command, pdir = data

    def mockwrite(data, fn, sort_cols, data_format):
        return None

    # On initialization, MultiQC sets up a configuration dictionary
    # with parameters required for locating files. We have to mock two
    # of these functions for the tests to run. Most importantly, set
    # the filelist attribute to point to the fixture directory
    filelist = [{'root': p.dirname, 'fn': p.basename} for p in pdir.visit()]
    monkeypatch.setattr(multiqc.utils.report, 'files', filelist)
    monkeypatch.setattr(multiqc.utils.util_functions, 'write_data_file', mockwrite)

    mod = "multiqc.modules.{}".format(module)
    loader = importlib.find_loader(mod)
    m = loader.load_module()
    multiqc_module = getattr(m, "MultiqcModule")
    multiqc_obj = multiqc_module()
    evaluate_parser(multiqc_obj, module, command)
