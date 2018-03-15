#!/usr/bin/env python

""" MultiQC lint helpers. Simple additional tests to run when
--lint is specified (outside scope of normal functions) """

from __future__ import print_function
import os
import yaml

from multiqc.utils import config, report
logger = config.logger

def run_tests():
    """ Run all lint tests """
    if config.lint:
        check_mods_docs_readme()


def check_mods_docs_readme():
    """ Check that all modules are listed in the YAML index
    at the top of docs/README.md """

    docs_mods = []

    readme_fn = os.path.join( os.path.dirname(config.MULTIQC_DIR), 'docs', 'README.md')
    if not os.path.isfile(readme_fn):
        if os.environ.get('TRAVIS_BUILD_DIR') is not None:
            readme_fn = os.path.join( os.environ.get('TRAVIS_BUILD_DIR'), 'docs', 'README.md')
        else:
            logger.warn("Can't check docs readme in lint test as file doesn't exist: {}".format(readme_fn))
            return None
    with open(readme_fn) as f:
        fm = next(yaml.load_all(f))

    for section in fm['MultiQC Modules']:
        for name, fn in fm['MultiQC Modules'][section].items():
            # remove modules/ and .md
            docs_mods.append(fn[8:-3])

    # Check that installed modules are listed in docs/README.md
    for m in config.avail_modules.keys():
        if m not in docs_mods and m != 'custom_content':
            errmsg = "LINT: Module '{}' found in installed modules, but not docs/README.md".format(m)
            logger.error(errmsg)
            report.lint_errors.append(errmsg)

    # Check that modules in docs/README.md are installed
    for m in docs_mods:
        if m not in config.avail_modules.keys() and m != 'custom_content':
            errmsg = "LINT: Module '{}' found in docs/README.md, but not installed modules".format(m)
            logger.error(errmsg)
            report.lint_errors.append(errmsg)
