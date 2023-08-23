#!/usr/bin/env python

""" MultiQC lint helpers. Simple additional tests to run when
--lint is specified (outside scope of normal functions) """


import glob
import os

import yaml

from multiqc.utils import config, report

logger = config.logger


def run_tests():
    """Run all lint tests"""
    if config.lint:
        check_mods_docs_readme()


def lint_error(msg):
    """Add a lint error to the report"""
    logger.error(msg)
    report.lint_errors.append(msg)


def check_mods_docs_readme():
    """Check that all modules are listed in the YAML index
    at the top of docs/README.md"""

    docs_mods = []

    docs_dir = os.path.join(os.path.dirname(config.MULTIQC_DIR), "docs", "modules")
    if not os.path.isdir(docs_dir) and os.environ.get("GITHUB_WORKSPACE"):
        docs_dir = os.path.join(os.environ.get("GITHUB_WORKSPACE"), "docs", "modules")
    if not os.path.isdir(docs_dir):
        logger.warning(f"Can't check docs readmes in lint test as directory doesn't exist: {docs_dir}")
        return None

    for fn in glob.glob(os.path.join(docs_dir, "*.md")):
        docs_mods.append(os.path.basename(fn)[:-3])
    logger.info(f"Checking docs readmes in '{docs_dir}' as --lint specified")

    # Check that installed modules are listed in docs/modules
    for m in config.avail_modules.keys():
        if m not in docs_mods and m != "custom_content":
            lint_error(f"LINT: Module '{m}' found in installed modules, but not docs/modules")

    # Check that modules in docs/modules are installed
    for m in docs_mods:
        if m not in config.avail_modules.keys() and m != "custom_content":
            lint_error(f"LINT: Module '{m}' found in docs/modules, but not installed modules")

    # Check that all modules have a YAML header conforming to the required structure
    for fn in glob.glob(os.path.join(docs_dir, "*.md")):
        with open(fn) as fh:
            header = yaml.safe_load(fh)
            if not header:
                lint_error(
                    f"LINT: the module doc does not have a YAML header with `name`, "
                    f"`url`, and `description` fields: {fn}"
                )
            else:
                for field in ["name", "url", "description"]:
                    if field not in header:
                        lint_error(f"LINT: the YAML header does not have a '{field}' " f"field: {fn}")
