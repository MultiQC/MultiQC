""" MultiQC lint helpers. Simple additional tests to run when
--strict is specified (outside scope of normal functions) """


import glob
import os

import yaml

from multiqc.utils import config, report

logger = config.logger


def run_tests():
    """Run all lint tests"""
    if config.strict:
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
    logger.info(f"Checking docs readmes in '{docs_dir}' as --strict specified")

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
            # Load the YAML header from the markdown file. YAML header should be placed between --- and --- in the beginning of the file
            try:
                header = fh.read().split("---")[1]
            except IndexError:
                lint_error(f"LINT: '{fn}' doesn't have a YAML header between '---'")
                continue
            try:
                header = yaml.safe_load(header)
            except yaml.YAMLError as e:
                lint_error(f"LINT: '{fn}' contains an incorrectly formatted YAML header: {e}")
                continue
            if header is None:
                lint_error(f"LINT: '{fn}' contains an empty YAML header")
                continue
            req_fields = ["name", "url", "description"]
            for field in req_fields:
                if field not in header:
                    lint_error(
                        f"LINT: the YAML header in '{fn}' does not have a '{field}' field. Required fields: {req_fields}"
                    )
