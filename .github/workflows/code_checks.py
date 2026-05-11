"""Crude code-quality checks for running in CI"""

import glob
import os
import sys

from rich import print

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODULES_DIR = os.path.join(BASE_DIR, "multiqc", "modules")
PACKAGE_DIR = os.path.join(BASE_DIR, "multiqc")
num_errors = 0

SUPER_INIT_TERMS = ("super(MultiqcModule, self).__init__(", "super().__init__(")

must_be_present_after = [
    ("self.add_data_source", "self.find_log_files"),
    ("self.write_data_file", "self.find_log_files"),
    ("doi=", SUPER_INIT_TERMS),
    ("self.add_software_version", "self.find_log_files"),
]

must_be_avoided_after = [
    ('f["contents_lines"]', "self.find_log_files", "Use 'f[\"f\"].splitlines()' instead"),
    ("super(MultiqcModule", "class MultiqcModule", "Use 'super().__init__()' instead"),
]

# Patterns that must never appear anywhere in the multiqc/ package. Each entry is
# (forbidden_substring, suggestion, exempt_relpaths). Files in exempt_relpaths are
# allowed to contain the pattern (e.g. legacy compatibility shims).
must_be_avoided_globally = [
    (
        "from __future__",
        "Drop __future__ imports — MultiQC is Python 3 only",
        (),
    ),
    (
        "OrderedDict",
        "Use a regular dict — Python 3.7+ preserves insertion order",
        ("plots/bargraph.py", "plots/box.py"),
    ),
]


def matches_search_term(contents, search_term):
    """Return True if any of the search term variants appears in contents."""
    if isinstance(search_term, tuple):
        return any(term in contents for term in search_term)
    return search_term in contents


# Check that add_data_source() is called for each module
for must_be_present, search_term in must_be_present_after:
    print(f"[bold black on yellow]  Checking for '{must_be_present}'  [/]")
    missing_files = []
    for fn in glob.glob(os.path.join(MODULES_DIR, "*", "*.py")):
        with open(fn, "r") as fh:
            contents = fh.read()
            if matches_search_term(contents, search_term) and must_be_present not in contents:
                relpath = os.path.relpath(fn, MODULES_DIR)
                missing_files.append((relpath, f"Can't find '{must_be_present}' in /{relpath}"))
                num_errors += 1

    for file in sorted(missing_files, key=lambda x: x[0]):
        print(file[1])
    print("\n\n")

for must_be_avoided, search_term, suggestion in must_be_avoided_after:
    print(f"[bold black on yellow]  Checking that '{must_be_avoided}' is not found  [/]")
    found_files = []
    for fn in glob.glob(os.path.join(MODULES_DIR, "*", "*.py")):
        with open(fn, "r") as fh:
            contents = fh.read()
            if matches_search_term(contents, search_term) and must_be_avoided in contents:
                relpath = os.path.relpath(fn, MODULES_DIR)
                message = f"Found '{must_be_avoided}' in /{relpath}"
                if suggestion:
                    message += f" ({suggestion})"
                found_files.append((relpath, message))
                num_errors += 1

    for file in sorted(found_files, key=lambda x: x[0]):
        print(file[1])
    print("\n\n")

for forbidden, suggestion, exempt in must_be_avoided_globally:
    print(f"[bold black on yellow]  Checking that '{forbidden}' is not found anywhere in multiqc/  [/]")
    found_files = []
    for fn in glob.glob(os.path.join(PACKAGE_DIR, "**", "*.py"), recursive=True):
        relpath = os.path.relpath(fn, PACKAGE_DIR)
        if relpath in exempt:
            continue
        with open(fn, "r") as fh:
            contents = fh.read()
            if forbidden in contents:
                message = f"Found '{forbidden}' in /{relpath}"
                if suggestion:
                    message += f" ({suggestion})"
                found_files.append((relpath, message))
                num_errors += 1

    for file in sorted(found_files, key=lambda x: x[0]):
        print(file[1])
    print("\n\n")

if num_errors > 0:
    print(f"[bold red]{num_errors} problems found")
    print(
        "\n\nNB: If a file shouldn't use a function for a legitimate reason, "
        "add a comment explaining why that is the case with the function name. "
        "Then the function name string will be found."
    )
    sys.exit(1)
else:
    print("\n\n[green]No problems found!")
