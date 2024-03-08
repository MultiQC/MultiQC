""" Crude code-quality checks for running in CI """

import glob
import os
import sys

from rich import print

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODULES_DIR = os.path.join(BASE_DIR, "multiqc", "modules")
num_errors = 0

checks = {
    "self.add_data_source": "self.find_log_files",
    "self.write_data_file": "self.find_log_files",
    "doi=": "super(MultiqcModule, self).__init__(",
    "self.add_software_version": "self.find_log_files",
}

# Check that add_data_source() is called for each module
for function, search_term in checks.items():
    print(f"[bold black on yellow]  Checking for function call '{function}'  [/]")
    missing_files = []
    for fn in glob.glob(os.path.join(MODULES_DIR, "*", "*.py")):
        with open(fn, "r") as fh:
            modfile = fh.read()
            if search_term in modfile and function not in modfile:
                relpath = os.path.relpath(fn, MODULES_DIR)
                missing_files.append((relpath, f"Can't find '{function}' in /{relpath}"))
                num_errors += 1

    for file in sorted(missing_files, key=lambda x: x[0]):
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
