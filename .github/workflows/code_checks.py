#!/usr/bin/env python
""" Crude code-quality checks for running in CI """

from rich import print
import glob
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODULES_DIR = os.path.join(BASE_DIR, "multiqc", "modules")

# Check that add_data_source() is called for each module
print("[bold black on yellow]  Checking for function call 'self.add_data_source'  [/]")
for fn in glob.glob(os.path.join(MODULES_DIR, "*", "*.py")):
    with open(fn, "r") as fh:
        modfile = fh.read()
        if "self.find_log_files" in modfile and "self.add_data_source" not in modfile:
            print(f"Can't find 'self.add_data_source' in /{os.path.relpath(fn, MODULES_DIR)}")

# Check that add_data_source() is called for each module
print("\n\n\n[bold black on yellow]  Checking for function call 'self.write_data_file'  [/]")
for fn in glob.glob(os.path.join(MODULES_DIR, "*", "*.py")):
    with open(fn, "r") as fh:
        modfile = fh.read()
        if "self.find_log_files" in modfile and "self.write_data_file" not in modfile:
            print(f"Can't find 'self.write_data_file' in /{os.path.relpath(fn, MODULES_DIR)}")
