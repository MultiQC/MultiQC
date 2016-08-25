# -*- coding: utf-8 -*-
import csv


def find_header(pattern, lines):
    """Parse out data from VariantEval tool."""
    all_lines = iter(lines)
    # scan until section header
    for line in all_lines:
        if line.startswith(pattern):
            break

    # collect lines until empty lines
    relevant_lines = []
    for raw_line in all_lines:
        line = raw_line.strip()
        if not line:
            break
        else:
            relevant_lines.append(line)

    # parse rows into dict
    values = csv.DictReader(relevant_lines, delimiter=' ', skipinitialspace=True)
    return values
