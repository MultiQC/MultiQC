#!/usr/bin/env python

import os
import re
from collections import OrderedDict


def read_sample_name(line_iter, clean_fn, program_name):
    """
    Consumes lines from the provided line_iter and parses those lines
    as a header for the picard base distribution file.  The header
    file is assumed to contain a line with both 'INPUT' and
    'BaseDistributionByCycle'.

    If the header parses correctly, the sample name is returned.  If
    the header does not parse correctly, None is returned.
    """
    try:
        while True:
            new_line = next(line_iter)
            new_line = new_line.strip()
            if program_name in new_line and "INPUT" in new_line:
                # Pull sample name from input
                fn_search = re.search(r"INPUT=?\s*(\[?[^\s]+\]?)", new_line, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip("[]"))
                    s_name = clean_fn(s_name)
                    return s_name
    except StopIteration:
        return None


def read_histogram(self, program_key, program_name, headers, formats):
    """
    Reads a Picard HISTOGRAM file.

    Args:
        self: the Picard QC module
        program_key: the key used to find the program (ex. picard/quality_by_cycle)
        program_name: the program key in the header to find the I/INPUT line
        headers: the list of expected headers for the histogram
        formats: the list of methods to apply to re-format each field (on a given row)
    """
    all_data = OrderedDict()

    assert len(formats) == len(headers)

    # Go through logs and find Metrics
    for f in self.find_log_files(program_key, filehandles=True):
        lines = iter(f["f"])

        # read through the header of the file to obtain the
        # sample name
        clean_fn = lambda n: self.clean_s_name(n, f)
        s_name = read_sample_name(lines, clean_fn, program_name)
        if s_name is None:
            continue

        sample_data = OrderedDict()

        try:
            # skip to the histogram
            line = next(lines)
            while not line.startswith("## HISTOGRAM"):
                line = next(lines)

            # check the header
            line = next(lines)
            if headers != line.strip().split("\t"):
                continue

            # slurp the data
            line = next(lines).rstrip()
            while line:
                fields = line.split("\t")
                assert len(fields) == len(headers)
                for i in range(len(fields)):
                    fields[i] = formats[i](fields[i])

                sample_data[fields[0]] = OrderedDict(zip(headers, fields))
                line = next(lines).rstrip()

        except StopIteration:
            pass

        # append the data
        if sample_data:
            all_data[s_name] = sample_data

    return self.ignore_samples(all_data)
