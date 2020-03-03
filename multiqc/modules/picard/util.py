#!/usr/bin/env python

import os
import re


def read_sample_name(line_iter, clean_fn):
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
            if 'BaseDistributionByCycle' in new_line and 'INPUT' in new_line:
                # Pull sample name from input
                fn_search = re.search(r"INPUT=?\s*(\[?[^\s]+\]?)", new_line, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip('[]'))
                    s_name = clean_fn(s_name)
                    return s_name
    except StopIteration:
        return None
