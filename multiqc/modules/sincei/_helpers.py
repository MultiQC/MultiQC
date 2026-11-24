from itertools import groupby

import numpy as np


def group_median_by_cell_prefix(dat, val):
    """
    Group cells by the prefix of their `Cell_ID` (everything before `::`)
    and return the median of column `val` per group, cast to float.
    """
    dc = groupby(
        sorted(dat.items(), key=lambda x: x[1]["Cell_ID"].split("::")[0]),
        lambda x: x[1]["Cell_ID"].split("::")[0],
    )
    return {i: {val: np.median([float(j[1][val]) for j in j])} for i, j in dc}
