#!/usr/bin/env python
from __future__ import print_function

import re
from collections import defaultdict


def parse_fastqc_metrics_file(f):
    """
    NA12878.fastqc_metrics.csv

    READ MEAN QUALITY,Read1,Q30 Reads,151104
    ...
    POSITIONAL BASE CONTENT,Read1,ReadPos 1 A Bases,2963930
    ...
    NUCLEOTIDE QUALITY,Read1,Q7 A Bases,777933
    ...
    READ LENGTHS,Read1,145-152bp Length Reads,10142922
    ...
    READ BASE CONTENT,Read1,0% A Reads,1997
    ...
    """
    f['s_name'] = re.search(r'(.*).fastqc_metrics.csv', f['fn']).group(1)
    r1_name = "{}_R1".format(f['s_name'])
    r2_name = "{}_R2".format(f['s_name'])

    data_by_sample = {}
    data_by_sample[r1_name] = defaultdict(lambda: defaultdict(int))
    data_by_sample[r2_name] = defaultdict(lambda: defaultdict(int))
    for line in f['f'].splitlines():
        group, mate, metric, value = line.split(',')
        try:
            value = int(value)
        except ValueError:
            pass

        # Store each value by group and by metric
        assert mate in ['Read1', 'Read2']
        if mate == "Read1":
            s_name = r1_name
        elif mate == "Read2":
            s_name = r2_name
        data_by_sample[s_name][group][metric] = value

    # Delete empty mate groups so we don't generate empty datasets
    for s_name in [r1_name, r2_name]:
        if len(data_by_sample[s_name]) == 0:
            del data_by_sample[s_name]

    return data_by_sample


def average_from_range(metric_range):
    if '+' in metric_range:
        metric_range = metric_range.replace('+', '')
    if '>=' in metric_range:
        metric_range = metric_range.replace('>=', '')
    if "-" in metric_range:
        start, end = metric_range.split('-')
        avg_pos = (int(end) + int(start)) / 2.0
    else:
        avg_pos = int(metric_range)
    return avg_pos


def average_pos_from_metric(metric):
    parts = metric.split()
    metric_range = parts[1]
    return average_from_range(metric_range)


def average_pos_from_size(metric):
    parts = metric.split()
    len_range = parts[0].split('bp')[0]
    return average_from_range(len_range)


def percentage_from_content_metric(metric):
    parts = metric.split()
    pct = int(parts[0].split("%")[0])
    return pct


def pos_qual_table_cmp(key):
    parts = key.split()
    pos = average_from_range(parts[1])
    pct = int(parts[2][:-1])
    return (pos * 1000 + pct)


def sortPosQualTableKeys(data_dict):
    return sorted(data_dict.keys(), key=pos_qual_table_cmp)
