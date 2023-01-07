import re
from collections import defaultdict

_R1 = "Read1"
_R2 = "Read2"
_VALID_MATES = [_R1, _R2]


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
    s_name = re.search(r"(.*).fastqc_metrics.csv", f["fn"]).group(1)
    data_by_sample = initialize_dataset(s_name)

    for line in f["f"].splitlines():
        # This conforms with the standard DRAGEN metrics format
        # Percentage is currently unused
        tokens = line.split(",")
        if len(tokens) == 4:
            group, mate, metric, value = line.split(",")
            percentage = None
        elif len(tokens) == 5:
            group, mate, metric, value, percentage = line.split(",")
        else:
            raise ValueError(f"Unexpected number of values in line {line}")

        try:
            value = int(value)
        except ValueError:
            pass

        # Store each value by group and by metric
        assert mate in _VALID_MATES
        data_by_sample[s_name][mate][group][metric] = value

    # Delete empty mate groups so we don't generate empty datasets
    # Altered to work on python3 which doesn't allow deleting from dict during iterating
    delete = [mate for mate, mate_data in data_by_sample[s_name].items() if len(mate_data) == 0]
    for key in delete:
        del data_by_sample[s_name][key]

    return s_name, data_by_sample


def initialize_dataset(s_name):
    data = dict()
    data[s_name] = dict()
    for mate in _VALID_MATES:
        data[s_name][mate] = defaultdict(lambda: defaultdict(int))
    return data


def average_from_range(metric_range):
    if "+" in metric_range:
        metric_range = metric_range.replace("+", "")
    if ">=" in metric_range:
        metric_range = metric_range.replace(">=", "")
    if "-" in metric_range:
        start, end = metric_range.split("-")
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
    len_range = parts[0].split("bp")[0]
    return average_from_range(len_range)


def percentage_from_content_metric(metric):
    parts = metric.split()
    pct = int(parts[0].split("%")[0])
    return pct


def pos_qual_table_cmp(key):
    parts = key.split()
    pos = average_from_range(parts[1])
    pct = int(parts[2][:-1])
    return pos * 1000 + pct


def sortPosQualTableKeys(data_dict):
    return sorted(data_dict.keys(), key=pos_qual_table_cmp)
