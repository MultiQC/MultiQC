from collections import OrderedDict
from multiqc import config


read_format = '{:,.1f}'
if config.read_count_multiplier == 1:
    read_format = '{:,.0f}'
read_format += '&nbsp;' + config.read_count_prefix

base_format = '{:,.1f}&nbsp;'
if config.base_count_multiplier == 1:
    base_format = '{:,.0f}'
elif config.base_count_multiplier == 0.000000001:
    base_format = '{:,.2f}'
base_format += '&nbsp;' + config.base_count_prefix


class Metric:
    def __init__(self, id, title, in_genstats=None, in_own_tabl=None, unit=None,
                 descr=None, fmt=None, modify=None, namespace=None, precision=None,
                 reversed=False):
        self.id = id
        self.title = title
        self.in_genstats = in_genstats
        self.in_own_tabl = in_own_tabl
        self.unit = unit
        self.descr = descr
        self.fmt = fmt
        self.modify = modify
        self.namespace = namespace
        self.precision = precision
        self.reversed = reversed


def make_headers(parsed_metric_ids, metrics):
    # Init general stats table
    genstats_headers = OrderedDict()

    # Init headers for an own separate table
    own_tabl_headers = OrderedDict()

    for metric in metrics:
        col = dict(
            title=metric.title,
            description=metric.descr,
            min=0,
        )

        # choosing color based on metric namespace, guessing namespace by unit
        if metric.unit == 'reads':
            col['scale'] = 'Greens'
        elif metric.unit == 'bases':
            col['scale'] = 'Greys'
        elif metric.unit == 'bp':
            col['scale'] = 'Purples'
        elif metric.unit == 'x' or metric.namespace == 'Coverage':
            col['scale'] = 'Blues'
            if metric.unit == '%':
                col['scale'] = 'RdBu'
                if metric.reversed:
                    col['scale'] += '-rev'
        elif metric.namespace == 'Variants':
            col['scale'] = 'Purples'

        if metric.id + ' pct' in parsed_metric_ids:
            # if % value is available, showing it instead of the number value; the number value will be hidden
            pct_col = dict(
                col,
                description=metric.descr.replace(', {}', '').replace('Number of ', '% of '),
                max=100,
                suffix='%',
            )
            if metric.unit == 'reads':
                pct_col['scale'] = 'RdYlGn'
            elif metric.unit == 'bases':
                pct_col['scale'] = 'RdGy'
            elif metric.unit == 'x' or metric.namespace == 'Coverage':
                pct_col['scale'] = 'RdBu'
            if metric.reversed:
                pct_col['scale'] += '-rev'

            if metric.precision is not None:
                pct_col['format'] = '{:,.' + str(metric.precision) + 'f}'

            if metric.in_genstats is not None:
                show = metric.in_genstats == '%'
                genstats_headers[metric.id + ' pct'] = dict(pct_col, hidden=not show)
            if metric.in_own_tabl is not None:
                show = metric.in_own_tabl == '%'
                own_tabl_headers[metric.id + ' pct'] = dict(pct_col, hidden=not show)

        if metric.unit == 'reads':
            col['description'] = col['description'].format(config.read_count_desc)
            col['modify'] = lambda x: x * config.read_count_multiplier
            col['shared_key'] = 'read_count'
            col['format'] = read_format
        elif metric.unit == 'bases':
            col['description'] = col['description'].format(config.base_count_desc)
            col['modify'] = lambda x: x * config.base_count_multiplier
            col['shared_key'] = 'base_count'
            col['format'] = base_format
        elif metric.unit == 'len':
            col['suffix'] = ' bp'
            col['format'] = '{:,.0f}'
        elif metric.unit == 'x':
            col['suffix'] = ' x'
            col['format'] = '{:,.1f}'
        elif metric.unit == '%':
            col['suffix'] = ' %'
            col['format'] = '{:,.1f}'
        elif any(metric.descr.startswith(pref) for pref in ('Total number of ', 'The number of ', 'Number of ')):
            col['format'] = '{:,.0f}'
        if metric.precision is not None:
            col['format'] = '{:,.' + str(metric.precision) + 'f}'

        if metric.in_genstats is not None:
            genstats_headers[metric.id] = dict(col, hidden=metric.in_genstats != '#')
        if metric.in_own_tabl is not None:
            own_tabl_headers[metric.id] = dict(col, hidden=metric.in_own_tabl != '#')

    return genstats_headers, own_tabl_headers
