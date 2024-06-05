"""MultiQC functions to plot a scatter plot"""

import logging
from typing import Union, Dict

from multiqc import config
from multiqc.plots.plotly import scatter
from multiqc.plots.plotly.scatter import ScatterConfig

logger = logging.getLogger(__name__)

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(
    data,
    pconfig: Union[Dict, ScatterConfig, None] = None,
) -> Union[scatter.ScatterPlot, str]:
    """Plot a scatter plot with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    assert pconfig is not None, "pconfig must be provided"
    if isinstance(pconfig, dict):
        pconfig = ScatterConfig(**pconfig)

    # Given one dataset - turn it into a list
    if not isinstance(data, list):
        data = [data]

    plotdata = list()
    for data_index, ds in enumerate(data):
        d = list()
        for s_name in ds:
            # Ensure any overwriting conditionals from data_labels (e.g. ymax) are taken in consideration
            series_config: ScatterConfig = pconfig.model_copy()
            if pconfig.data_labels and isinstance(pconfig.data_labels[data_index], dict):
                # if not a dict: only dataset name is provided
                for k, v in pconfig.data_labels[data_index].items():
                    if k in series_config.model_fields:
                        setattr(series_config, k, v)

            if not isinstance(ds[s_name], list):
                ds[s_name] = [ds[s_name]]
            for point in ds[s_name]:
                if point["x"] is not None:
                    if series_config.xmax is not None and float(point["x"]) > float(series_config.xmax):
                        continue
                    if series_config.xmin is not None and float(point["x"]) < float(series_config.xmin):
                        continue
                if point["y"] is not None:
                    if series_config.ymax is not None and float(point["y"]) > float(series_config.ymax):
                        continue
                    if series_config.ymin is not None and float(point["y"]) < float(series_config.ymin):
                        continue
                if "name" in point:
                    point["name"] = f'{s_name}: {point["name"]}'
                else:
                    point["name"] = s_name

                for k in ["color", "opacity", "marker_size", "marker_line_width"]:
                    if k not in point:
                        v = getattr(series_config, k)
                        if v is not None:
                            if isinstance(v, dict) and s_name in v:
                                point[k] = v[s_name]
                            else:
                                point[k] = v
                d.append(point)
        plotdata.append(d)

    if pconfig.square:
        if pconfig.ymax is None and pconfig.xmax is None:
            # Find the max value
            max_val = 0
            for d in plotdata:
                for s in d:
                    max_val = max(max_val, s["x"], s["y"])
            max_val = 1.02 * max_val  # add 2% padding
            pconfig.xmax = pconfig.xmax if pconfig.xmax is not None else max_val
            pconfig.ymax = pconfig.ymax if pconfig.ymax is not None else max_val

    # Add on annotation data series
    # noinspection PyBroadException
    try:
        if pconfig.extra_series:
            extra_series = pconfig.extra_series
            if isinstance(pconfig.extra_series, dict):
                extra_series = [[pconfig.extra_series]]
            elif isinstance(pconfig.extra_series, list) and isinstance(pconfig.extra_series[0], dict):
                extra_series = [pconfig.extra_series]
            for i, es in enumerate(extra_series):
                for s in es:
                    plotdata[i].append(s)
    except Exception:
        pass

    # Make a plot
    mod = get_template_mod()
    if "scatter" in mod.__dict__ and callable(mod.scatter):
        # noinspection PyBroadException
        try:
            return mod.scatter(plotdata, pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return scatter.plot(plotdata, pconfig)
