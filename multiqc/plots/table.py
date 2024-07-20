import logging
from typing import List, Dict, Union, Optional, Sequence

from multiqc.plots import table_object
from multiqc.plots.plotly.plot import Plot
from multiqc import config
from multiqc.plots.plotly import table
from multiqc.plots.table_object import TableConfig, DatasetT

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
    data: Union[DatasetT, Sequence[DatasetT]],
    headers: Optional[Union[List[Dict], Dict]] = None,
    pconfig: Union[Dict, TableConfig, None] = None,
) -> Union[str, Plot]:
    """Return HTML for a MultiQC table.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param headers: list of optional dicts with column config in key:value pairs.
    :param pconfig: plot config dict
    :return: HTML ready to be inserted into the page
    """
    assert pconfig is not None, "pconfig must be provided"
    if isinstance(pconfig, dict):
        pconfig = TableConfig(**pconfig)

    # Make a datatable object
    dt = table_object.DataTable.create(data, pconfig, headers)

    return plot_dt(dt)


def plot_dt(dt: table_object.DataTable) -> Union[str, Plot]:
    mod = get_template_mod()
    if "table" in mod.__dict__ and callable(mod.table):
        # Collect unique sample names
        s_names = set()
        for d in dt.raw_data:
            for s_name in d.keys():
                s_names.add(s_name)

        # noinspection PyBroadException
        try:
            return mod.table(dt, s_names, dt.pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return table.plot([dt])
