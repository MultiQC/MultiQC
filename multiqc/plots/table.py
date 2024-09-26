import logging
from typing import Any, Dict, List, Optional, Set, Union, cast

from importlib_metadata import EntryPoint

from multiqc import config, report
from multiqc.plots import table_object
from multiqc.plots.plotly import table
from multiqc.plots.plotly.violin import ViolinPlot
from multiqc.plots.table_object import ColumnDict, ColumnKey, SectionT, TableConfig
from multiqc.types import Anchor

logger = logging.getLogger(__name__)

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod: Optional[EntryPoint] = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(
    data: Union[SectionT, List[SectionT]],
    headers: Optional[
        Union[
            List[Dict[str, ColumnDict]],
            Dict[str, ColumnDict],
        ]
    ] = None,
    pconfig: Union[Dict[str, Any], TableConfig, None] = None,
) -> "ViolinPlot":
    """Return HTML for a MultiQC table.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param headers: list of optional dicts with column config in key:value pairs.
    :param pconfig: plot config dict
    :return: HTML ready to be inserted into the page
    """
    pconf = cast(TableConfig, TableConfig.from_pconfig_dict(pconfig))

    # Make a datatable object
    table_id = pconf.id
    table_anchor = Anchor(f"{pconf.anchor or table_id}_table")
    table_anchor = Anchor(report.save_htmlid(table_anchor))  # make sure it's unique
    dt = table_object.DataTable.create(
        data, table_id=table_id, table_anchor=table_anchor, pconfig=pconf.model_copy(), headers=headers
    )

    return plot_dt(dt)


def plot_dt(dt: table_object.DataTable) -> "ViolinPlot":
    mod = get_template_mod()
    if "table" in mod.__dict__ and callable(mod.__dict__["table"]):
        # Collect unique sample names
        s_names: Set[str] = set()
        for section in dt.sections:
            for s_name in section.rows_by_sgroup.keys():
                s_names.add(s_name)

        # noinspection PyBroadException
        try:
            return mod.__dict__["table"](dt, s_names, dt.pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return table.plot([dt])
