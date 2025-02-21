import logging
from typing import Any, Dict, List, Optional, Union, cast

from importlib_metadata import EntryPoint

from multiqc import config, report
from multiqc.plots import table_object
from multiqc.plots.plotly import violin
from multiqc.plots.table_object import ColumnDict, ColumnKeyT, SectionT, TableConfig
from multiqc.types import Anchor

logger = logging.getLogger(__name__)


# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod: Optional[EntryPoint] = None


def get_template_mod() -> EntryPoint:
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    assert _template_mod is not None
    return _template_mod


def plot(
    data: Union[List[SectionT], SectionT],
    headers: Optional[Union[List[Dict[ColumnKeyT, ColumnDict]], Dict[ColumnKeyT, ColumnDict]]] = None,
    pconfig: Union[Dict[str, Any], TableConfig, None] = None,
) -> Union[str, violin.ViolinPlot]:
    """
    Helper HTML for a violin plot.
    :param data: A list of data dicts
    :param headers: A list of dicts with information
                    for the series, such as colour scales, min and
                    max values etc.
    :param pconfig: plot config dict
    :return: HTML string
    """
    pconf = cast(TableConfig, TableConfig.from_pconfig_dict(pconfig))

    if not isinstance(data, list):
        data = [data]
    if headers is not None and not isinstance(headers, list):
        headers = [headers]

    # Make datatable objects
    dts = []
    for i, d in enumerate(data):
        h = headers[i] if headers and len(headers) > i else None
        table_id = pconf.id
        table_anchor = Anchor(f"{pconf.anchor or table_id}_table")
        if len(data) > 0:
            table_anchor = Anchor(f"{table_anchor}-{i + 1}")
        table_anchor = Anchor(report.save_htmlid(table_anchor))  # make sure it's unique
        dt = table_object.DataTable.create(d, table_id, table_anchor, pconf.model_copy(), h)
        dts.append(dt)

    mod = get_template_mod()
    if "violin" in mod.__dict__ and callable(mod.__dict__["violin"]):
        # noinspection PyBroadException
        try:
            return mod.__dict__["violin"](dts)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return violin.plot(dts)
