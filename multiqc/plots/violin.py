import logging
import random
import string
from typing import List, Dict, Optional, Union

from multiqc import config, report
from multiqc.plots import table_object
from multiqc.plots.plotly import violin
from multiqc.plots.table_object import TableConfig, InputSectionT, InputHeaderT
from multiqc.types import AnchorT

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
    data: Union[List[InputSectionT], InputSectionT],
    headers: Optional[Union[List[InputHeaderT], InputHeaderT]] = None,
    pconfig: Union[Dict, TableConfig, None] = None,
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
    assert pconfig is not None, "pconfig must be provided"
    if isinstance(pconfig, dict):
        pconfig = TableConfig(**pconfig)

    if not isinstance(data, list):
        data = [data]
    if headers is not None and not isinstance(headers, list):
        headers = [headers]

    # Make datatable objects
    dts = []
    for i, d in enumerate(data):
        h = headers[i] if headers and len(headers) > i else None
        table_id = pconfig.id
        table_anchor = AnchorT(f"{pconfig.anchor or table_id}-table")
        if len(data) > 0:
            table_anchor = AnchorT(f"{table_anchor}-{i + 1}")
        table_anchor = report.save_htmlid(table_anchor)  # make sure it's unique
        dt = table_object.DataTable.create(d, table_id, table_anchor, pconfig.model_copy(), h)
        dts.append(dt)

    mod = get_template_mod()
    if "violin" in mod.__dict__ and callable(mod.violin):
        # noinspection PyBroadException
        try:
            return mod.violin(dts)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return violin.plot(dts)
