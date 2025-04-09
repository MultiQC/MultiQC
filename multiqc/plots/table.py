import logging
from typing import Any, Dict, Optional, Union


from multiqc.plots.table_object import (
    ColumnDict,
    ColumnKeyT,
    SectionT,
    TableConfig,
)
from multiqc.plots.violin import ViolinPlotInputData, ViolinPlot
from multiqc.types import ColumnKey, SectionKey

logger = logging.getLogger(__name__)


def plot_with_sections(
    data: Dict[SectionKey, SectionT],
    headers: Dict[SectionKey, Dict[ColumnKey, ColumnDict]],
    pconfig: Union[Dict[str, Any], TableConfig, None] = None,
) -> Union["ViolinPlot", str]:
    """
    Helper HTML for a violin plot.
    """
    inputs = ViolinPlotInputData.create_from_sections(data, headers, pconfig)
    inputs = ViolinPlotInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        logger.warning(f"Tried to make table/violin plot, but had no data. pconfig: {pconfig}")
        return '<p class="text-danger">Error - was not able to plot data.</p>'

    return ViolinPlot.create(
        inputs.dts,
        anchor=inputs.anchor,
        show_table_by_default=True,
    )


def plot(
    data: SectionT,
    headers: Optional[Dict[ColumnKeyT, ColumnDict]] = None,
    pconfig: Union[Dict[str, Any], TableConfig, None] = None,
) -> Union["ViolinPlot", str, None]:
    """Return HTML for a MultiQC table.
    :param data: A list of data dicts
    :param headers: A list of dicts with information
                    for the series, such as colour scales, min and
                    max values etc.
    :param pconfig: plot config dict
    :return: plot object
    """
    from multiqc.plots.violin import plot as violin_plot

    return violin_plot(data, headers, pconfig, show_table_by_default=True)
