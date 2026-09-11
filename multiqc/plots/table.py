import logging
from typing import Any, Union

from multiqc.plots.table_object import (
    ColumnDict,
    ColumnKeyT,
    SectionT,
    TableConfig,
)
from multiqc.plots.violin import ViolinPlot, ViolinPlotInputData
from multiqc.types import ColumnKey, SectionKey

logger = logging.getLogger(__name__)


def plot_with_sections(
    data: dict[SectionKey, SectionT],
    headers: dict[SectionKey, dict[ColumnKey, ColumnDict]],
    pconfig: dict[str, Any] | TableConfig | None = None,
) -> Union["ViolinPlot", str, None]:
    """
    Helper HTML for a violin plot.
    """
    inputs = ViolinPlotInputData.create_from_sections(data, headers, pconfig)
    inputs = ViolinPlotInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        return None

    return ViolinPlot.from_inputs(inputs)


def plot(
    data: SectionT,
    headers: dict[ColumnKeyT, ColumnDict] | None = None,
    pconfig: dict[str, Any] | TableConfig | None = None,
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

    assert isinstance(data, dict)
    if headers is not None:
        assert isinstance(headers, dict)
    return violin_plot(data, headers, pconfig, show_table_by_default=True)
