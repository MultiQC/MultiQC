import io
from enum import Enum
from typing import Generic, NewType, Optional, TypeVar, Union
from typing_extensions import TypedDict

Anchor = NewType("Anchor", str)
ModuleId = NewType("ModuleId", str)
SectionId = NewType("SectionId", str)

ColumnKey = NewType("ColumnKey", str)
SampleName = NewType("SampleName", str)
SampleGroup = NewType("SampleGroup", str)


class FileDict(TypedDict):
    fn: str
    root: str
    sp_key: str


FT = TypeVar("FT", bound=Union[str, io.IOBase, None])


class LoadedFileDict(TypedDict, Generic[FT]):
    fn: str
    root: str
    sp_key: str
    s_name: str
    f: FT


class PlotType(Enum):
    """
    Plot labels used in custom content, as well as in JS to load plot into Plotly-js
    """

    BAR = "bar_graph"
    LINE = "xy_line"
    BOX = "box"
    SCATTER = "scatter"
    VIOLIN = "violin"
    TABLE = "table"
    HEATMAP = "heatmap"
    HTML = "html"
    IMAGE = "image"
    GENERALSTATS = "generalstats"

    @staticmethod
    def from_str(val: Union[str, None, "PlotType"]) -> Optional["PlotType"]:
        if val is None:
            return None
        if isinstance(val, PlotType):
            return val
        if val in ["bar", "bargraph", "bar_graph"]:
            return PlotType.BAR
        elif val in ["line", "linegraph", "xy_line"]:
            return PlotType.LINE
        elif val in ["box", "boxplot", "box_plot"]:
            return PlotType.BOX
        elif val in ["scatter", "scatterplot"]:
            return PlotType.SCATTER
        elif val in ["violin", "beeswarm", "violin_plot"]:
            return PlotType.VIOLIN
        elif val in ["heatmap"]:
            return PlotType.HEATMAP
        elif val in ["table"]:
            return PlotType.TABLE
        elif val in ["html"]:
            return PlotType.HTML
        elif val in ["image"]:
            return PlotType.IMAGE
        elif val in ["generalstats"]:
            return PlotType.GENERALSTATS
        return None
