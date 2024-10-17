import dataclasses
import io
from enum import Enum
from typing import Generic, List, NewType, Optional, TypeVar, Union, TypedDict

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


@dataclasses.dataclass
class SampleNameMeta:
    original_name: SampleName
    trimmed_name: Optional[SampleName] = None
    trimmed_suffixes: List[str] = dataclasses.field(default_factory=list)
    group: Optional[SampleGroup] = None
    labels: List[str] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class Section:
    name: str
    anchor: Anchor
    id: SectionId  # unlike anchor, doesn't have to be different from the module or plot ids
    description: str
    module: str
    module_info: str = ""
    comment: str = ""
    helptext: str = ""
    content_before_plot: str = ""
    content: str = ""
    plot: str = ""
    print_section: bool = True
    plot_anchor: Optional[Anchor] = None
    ai_summary: str = ""
