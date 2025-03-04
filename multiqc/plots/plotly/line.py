import io
import logging
import os
import random
from typing import Any, Dict, Generic, List, Literal, Mapping, Optional, Tuple, Type, TypeVar, Union

import plotly.graph_objects as go  # type: ignore
from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.plots.plotly.plot import BaseDataset, PConfig, Plot, PlotType, convert_dash_style
from multiqc.types import SampleName
from multiqc.utils.util_functions import update_dict
from multiqc.validation import ValidatedConfig, add_validation_warning

logger = logging.getLogger(__name__)


KeyT = TypeVar("KeyT", int, str, float)
ValT = TypeVar("ValT", int, str, float, None)
XToYDictT = Mapping[KeyT, ValT]
DatasetT = Mapping[Union[str, SampleName], XToYDictT[KeyT, ValT]]


class Marker(ValidatedConfig):
    symbol: Optional[str] = None
    color: Optional[str] = None
    line_color: Optional[str] = None
    fill_color: Optional[str] = None
    width: int = 1

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("Marker",), **data)


class Series(ValidatedConfig, Generic[KeyT, ValT]):
    name: str = Field(default_factory=lambda: f"series-{random.randint(1000000, 9999999)}")
    pairs: List[Tuple[KeyT, ValT]]
    color: Optional[str] = None
    width: int = 2
    dash: Optional[str] = None
    showlegend: bool = True
    marker: Optional[Marker] = None

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        path_in_cfg = path_in_cfg or ("Series",)

        if "dashStyle" in data:
            add_validation_warning(path_in_cfg, "'dashStyle' field is deprecated. Please use 'dash' instead")
            data["dash"] = data.pop("dashStyle")

        tuples: List[Tuple[KeyT, ValT]] = []
        if "data" in data:
            add_validation_warning(path_in_cfg + ("data",), "'data' field is deprecated. Please use 'pairs' instead")
        for p in data.pop("data") if "data" in data else data.get("pairs", []):
            if isinstance(p, list):
                tuples.append(tuple(p))
            else:
                tuples.append(p)
        data["pairs"] = tuples

        super().__init__(**data, path_in_cfg=path_in_cfg)

        if self.dash is not None:
            self.dash = convert_dash_style(self.dash, path_in_cfg=path_in_cfg + ("dash",))

    def get_x_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        xs = [x[0] for x in self.pairs]
        if len(xs) > 0:
            return min(xs), max(xs)  # type: ignore
        return None, None

    def get_y_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        ys = [x[1] for x in self.pairs if x[1] is not None]
        if len(ys) > 0:
            return min(ys), max(ys)  # type: ignore
        return None, None


SeriesT = Union[Series, Dict[str, Any]]


class LinePlotConfig(PConfig):
    xlab: Optional[str] = None
    ylab: Optional[str] = None
    categories: bool = False
    smooth_points: Optional[int] = 500
    smooth_points_sumcounts: Union[bool, List[bool], None] = None
    extra_series: Optional[Union[Series, List[Series], List[List[Series]]]] = None
    style: Optional[Literal["lines", "lines+markers"]] = None
    hide_zero_cats: Optional[bool] = Field(False, deprecated="hide_empty")
    hide_empty: bool = False
    colors: Dict[str, str] = {}

    @classmethod
    def parse_extra_series(
        cls,
        data: Union[SeriesT, List[SeriesT], List[List[SeriesT]]],
        path_in_cfg: Tuple[str, ...],
    ) -> Union[Series, List[Series], List[List[Series]]]:
        if isinstance(data, list):
            if isinstance(data[0], list):
                return [[Series(path_in_cfg=path_in_cfg, **d) if isinstance(d, dict) else d for d in ds] for ds in data]  # type: ignore
            return [Series(path_in_cfg=path_in_cfg, **d) if isinstance(d, dict) else d for d in data]  # type: ignore
        return Series(path_in_cfg=path_in_cfg, **data) if isinstance(data, dict) else data  # type: ignore

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("lineplot",), **data)


def plot(
    lists_of_lines: List[List[Series[KeyT, ValT]]], pconfig: LinePlotConfig, sample_names: List[SampleName]
) -> "LinePlot[KeyT, ValT]":
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param lists_of_lines: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """

    # if self.n_samples >= config.max_table_rows:
    # Get a line of median values at each point with an interval of max and
    # min values
    # Create a violin of median values in each sample, showing dots for outliers
    # Clicking on a dot of a violin will show the line plot for that sample

    return LinePlot.create(pconfig, lists_of_lines, sample_names)


class Dataset(BaseDataset, Generic[KeyT, ValT]):
    lines: List[Series[KeyT, ValT]]

    def get_x_range(self) -> Tuple[Optional[KeyT], Optional[KeyT]]:
        if not self.lines:
            return None, None
        xmax, xmin = None, None
        for line in self.lines:
            _xmin, _xmax = line.get_x_range()
            if _xmin is not None:
                xmin = min(xmin, _xmin) if xmin is not None else _xmin  # type: ignore
            if _xmax is not None:
                xmax = max(xmax, _xmax) if xmax is not None else _xmax  # type: ignore
        return xmin, xmax

    def get_y_range(self) -> Tuple[Optional[ValT], Optional[ValT]]:
        if not self.lines:
            return None, None
        ymax, ymin = None, None
        for line in self.lines:
            _ymin, _ymax = line.get_y_range()
            if _ymin is not None:
                ymin = min(ymin, _ymin) if ymin is not None else _ymin  # type: ignore
            if _ymax is not None:
                ymax = max(ymax, _ymax) if ymax is not None else _ymax  # type: ignore
        return ymin, ymax

    @staticmethod
    def create(
        base_dataset: BaseDataset,
        lines: List[Series[KeyT, ValT]],
        pconfig: LinePlotConfig,
    ) -> "Dataset[KeyT, ValT]":
        dataset: Dataset[KeyT, ValT] = Dataset(**base_dataset.model_dump(), lines=lines)

        # Prevent Plotly-JS from parsing strings as numbers
        if pconfig.categories or dataset.dconfig.get("categories"):
            dataset.layout["xaxis"]["type"] = "category"

        if pconfig.style is not None:
            mode = pconfig.style
        else:
            num_data_points = sum(len(x.pairs) for x in lines)
            if num_data_points < config.lineplot_number_of_points_to_hide_markers:
                mode = "lines+markers"
            else:
                mode = "lines"

        dataset.trace_params.update(
            mode=mode,
            line={"width": 2},
        )
        if mode == "lines+markers":
            dataset.trace_params.update(
                line={"width": 0.6},
                marker={"size": 5},
            )
        return dataset

    def create_figure(
        self,
        layout: go.Layout,
        is_log: bool = False,
        is_pct: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        if layout.showlegend is True:
            # Extra space for legend
            if hasattr(layout, "height") and isinstance(layout.height, int):
                layout.height += len(self.lines) * 5

        fig = go.Figure(layout=layout)
        for series in self.lines:
            xs = [x[0] for x in series.pairs]
            ys = [x[1] for x in series.pairs]
            params: Dict[str, Any] = {
                "showlegend": series.showlegend,
                "line": {
                    "color": series.color,
                    "dash": series.dash,
                    "width": series.width,
                },
            }
            if series.marker:
                params["mode"] = "lines+markers"
                params["marker"] = {
                    "symbol": series.marker.symbol,
                    "color": series.marker.fill_color or series.marker.color or series.color,
                    "line": {
                        "width": series.marker.width,
                        "color": series.marker.line_color or series.marker.color or "black",
                    },
                }
            params = update_dict(params, self.trace_params, none_only=True)
            if len(series.pairs) == 1:
                params["mode"] = "lines+markers"  # otherwise it's invisible

            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    name=series.name,
                    text=[series.name] * len(xs),
                    **params,
                )
            )
        return fig

    def save_data_file(self) -> None:
        y_by_x_by_sample: Dict[str, Dict[Union[float, str], Any]] = dict()
        last_cats = None
        shared_cats = True
        for series in self.lines:
            y_by_x_by_sample[series.name] = dict()

            # Check to see if all categories are the same
            if len(series.pairs) > 0 and isinstance(series.pairs[0], list):
                if last_cats is None:
                    last_cats = [x[0] for x in series.pairs]
                elif last_cats != [x[0] for x in series.pairs]:
                    shared_cats = False

            for i, x in enumerate(series.pairs):
                if isinstance(x, list):
                    y_by_x_by_sample[series.name][x[0]] = x[1]
                else:
                    try:
                        y_by_x_by_sample[series.name][self.dconfig["categories"][i]] = x
                    except (ValueError, KeyError, IndexError):
                        y_by_x_by_sample[series.name][str(i)] = x

        # Custom tsv output if the x-axis varies
        if not shared_cats and config.data_format in ["tsv", "csv"]:
            sep = "\t" if config.data_format == "tsv" else ","
            fout = ""
            for series in self.lines:
                fout += series.name + sep + "X" + sep + sep.join([str(x[0]) for x in series.pairs]) + "\n"
                fout += series.name + sep + "Y" + sep + sep.join([str(x[1]) for x in series.pairs]) + "\n"

            fn = f"{self.uid}.{config.data_format_extensions[config.data_format]}"
            fpath = os.path.join(report.data_tmp_dir(), fn)
            with io.open(fpath, "w", encoding="utf-8") as f:
                f.write(fout.encode("utf-8", "ignore").decode("utf-8"))
        else:
            report.write_data_file(y_by_x_by_sample, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: PConfig, keep_hidden: bool = True) -> str:
        xsuffix = self.layout.get("xaxis", {}).get("ticksuffix", "")
        ysuffix = self.layout.get("yaxis", {}).get("ticksuffix", "")

        # Use pseudonyms for sample names if available
        pseudonyms = [report.anonymize_sample_name(series.name) for series in self.lines]

        # Create header with axis information and common suffixes
        result = "Samples: " + ", ".join(pseudonyms) + "\n\n"

        # If all y-values have the same suffix (like %), mention it in the header
        if ysuffix:
            result += f"Y values are in {ysuffix}\n\n"
        if xsuffix:
            result += f"X values are in {xsuffix}\n\n"

        for pseudonym, series in zip(pseudonyms, self.lines):
            # For other plots or fewer points, use the original format but without redundant suffixes
            pairs = [f"{self.fmt_value_for_llm(x[0])}: {self.fmt_value_for_llm(x[1])}" for x in series.pairs]
            result += f"{pseudonym} {', '.join(pairs)}\n\n"

        return result


class LinePlot(Plot[Dataset[KeyT, ValT], LinePlotConfig], Generic[KeyT, ValT]):
    datasets: List[Dataset[KeyT, ValT]]
    sample_names: List[SampleName]

    def samples_names(self) -> List[SampleName]:
        return self.sample_names

    def _plot_ai_header(self) -> str:
        result = super()._plot_ai_header()
        if self.pconfig.xlab:
            result += f"X axis: {self.pconfig.xlab}\n"
        if self.pconfig.ylab:
            result += f"Y axis: {self.pconfig.ylab}\n"
        return result

    @staticmethod
    def create(
        pconfig: LinePlotConfig,
        lists_of_lines: List[List[Series[KeyT, ValT]]],
        sample_names: List[SampleName],
    ) -> "LinePlot[KeyT, ValT]":
        n_samples_per_dataset = [len(x) for x in lists_of_lines]

        model: Plot[Dataset[KeyT, ValT], LinePlotConfig] = Plot.initialize(
            plot_type=PlotType.LINE,
            pconfig=pconfig,
            n_samples_per_dataset=n_samples_per_dataset,
            axis_controlled_by_switches=["yaxis"],
            default_tt_label="<br>%{x}: %{y}",
        )

        # Very large legend for automatically enabled flat plot mode is not very helpful
        max_n_samples = max(len(x) for x in lists_of_lines) if len(lists_of_lines) > 0 else 0
        if pconfig.showlegend is None and max_n_samples > 250:
            model.layout.showlegend = False

        model.datasets = [Dataset.create(d, lines, pconfig) for d, lines in zip(model.datasets, lists_of_lines)]

        # Make a tooltip always show on hover over any point on plot
        model.layout.hoverdistance = -1

        return LinePlot(**model.__dict__, sample_names=sample_names)


def remove_nones_and_empty_dicts(d: Mapping[Any, Any]) -> Dict[Any, Any]:
    """Remove None and empty dicts from a dict recursively."""
    return {k: remove_nones_and_empty_dicts(v) for k, v in d.items() if v is not None and v != {}}
