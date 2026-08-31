"""MultiQC functions to plot a per-base sequence content heatmap"""

import logging
import math
import re
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union, cast

import numpy as np
import plotly.graph_objects as go  # type: ignore
import polars as pl
from natsort import natsorted
from pydantic import BaseModel

from multiqc import report
from multiqc.plots.plot import (
    BaseDataset,
    NormalizedPlotInputData,
    PConfig,
    Plot,
    PlotType,
    plot_anchor,
)
from multiqc.types import Anchor, SampleName

logger = logging.getLogger(__name__)


# Bin labels are FastQC's "base" field: a single position ("7") or an inclusive
# range ("10-14").
_BIN_LABEL_RE = re.compile(r"^(\d+)(?:-(\d+))?$")

_BASES: Tuple[str, str, str, str] = ("a", "c", "g", "t")


def _parse_bin_label(label: str) -> Tuple[int, int]:
    """
    Parse a sequence-content bin label into a (start, end) 1-based, inclusive range.
    "7" -> (7, 7); "10-14" -> (10, 14). Crashes loudly on anything else.
    """
    m = _BIN_LABEL_RE.match(label.strip())
    if not m:
        raise ValueError(f"Cannot parse seqcontent bin label: {label!r}")
    start = int(m.group(1))
    end = int(m.group(2)) if m.group(2) else start
    return start, end


def _make_bin(label: str, start: int, end: int, values: Mapping[str, float]) -> "SeqContentBin":
    """
    Build a normalized SeqContentBin from raw a/c/g/t values, applying the golden
    parity rules (see BUILD_PLAN.md section 1.4, preserved from
    multiqc_fastqc.js:181-194 and fastqc.py:754-760):

    1. NaN or missing base -> 0.0.
    2. Very old FastQC gives counts instead of percentages: if t > 100, renormalize
       all four values to percentages of their sum.
    3. Round to 2 decimal places.
    """
    raw: Dict[str, float] = {}
    for base in _BASES:
        if base not in values:
            raw[base] = 0.0
        else:
            v = float(values[base])
            raw[base] = 0.0 if math.isnan(v) else v

    if raw["t"] > 100:
        total = raw["t"] + raw["a"] + raw["c"] + raw["g"]
        if total > 0:
            for base in _BASES:
                raw[base] = raw[base] / total * 100

    return SeqContentBin(
        label=label,
        start=start,
        end=end,
        a=round(raw["a"], 2),
        c=round(raw["c"], 2),
        g=round(raw["g"], 2),
        t=round(raw["t"], 2),
    )


def bin_rgb(b: "SeqContentBin") -> Tuple[int, int, int]:
    """
    Golden RGB contract, must stay identical to its JS twins in
    templates/default/src/js/plots/seqcontent.js and
    templates/echarts/src/js/plots/seqcontent.js / plots/echarts/seqcontent.py:
    R = %T, G = %A, B = %C (%G implied by the complement of the other three).
    """
    r = round(b.t / 100 * 255)
    g = round(b.a / 100 * 255)
    bl = round(b.c / 100 * 255)
    return r, g, bl


def _static_scaleratio(layout: go.Layout, x_range: float, y_range: float) -> Optional[float]:
    """
    Static-export (kaleido/flat) twin of templates/default/src/js/plots/seqcontent.js
    ::fixAspectRatio(). A go.Image trace forces an equal-aspect y axis via
    yaxis.scaleanchor, which cannot be cleared; yaxis.scaleratio is the supported
    knob to distort that forced aspect so the image stretches to fill the plot area
    instead of collapsing to a thin strip. The JS twin computes this from the live
    rendered size (_fullLayout._size); here, ahead of render, we compute the
    equivalent from the export's configured height/width minus its margins, using
    the same formula: scaleratio = (plot_area_h / plot_area_w) * (x_range / y_range).
    """
    height = layout.height
    width = layout.width
    if not height or not width or not x_range or not y_range:
        return None
    margin = layout.margin
    plot_area_h = height - (margin.t or 0) - (margin.b or 0)
    plot_area_w = width - (margin.l or 0) - (margin.r or 0)
    if plot_area_h <= 0 or plot_area_w <= 0:
        return None
    return (plot_area_h / plot_area_w) * (x_range / y_range)


class SeqContentConfig(PConfig):
    """Configuration for a per-base sequence content plot"""

    xlab: str = "Position (bp)"

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data: Any):
        super().__init__(path_in_cfg=path_in_cfg or ("seqcontent",), **data)


class SeqContentBin(BaseModel):
    """One bin (a single position, or an inclusive range of positions) of
    per-base sequence content. a/c/g/t are percentages in [0, 100], already
    normalized and rounded."""

    label: str
    start: int
    end: int
    a: float
    c: float
    g: float
    t: float


class SeqContentNormalizedInputData(NormalizedPlotInputData[SeqContentConfig]):
    """
    Represents normalized input data for a seqcontent plot.
    """

    data: Dict[str, List[SeqContentBin]]

    def is_empty(self) -> bool:
        return len(self.data) == 0 or all(len(bins) == 0 for bins in self.data.values())

    def to_df(self) -> pl.DataFrame:
        """
        Convert the seqcontent data to a long-format polars DataFrame for storage
        and reloading: one row per (sample, bin).
        """
        records: List[Dict[str, Any]] = []
        for s_name, bins in self.data.items():
            for idx, b in enumerate(bins):
                records.append(
                    {
                        "sample": s_name,
                        "bin_idx": idx,
                        "label": b.label,
                        "start": b.start,
                        "end": b.end,
                        "a": b.a,
                        "c": b.c,
                        "g": b.g,
                        "t": b.t,
                    }
                )
        df = pl.DataFrame(records)
        return self.finalize_df(df)

    @classmethod
    def from_df(
        cls, df: pl.DataFrame, pconfig: Union[Dict[str, Any], SeqContentConfig], anchor: Anchor
    ) -> "SeqContentNormalizedInputData":
        """
        Create a SeqContentNormalizedInputData object from a polars DataFrame.
        """
        if "anchor" in df.columns:
            df = df.filter(pl.col("anchor") == str(anchor))

        if df.is_empty():
            pconf = (
                pconfig
                if isinstance(pconfig, SeqContentConfig)
                else cast(SeqContentConfig, SeqContentConfig.from_pconfig_dict(pconfig))
            )
            return cls(
                anchor=anchor,
                data={},
                pconfig=pconf,
                plot_type=PlotType.SEQCONTENT,
                creation_date=report.creation_date,
            )

        pconf = cast(SeqContentConfig, SeqContentConfig.from_df(df))

        data: Dict[str, List[SeqContentBin]] = {}
        for r in df.sort(["sample", "bin_idx"]).iter_rows(named=True):
            s_name = str(r["sample"])
            data.setdefault(s_name, []).append(
                SeqContentBin(
                    label=str(r["label"]),
                    start=int(r["start"]),
                    end=int(r["end"]),
                    a=float(r["a"]),
                    c=float(r["c"]),
                    g=float(r["g"]),
                    t=float(r["t"]),
                )
            )

        return cls(
            anchor=anchor,
            data=data,
            pconfig=pconf,
            plot_type=PlotType.SEQCONTENT,
            creation_date=cls.creation_date_from_df(df),
        )

    @classmethod
    def merge(
        cls, old_data: "SeqContentNormalizedInputData", new_data: "SeqContentNormalizedInputData"
    ) -> "SeqContentNormalizedInputData":
        """
        Merge two SeqContentNormalizedInputData objects: union of samples, with
        the new run winning per-sample (whole-sample granularity).
        """
        if old_data.is_empty():
            return new_data
        if new_data.is_empty():
            return old_data

        merged_data: Dict[str, List[SeqContentBin]] = dict(old_data.data)
        merged_data.update(new_data.data)

        return SeqContentNormalizedInputData(
            anchor=new_data.anchor,
            data=merged_data,
            pconfig=new_data.pconfig,
            plot_type=PlotType.SEQCONTENT,
            creation_date=report.creation_date,
        )

    @staticmethod
    def create(
        data_by_sample: Mapping[str, Mapping[str, Mapping[str, float]]],
        pconfig: Union[Dict[str, Any], SeqContentConfig, None] = None,
    ) -> "SeqContentNormalizedInputData":
        pconf = cast(SeqContentConfig, SeqContentConfig.from_pconfig_dict(pconfig))

        data: Dict[str, List[SeqContentBin]] = {}
        for s_name in natsorted(data_by_sample.keys()):
            bins: List[SeqContentBin] = []
            for label, values in data_by_sample[s_name].items():
                start, end = _parse_bin_label(str(label))
                bins.append(_make_bin(str(label), start, end, values))
            bins.sort(key=lambda b: b.start)
            data[str(s_name)] = bins

        return SeqContentNormalizedInputData(
            anchor=plot_anchor(pconf),
            plot_type=PlotType.SEQCONTENT,
            data=data,
            pconfig=pconf,
            creation_date=report.creation_date,
        )


def plot(
    data_by_sample: Mapping[str, Mapping[str, Mapping[str, float]]],
    pconfig: Union[Dict[str, Any], SeqContentConfig, None] = None,
) -> Union["SeqContentPlot", str, None]:
    """
    Plot a per-base sequence content heatmap.
    :param data_by_sample: sample -> bin label -> {"a": float, "c": float, "g": float, "t": float}.
        Bin labels are FastQC's "base" field ("1", "2-3", "10-14"). Values are
        percentages, except very old FastQC data which gives raw counts.
    :param pconfig: optional dict / SeqContentConfig with the usual id/title.
    :return: SeqContentPlot, ready to be inserted into the report.
    """
    inputs: SeqContentNormalizedInputData = SeqContentNormalizedInputData.create(data_by_sample, pconfig)
    inputs = SeqContentNormalizedInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        return None

    return SeqContentPlot.from_inputs(inputs)


class Dataset(BaseDataset):
    samples: List[str]
    rows: List[List[SeqContentBin]]  # parallel to samples
    max_bp: int

    def sample_names(self) -> List[SampleName]:
        return [SampleName(s) for s in self.samples]

    @staticmethod
    def create(
        dataset: BaseDataset,
        samples: List[str],
        rows: List[List[SeqContentBin]],
    ) -> "Dataset":
        max_bp = max((b.end for bins in rows for b in bins), default=0)
        return Dataset(
            **dataset.__dict__,
            samples=samples,
            rows=rows,
            max_bp=max_bp,
        )

    def create_figure(
        self,
        layout: Optional[go.Layout] = None,
        is_log: bool = False,
        is_pct: bool = False,
        **kwargs: Any,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset: expand the (possibly non-uniform)
        bins into a uniform (n_samples, max_bp, 3) uint8 RGB grid. Positions never
        covered by any bin for a sample stay black.
        """
        n_samples = len(self.samples)
        n_cols = max(self.max_bp, 1)
        arr = np.zeros((n_samples, n_cols, 3), dtype=np.uint8)
        for row_idx, bins in enumerate(self.rows):
            for b in bins:
                r, g, bl = bin_rgb(b)
                # Bin covers columns start..end inclusive, 1-based.
                arr[row_idx, b.start - 1 : b.end, 0] = r
                arr[row_idx, b.start - 1 : b.end, 1] = g
                arr[row_idx, b.start - 1 : b.end, 2] = bl

        # Copy, never mutate the shared layout: this figure feeds only the
        # static/kaleido export path (see plot.py Plot.get_figure -> create_figure),
        # the interactive report dumps self.layout as-is via model_dump().
        fig_layout = go.Layout(layout if layout is not None else self.layout)
        scaleratio = _static_scaleratio(fig_layout, x_range=n_cols, y_range=n_samples)
        if scaleratio is not None:
            fig_layout.yaxis.scaleratio = scaleratio

        return go.Figure(
            data=go.Image(z=arr, x0=1, dx=1, colormodel="rgb", **self.trace_params),
            layout=fig_layout,
        )

    def save_data_file(self) -> None:
        data: Dict[str, Dict[str, Dict[str, float]]] = {
            sample: {b.label: {"a": b.a, "c": b.c, "g": b.g, "t": b.t} for b in bins}
            for sample, bins in zip(self.samples, self.rows)
        }
        report.write_data_file(data, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: SeqContentConfig, keep_hidden: bool = True) -> str:
        """
        Format seqcontent as a compact per-sample "pos: T/C/A/G" listing.
        """
        prompt = ""
        for sample, bins in zip(self.samples, self.rows):
            s_name = report.anonymize_sample_name(sample)
            positions = ", ".join(f"{b.label}: T={b.t}/C={b.c}/A={b.a}/G={b.g}" for b in bins)
            prompt += f"{s_name}: {positions}\n"
        return prompt


class SeqContentPlot(Plot[Dataset, SeqContentConfig]):
    datasets: List[Dataset]

    def sample_names(self) -> List[SampleName]:
        names: List[SampleName] = []
        for ds in self.datasets:
            names.extend(SampleName(s) for s in ds.samples)
        return names

    @staticmethod
    def from_inputs(inputs: SeqContentNormalizedInputData) -> Union["SeqContentPlot", str, None]:
        plot_obj = SeqContentPlot.create(
            data=inputs.data,
            pconfig=inputs.pconfig,
            anchor=inputs.anchor,
        )
        inputs.save_to_parquet()
        return plot_obj

    @staticmethod
    def create(
        data: Dict[str, List[SeqContentBin]],
        pconfig: SeqContentConfig,
        anchor: Anchor,
    ) -> "SeqContentPlot":
        samples = list(data.keys())
        rows = [data[s] for s in samples]
        n_samples = len(samples)

        model: Plot[Dataset, SeqContentConfig] = Plot.initialize(
            plot_type=PlotType.SEQCONTENT,
            pconfig=pconfig,
            anchor=anchor,
            n_series_per_dataset=[n_samples],
            n_samples_per_dataset=[n_samples],
            defer_render_if_large=False,
        )

        model.layout.update(
            yaxis=dict(
                tickmode="array",
                tickvals=list(range(n_samples)),
                ticktext=samples,
            ),
        )

        # Number of samples to a reasonable row height in pixels, similar spirit
        # to heatmap.py's n_elements_to_size, capped to feel like the old canvas's
        # 300px default.
        def n_elements_to_size(n: int) -> int:
            if n >= 50:
                return 4
            if n >= 20:
                return 8
            if n >= 10:
                return 15
            if n >= 5:
                return 25
            return 40

        px_per_sample = n_elements_to_size(n_samples)
        model.layout.height = pconfig.height or min(600, max(100, n_samples * px_per_sample) + 100)

        model.datasets = [Dataset.create(model.datasets[0], samples=samples, rows=rows)]

        return SeqContentPlot(**model.__dict__)

    def _plot_ai_header(self) -> str:
        result = super()._plot_ai_header()
        if self.pconfig.xlab:
            result += f"X axis: {self.pconfig.xlab}\n"
        return result
