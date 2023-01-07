import copy
import logging
from collections import OrderedDict, defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import boxplot, linegraph

from .util import average_from_range, sortPosQualTableKeys

log = logging.getLogger(__name__)

N_QV = 2


class DragenBaseMetrics(BaseMultiqcModule):
    """
    Rendering ALL THE THINGS!
    """

    def add_base_metrics(self):
        num_samples = len(self.dragen_fastqc_data.keys())
        # Add each section in order
        dragen_fastqc_config = getattr(config, "dragen_fastqc", {})
        quality_range_boxplots_max_samples = dragen_fastqc_config.get("quality_range_boxplots_max_samples", 2)
        force_boxplots = dragen_fastqc_config.get("force_quality_range_boxplots", False)
        if num_samples <= quality_range_boxplots_max_samples or force_boxplots:
            self.positional_quality_range_plot()
        else:
            log.debug(
                f"Skipping DRAGEN-FastQC quality range boxplots as > {quality_range_boxplots_max_samples} samples ({num_samples}). See docs for how to override this behaviour using a config."
            )
        self.positional_mean_quality_plot()

        return self.dragen_fastqc_data.keys()

    def positional_quality_range_plot(self):
        """Box plot image showing range of quality scores at each position"""

        data = OrderedDict()
        data_labels = []
        GROUP = "POSITIONAL QUALITY"
        for s_name in sorted(self.dragen_fastqc_data):
            for mate in sorted(self.dragen_fastqc_data[s_name]):
                r_name = "{}_{}".format(s_name, mate)
                data[r_name] = defaultdict(float)
                data_labels.append({"name": r_name})

                sorted_keys = sortPosQualTableKeys(self.dragen_fastqc_data[s_name][mate][GROUP])
                for key in sorted_keys:
                    value = int(self.dragen_fastqc_data[s_name][mate][GROUP][key])
                    parts = key.split()
                    pos = average_from_range(parts[1])
                    quantile = int(parts[2][:-1])
                    qv = int(value)

                    try:
                        data[r_name][pos][quantile] = qv
                    except:
                        data[r_name][pos] = OrderedDict()
                        data[r_name][pos][quantile] = qv

        pconfig = {
            "id": "dragen_fastqc_per_base_sequence_quality_range_plot",
            "title": "DRAGEN-QC: Per-Position Quality Range",
            "ylab": "Phred Quality Score",
            "xlab": "Position (bp)",
            "ymin": 0,
            "ymax": 43,
            "tt_label": "<b>Base {point.x}</b>: {point.y:.2f}",
            # 'colors': self.get_status_cols('per_base_sequence_quality'),
            "yPlotBands": [
                {"from": 28, "to": 100, "color": "#c3e6c3"},
                {"from": 20, "to": 28, "color": "#e6dcc3"},
                {"from": 0, "to": 20, "color": "#e6c3c3"},
            ],
            "data_labels": data_labels,
        }

        self.add_section(
            name="Per-Position Quality Score Ranges",
            anchor="dragen_fastqc_pos_qual_ranges",
            description="The range of quality value across each base position in each sample or read",
            plot=boxplot.plot(data, pconfig),
        )

    def positional_mean_quality_plot(self):
        """Create the HTML for the positional mean-quality score plot"""

        base_dict = {"A": {}, "C": {}, "G": {}, "T": {}, "N": {}}

        AVG_GROUP = "POSITIONAL BASE MEAN QUALITY"
        COUNT_GROUP = "POSITIONAL BASE CONTENT"

        data = dict()
        for s_name in sorted(self.dragen_fastqc_data):
            for mate in sorted(self.dragen_fastqc_data[s_name]):
                # Parse our per-base, per-position average qualities into a dictionary
                avgs = copy.deepcopy(base_dict)
                for key, value in self.dragen_fastqc_data[s_name][mate][AVG_GROUP].items():
                    if value == "NA":
                        continue
                    parts = key.split()
                    pos = average_from_range(parts[1])
                    base = parts[2].upper()
                    avgs[base][pos] = float(value)

                # Parse matching per-base and total counts by position
                counts = copy.deepcopy(base_dict)
                totals = defaultdict(int)
                for key, value in self.dragen_fastqc_data[s_name][mate][COUNT_GROUP].items():
                    parts = key.split()
                    pos = average_from_range(parts[1])
                    base = parts[2].upper()
                    counts[base][pos] = int(value)
                    totals[pos] += int(value)

                # Use the the count and averages to recompute total QVs
                qv_sums = defaultdict(int)
                for base, pos_data in counts.items():
                    for pos, count in pos_data.items():
                        if count == 0:
                            continue
                        elif base == "N":
                            qv_sums[pos] += count * N_QV
                        else:
                            qv_sums[pos] += int(round(count * avgs[base][pos]))

                # Compute the positional, base-agnostic mean QV
                r_name = "{}_{}".format(s_name, mate)
                data[r_name] = dict()
                for pos, qv_sum in qv_sums.items():
                    total = totals[pos]
                    if total > 0:
                        data[r_name][int(pos)] = qv_sum / total

        pconfig = {
            "id": "dragen_fastqc_per_base_sequence_quality_plot",
            "title": "DRAGEN-QC: Per-Position Quality Scores",
            "ylab": "Phred Quality Score",
            "xlab": "Position (bp)",
            "ymin": 0,
            "xDecimals": False,
            "tt_label": "<b>Base {point.x}</b>: {point.y:.2f}",
            # 'colors': self.get_status_cols('per_base_sequence_quality'),
            "yPlotBands": [
                {"from": 28, "to": 100, "color": "#c3e6c3"},
                {"from": 20, "to": 28, "color": "#e6dcc3"},
                {"from": 0, "to": 20, "color": "#e6c3c3"},
            ],
        }

        self.add_section(
            name="Per-Position Mean Quality Scores",
            anchor="dragen_fastqc_per_base_sequence_quality",
            description="The mean quality value across each base position in the read.",
            helptext="""
            To enable multiple samples to be plotted on the same graph, only the mean quality
            scores are plotted (unlike the box plots seen in FastQC reports).

            Taken from the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html):

            _The y-axis on the graph shows the quality scores. The higher the score, the better
            the base call. The background of the graph divides the y axis into very good quality
            calls (green), calls of reasonable quality (orange), and calls of poor quality (red).
            The quality of calls on most platforms will degrade as the run progresses, so it is
            common to see base calls falling into the orange area towards the end of a read._
            """,
            plot=linegraph.plot(data, pconfig),
        )
