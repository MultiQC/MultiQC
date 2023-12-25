import logging

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

from .util import average_pos_from_size

log = logging.getLogger(__name__)

N_QV = 2
ADAPTER_SEQS = ["AGATCGGAAGAG", "ATGGAATTCTCG", "CTGTCTCTTATA"]


class DragenReadMetrics(BaseMultiqcModule):
    """
    Rendering ALL THE THINGS!
    """

    def add_read_metrics(self):
        # Add each section in order
        self.per_seq_quality_plot()
        self.seq_length_dist_plot()

        return self.dragen_fastqc_data.keys()

    def per_seq_quality_plot(self):
        """Create the HTML for the per sequence quality score plot"""

        data = dict()
        GROUP = "READ MEAN QUALITY"
        MAX_QV = 64
        max_non_zero = 0
        for s_name in sorted(self.dragen_fastqc_data):
            for mate in sorted(self.dragen_fastqc_data[s_name]):
                r_name = f"{s_name}_{mate}"
                data[r_name] = dict()
                group_data = self.dragen_fastqc_data[s_name][mate][GROUP]
                for qv in range(MAX_QV):
                    metric = f"Q{qv} Reads"
                    count = group_data[metric]
                    if count > 0:
                        max_non_zero = max(qv, max_non_zero)
                    data[r_name][qv] = count

        for r_name in data:
            for qv in range(max_non_zero + 2, MAX_QV):
                del data[r_name][qv]

        if len(data) == 0:
            log.debug("per_seq_quality not found in DRAGEN FastQC reports")
            return None

        pconfig = {
            "id": "dragenqc_per_sequence_quality_scores_plot",
            "title": "DRAGEN-QC: Per-Sequence Quality Scores",
            "ylab": "Count",
            "xlab": "Mean Sequence Quality (Phred Quality Score)",
            "ymin": 0,
            "xmin": 0,
            "xDecimals": False,
            # 'colors': self.get_status_cols('per_sequence_quality_scores'),
            "tt_label": "<b>Phred {point.x}</b>: {point.y} reads",
            "xPlotBands": [
                {"from": 28, "to": 100, "color": "#c3e6c3"},
                {"from": 20, "to": 28, "color": "#e6dcc3"},
                {"from": 0, "to": 20, "color": "#e6c3c3"},
            ],
        }
        self.add_section(
            name="Per-Sequence Quality Scores",
            anchor="dragenqc_per_sequence_quality_scores",
            description="The number of reads with average quality scores. Shows if a subset of reads has poor quality.",
            helptext="""
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html):
            _The per sequence quality score report allows you to see if a subset of your
            sequences have universally low quality values. It is often the case that a
            subset of sequences will have universally poor quality, however these should
            represent only a small percentage of the total sequences._
            """,
            plot=linegraph.plot(data, pconfig),
        )

    def seq_length_dist_plot(self):
        """Create the HTML for the Sequence Length Distribution plot"""

        data = dict()
        seq_lengths = set()
        multiple_lenths = False
        avg_to_range = dict()
        GROUP = "READ LENGTHS"
        for s_name in sorted(self.dragen_fastqc_data):
            for mate in sorted(self.dragen_fastqc_data[s_name]):
                r_name = f"{s_name}_{mate}"
                data[r_name] = dict()

                group_data = self.dragen_fastqc_data[s_name][mate][GROUP]
                for metric, value in group_data.items():
                    if value > 0:
                        avg_pos = average_pos_from_size(metric)
                        data[r_name][avg_pos] = value
                        avg_to_range[avg_pos] = metric.split("bp")[0]

                seq_lengths.update([avg_to_range[k] for k in data[r_name].keys()])

                if len(set(data[r_name].keys())) > 1:
                    multiple_lenths = True

        if len(data) == 0:
            log.debug("sequence_length_distribution not found in FastQC reports")
            return None

        if not multiple_lenths:
            lengths = "bp , ".join([str(line) for line in list(seq_lengths)])
            desc = f"All samples have sequences within a single length bin ({lengths}bp)."
            if len(seq_lengths) > 1:
                desc += ' See the <a href="#general_stats">General Statistics Table</a>.'
            self.add_section(
                name="Sequence Length Distribution",
                anchor="dragenqc_sequence_length_distribution",
                description=f'<div class="alert alert-info">{desc}</div>',
            )
        else:
            pconfig = {
                "id": "dragenqc_sequence_length_distribution_plot",
                "title": "DRAGEN-QC: Sequence Length Distribution",
                "ylab": "Read Count",
                "xlab": "Sequence Length (bp)",
                "ymin": 0,
                "yMinTickInterval": 0.1,
                "xDecimals": False,
                # 'colors': self.get_status_cols('sequence_length_distribution'),
                "tt_label": "<b>{point.x} bp</b>: {point.y}",
            }
            self.add_section(
                name="Sequence Length Distribution",
                anchor="dragenqc_sequence_length_distribution",
                description="""The distribution of fragment sizes (read lengths) found.
                    See the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)""",
                plot=linegraph.plot(data, pconfig),
            )
