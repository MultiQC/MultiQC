import json

# Initialise the logger
import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph
from multiqc import report

from .util import average_from_range, average_pos_from_metric

log = logging.getLogger(__name__)

N_QV = 2
ADAPTER_SEQS = ["AGATCGGAAGAG", "ATGGAATTCTCG", "CTGTCTCTTATA"]


class DragenContentMetrics(BaseMultiqcModule):
    """
    Rendering ALL THE THINGS!
    """

    def add_content_metrics(self):
        # Now add each section in order
        self.n_content_plot()
        self.sequence_content_plot()
        self.adapter_content_plot()

        return self.dragen_fastqc_data.keys()

    def n_content_plot(self):
        """Create the HTML for the per base N content plot"""
        data = dict()
        GROUP = "POSITIONAL BASE CONTENT"
        for s_name in sorted(self.dragen_fastqc_data):
            for mate in sorted(self.dragen_fastqc_data[s_name]):
                totals = defaultdict(int)
                non_n = defaultdict(int)

                # Count total bases
                total_group_data = self.dragen_fastqc_data[s_name][mate][GROUP]
                for metric, value in total_group_data.items():
                    avg_pos = average_pos_from_metric(metric)
                    totals[avg_pos] += value
                    base = metric.split()[2]
                    if base != "N":
                        non_n[avg_pos] += value

                # Convert Total and Non-N counts into N%
                r_name = f"{s_name}_{mate}"
                data[r_name] = dict()
                for pos, count in totals.items():
                    if count == 0:
                        continue
                    non_n_count = non_n[pos]
                    n_count = count - non_n_count
                    n_frac = 100.0 * n_count / float(count)
                    data[r_name][pos] = n_frac

        if len(data) == 0:
            log.debug("per_base_n_content not found in DRAGEN FastQC reports")
            return None

        pconfig = {
            "id": "dragenqc_per_base_n_content_plot",
            "title": "DRAGEN-QC: Per-Position N Content",
            "ylab": "Percentage N-Count",
            "xlab": "Position in Read (bp)",
            "y_clipmax": 100,
            "y_minrange": 5,
            "ymin": 0,
            "xmin": 0,
            "tt_label": "<b>Base {point.x}</b>: {point.y:.2f}%",
            "y_bands": [
                {"from": 20, "to": 100, "color": "#e6c3c3"},
                {"from": 5, "to": 20, "color": "#e6dcc3"},
                {"from": 0, "to": 5, "color": "#c3e6c3"},
            ],
        }

        self.add_section(
            name="Per-Position N Content",
            anchor="dragenqc_per_base_n_content",
            description="The percentage of base calls at each position for which an `N` was called.",
            helptext="""
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html):
            _If a sequencer is unable to make a base call with sufficient confidence then it will
            normally substitute an `N` rather than a conventional base call. This graph shows the
            percentage of base calls at each position for which an `N` was called._
            _It's not unusual to see a very low proportion of Ns appearing in a sequence, especially
            nearer the end of a sequence. However, if this proportion rises above a few percent
            it suggests that the analysis pipeline was unable to interpret the data well enough to
            make valid base calls._
            """,
            plot=linegraph.plot(data, pconfig),
        )

    def sequence_content_plot(self):
        """Create the epic HTML for the FastQC sequence content heatmap"""

        # Prep the data
        data = dict()
        GROUP = "POSITIONAL BASE CONTENT"
        for s_name in sorted(self.dragen_fastqc_data):
            for mate in sorted(self.dragen_fastqc_data[s_name]):
                r_name = f"{s_name}_{mate}"
                data[r_name] = dict()
                group_data = self.dragen_fastqc_data[s_name][mate][GROUP]

                totals = defaultdict(int)
                for metric, value in group_data.items():
                    parts = metric.split()
                    # avg_pos = average_from_range(parts[1])
                    avg_pos = parts[1]
                    base = parts[-2].lower()

                    if avg_pos not in data[r_name]:
                        data[r_name][avg_pos] = dict()

                    # Store the current count and add it to the total
                    data[r_name][avg_pos][base] = value
                    totals[avg_pos] += value

                # Use the accumulated totals to normalize each bin to a percentage
                for pos, total in totals.items():
                    if total == 0:
                        del data[r_name][pos]
                        continue
                    for base in "acgt":
                        try:
                            data[r_name][pos][base] = (float(data[r_name][pos][base]) / float(total)) * 100.0
                        except Exception:
                            pass
                    data[r_name][pos]["base"] = pos

        if len(data) == 0:
            log.debug("sequence_content not found in FastQC reports")
            return None

        html = """<div id="dragen_fastqc_per_base_sequence_content_plot_div">
            <div class="alert alert-info">
               <span class="glyphicon glyphicon-hand-up"></span>
               Click a sample row to see a line plot for that dataset.
            </div>
            <h5><span class="s_name text-primary"><span class="glyphicon glyphicon-info-sign"></span> Rollover for sample name</span></h5>
            <div class="fastqc_seq_heatmap_key">
                Position: <span id="fastqc_seq_heatmap_key_pos">-</span>
                <div><span id="fastqc_seq_heatmap_key_t"> %T: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_c"> %C: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_a"> %A: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_g"> %G: <span>-</span></span></div>
            </div>
            <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">
                <div id="{id}" class="dragen_fastqc_per_base_sequence_content_plot hc-plot has-custom-export">
                    <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas>
                </div>
            </div>
            <div class="clearfix"></div>
        </div>
        <script type="application/json" class="fastqc_seq_content">{d}</script>
        """.format(
            # Generate unique plot ID, needed in mqc_export_selectplots
            id=report.save_htmlid("dragen_fastqc_per_base_sequence_content_plot"),
            d=json.dumps([self.anchor.replace("-", "_"), data]),
        )

        self.add_section(
            name="Per-Position Sequence Content",
            anchor="dragen_fastqc_per_base_sequence_content",
            description="The proportion of each base position for which each of the four normal DNA bases has been called.",
            helptext="""
            To enable multiple samples to be shown in a single plot, the base composition data
            is shown as a heatmap. The colours represent the balance between the four bases:
            an even distribution should give an even muddy brown colour. Hover over the plot
            to see the percentage of the four bases under the cursor.
            **To see the data as a line plot, as in the original FastQC graph, click on a sample track.**
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html):
            _Per Base Sequence Content plots out the proportion of each base position in a
            file for which each of the four normal DNA bases has been called._
            _In a random library you would expect that there would be little to no difference
            between the different bases of a sequence run, so the lines in this plot should
            run parallel with each other. The relative amount of each base should reflect
            the overall amount of these bases in your genome, but in any case they should
            not be hugely imbalanced from each other._
            _It's worth noting that some types of library will always produce biased sequence
            composition, normally at the start of the read. Libraries produced by priming
            using random hexamers (including nearly all RNA-Seq libraries) and those which
            were fragmented using transposases inherit an intrinsic bias in the positions
            at which reads start. This bias does not concern an absolute sequence, but instead
            provides enrichement of a number of different K-mers at the 5' end of the reads.
            Whilst this is a true technical bias, it isn't something which can be corrected
            by trimming and in most cases doesn't seem to adversely affect the downstream
            analysis._
            """,
            content=html,
        )

    def adapter_content_plot(self):
        """Create the epic HTML for the FastQC adapter content plot"""

        # Prep the data
        data = dict()
        COUNT_GROUP = "POSITIONAL BASE CONTENT"
        ADP_GROUP = "SEQUENCE POSITIONS"
        for s_name in sorted(self.dragen_fastqc_data):
            for mate in sorted(self.dragen_fastqc_data[s_name]):
                totals = defaultdict(int)
                for key, value in self.dragen_fastqc_data[s_name][mate][COUNT_GROUP].items():
                    parts = key.split()
                    pos = average_from_range(parts[1])
                    totals[pos] += int(value)

                adps = defaultdict(int)
                for key, value in self.dragen_fastqc_data[s_name][mate][ADP_GROUP].items():
                    parts = key.split()
                    seq = parts[0].split("'")[1]
                    if seq not in ADAPTER_SEQS:
                        continue
                    # avoid issues with metrics like "'AGATCGGAAGAG' Total Sequence Starts" where
                    # code attempts to parse 'Tot' as an integer
                    if not parts[1][:-2].isnumeric():
                        continue
                    pos = average_from_range(parts[1][:-2])
                    adps[pos] += int(value)

                r_name = f"{s_name}_{mate}"
                data[r_name] = dict()
                cumsum = 0
                for pos, adp_count in sorted(adps.items()):
                    total = totals[pos]
                    cumsum += adp_count
                    if total > 0 and cumsum > 0:
                        data[r_name][pos] = 100.0 * cumsum / total

        pconfig = {
            "id": "dragen_fastqc_adapter_content_plot",
            "title": "FastQC: Adapter Content",
            "ylab": "% of Sequences",
            "xlab": "Position (bp)",
            "y_clipmax": 100,
            "y_minrange": 5,
            "ymin": 0,
            "tt_label": "<b>Base {point.x}</b>: {point.y:.2f}%",
            "hide_empty": True,
            "y_bands": [
                {"from": 20, "to": 100, "color": "#e6c3c3"},
                {"from": 5, "to": 20, "color": "#e6dcc3"},
                {"from": 0, "to": 5, "color": "#c3e6c3"},
            ],
        }

        self.add_section(
            name="Adapter Content",
            anchor="dragen_fastqc_adapter_content",
            description="""The cumulative percentage count of the proportion of your
            library which has seen each of the adapter sequences at each position.""",
            helptext="""
            Note that only samples with ≥ 0.1% adapter contamination are shown.

            There may be several lines per sample, as one is shown for each adapter
            detected in the file.

            From the [FastQC Help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html):

            _The plot shows a cumulative percentage count of the proportion
            of your library which has seen each of the adapter sequences at each position.
            Once a sequence has been seen in a read it is counted as being present
            right through to the end of the read so the percentages you see will only
            increase as the read length goes on._
            """,
            plot=linegraph.plot(data, pconfig),
        )
