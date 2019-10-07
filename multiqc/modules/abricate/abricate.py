#!/usr/bin/env python

""" MultiQC module to parse output from ABRicate """

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import heatmap
import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='ABRicate', anchor='abricate',
        href="https://github.com/tseemann/abricate",
        info="Mass screening of contigs for antimicrobial resistance or virulence genes.")

        # Find all files for mymod
        for myfile in self.find_log_files('abricate'):
            print( "1" )
            print( myfile['f'] )       # File contents
            print( "2" )
            print( myfile['s_name'] )  # Sample name (from cleaned filename)
            print( "3" )
            print( myfile['fn'] )      # Filename
            print( "4" )
            print( myfile['root'] )    # Directory file was in
            print( "5" )

    def plot_gene_presence_heatmap(self):
        """ Return HTML for correlation heatmap """
        if data is not None:
            pconfig = {
                'id': 'gene presence and absence',
                'title': '{} def title'
            }
            self.add_section (
                name = '{} add section title',
                anchor = 'abricate',
                plot = heatmap.plot(data[1], data[0], data[0], pconfig)
            )

    def gene_presence_heatmap (self):
        """ Create the HTML for the abricate summary heatmap """

        # Prep the data
        data = OrderedDict()
        for s_name in sorted(self.fastqc_data.keys()):
            try:
                data[s_name] = {self.avg_bp_from_range(d['base']): d for d in self.fastqc_data[s_name]['per_base_sequence_content']}
            except KeyError:
                pass
            # Old versions of FastQC give counts instead of percentages
            for b in data[s_name]:
                tot = sum([data[s_name][b][base] for base in ['a','c','t','g']])
                if tot == 100.0:
                    break
                else:
                    for base in ['a','c','t','g']:
                        data[s_name][b][base] = (float(data[s_name][b][base])/float(tot)) * 100.0
        if len(data) == 0:
            log.debug('sequence_content not found in FastQC reports')
            return None

        html = '''<div id="fastqc_per_base_sequence_content_plot_div">
            <div class="alert alert-info">
               <span class="glyphicon glyphicon-hand-up"></span>
               Click a sample row to see a line plot for that dataset.
            </div>
            <h5><span class="s_name text-primary"><span class="glyphicon glyphicon-info-sign"></span> Rollover for sample name</span></h5>
            <button id="fastqc_per_base_sequence_content_export_btn"><span class="glyphicon glyphicon-download-alt"></span> Export Plot</button>
            <div class="fastqc_seq_heatmap_key">
                Position: <span id="fastqc_seq_heatmap_key_pos">-</span>
                <div><span id="fastqc_seq_heatmap_key_t"> %T: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_c"> %C: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_a"> %A: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_g"> %G: <span>-</span></span></div>
            </div>
            <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">
                <div id="fastqc_per_base_sequence_content_plot" class="hc-plot has-custom-export">
                    <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas>
                </div>
            </div>
            <div class="clearfix"></div>
        </div>
        <script type="application/json" class="fastqc_seq_content">{d}</script>
        '''.format(d=json.dumps([self.anchor.replace('-', '_'), data]))

        self.add_section (
            name = 'Per Base Sequence Content',
            anchor = 'fastqc_per_base_sequence_content',
            description = 'The proportion of each base position for which each of the four normal DNA bases has been called.',
            helptext = '''
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
            ''',
            content = html
        )

log.info('Hello World!')
