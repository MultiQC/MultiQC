#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

############################################################
######  LOOKING FOR AN EXAMPLE OF HOW MULTIQC WORKS?  ######
############################################################
#### Stop! This is one of the most complicated modules. ####
#### Have a look at Kallisto for a simpler example.     ####
############################################################

from __future__ import print_function
from collections import OrderedDict
import io
import json
import logging
import os
import re
import zipfile

from multiqc import config
from multiqc.plots import linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='FastQC', anchor='fastqc',
        href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
        info="is a quality control tool for high throughput sequence data,"\
        " written by Simon Andrews at the Babraham Institute in Cambridge.")

        self.fastqc_data = dict()

        # Find and parse unzipped FastQC reports
        for f in self.find_log_files('fastqc/data'):
            s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
            self.parse_fastqc_report(f['f'], s_name, f)

        # Find and parse zipped FastQC reports
        for f in self.find_log_files('fastqc/zip', filecontents=False):
            s_name = f['fn']
            if s_name.endswith('_fastqc.zip'):
                s_name = s_name[:-11]
            # Skip if we already have this report - parsing zip files is slow..
            if s_name in self.fastqc_data.keys():
                log.debug("Skipping '{}' as already parsed '{}'".format(f['fn'], s_name))
                continue
            try:
                fqc_zip = zipfile.ZipFile(os.path.join(f['root'], f['fn']))
            except Exception as e:
                log.warn("Couldn't read '{}' - Bad zip file".format(f['fn']))
                log.debug("Bad zip file error:\n{}".format(e))
                continue
            # FastQC zip files should have just one directory inside, containing report
            d_name = fqc_zip.namelist()[0]
            try:
                with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as fh:
                    r_data = fh.read().decode('utf8')
                    self.parse_fastqc_report(r_data, s_name, f)
            except KeyError:
                log.warning("Error - can't find fastqc_raw_data.txt in {}".format(f))

        # Filter to strip out ignored sample names
        self.fastqc_data = self.ignore_samples(self.fastqc_data)

        if len(self.fastqc_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.fastqc_data)))

        # Write the summary stats to a file
        data = dict()
        for s_name in self.fastqc_data:
            data[s_name] = self.fastqc_data[s_name]['basic_statistics']
            data[s_name].update(self.fastqc_data[s_name]['statuses'])
        self.write_data_file(data, 'multiqc_fastqc')

        # Add to self.css and self.js to be included in template
        self.css = { 'assets/css/multiqc_fastqc.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_fastqc.css') }
        self.js = { 'assets/js/multiqc_fastqc.js' : os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_fastqc.js') }

        # Colours to be used for plotting lines
        self.status_colours = { 'pass': '#5cb85c', 'warn': '#f0ad4e', 'fail': '#d9534f', 'default': '#999' }

        # Add to the general statistics table
        self.fastqc_general_stats()

        # Add the statuses to the intro for multiqc_fastqc.js JavaScript to pick up
        statuses = dict()
        for s_name in self.fastqc_data:
            for section, status in self.fastqc_data[s_name]['statuses'].items():
                try:
                    statuses[section][s_name] = status
                except KeyError:
                    statuses[section] = {s_name: status}
        self.intro += '<script type="text/javascript">fastqc_passfails = {};</script>'.format(json.dumps(statuses))

        # Now add each section in order
        self.sequence_quality_plot()
        self.per_seq_quality_plot()
        self.sequence_content_plot()
        self.gc_content_plot()
        self.n_content_plot()
        self.seq_length_dist_plot()
        self.seq_dup_levels_plot()
        self.overrepresented_sequences()
        self.adapter_content_plot()

    def parse_fastqc_report(self, file_contents, s_name=None, f=None):
        """ Takes contents from a fastq_data.txt file and parses out required
        statistics and data. Returns a dict with keys 'stats' and 'data'.
        Data is for plotting graphs, stats are for top table. """

        # Make the sample name from the input filename if we find it
        fn_search = re.search(r"Filename\s+(.+)", file_contents)
        if fn_search:
            s_name = self.clean_s_name(fn_search.group(1) , f['root'])

        if s_name in self.fastqc_data.keys():
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        self.add_data_source(f, s_name)
        self.fastqc_data[s_name] = { 'statuses': dict() }

        # Parse the report
        section = None
        s_headers = None
        self.dup_keys = []
        for l in file_contents.splitlines():
            if l == '>>END_MODULE':
                section = None
                s_headers = None
            elif l.startswith('>>'):
                (section, status) = l[2:].split("\t", 1)
                section = section.lower().replace(' ', '_')
                self.fastqc_data[s_name]['statuses'][section] = status
            elif section is not None:
                if l.startswith('#'):
                    s_headers = l[1:].split("\t")
                    # Special case: Total Deduplicated Percentage header line
                    if s_headers[0] == 'Total Deduplicated Percentage':
                        self.fastqc_data[s_name]['basic_statistics'].append({
                            'measure': 'total_deduplicated_percentage',
                            'value': float(s_headers[1])
                        })
                    else:
                        # Special case: Rename dedup header in old versions of FastQC (v10)
                        if s_headers[1] == 'Relative count':
                            s_headers[1] = 'Percentage of total'
                        s_headers = [s.lower().replace(' ', '_') for s in s_headers]
                        self.fastqc_data[s_name][section] = list()

                elif s_headers is not None:
                    s = l.split("\t")
                    row = dict()
                    for (i, v) in enumerate(s):
                        v.replace('NaN','0')
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                        row[s_headers[i]] = v
                    self.fastqc_data[s_name][section].append(row)
                    # Special case - need to remember order of duplication keys
                    if section == 'sequence_duplication_levels':
                        try:
                            self.dup_keys.append(float(s[0]))
                        except ValueError:
                            self.dup_keys.append(s[0])

        # Tidy up the Basic Stats
        self.fastqc_data[s_name]['basic_statistics'] = {d['measure']: d['value'] for d in self.fastqc_data[s_name]['basic_statistics']}

        # Calculate the average sequence length (Basic Statistics gives a range)
        length_bp = 0
        total_count = 0
        for d in self.fastqc_data[s_name].get('sequence_length_distribution', {}):
            length_bp += d['count'] * self.avg_bp_from_range(d['length'])
            total_count += d['count']
        if total_count > 0:
            self.fastqc_data[s_name]['basic_statistics']['avg_sequence_length'] = length_bp / total_count

    def fastqc_general_stats(self):
        """ Add some single-number stats to the basic statistics
        table at the top of the report """

        # Prep the data
        data = dict()
        for s_name in self.fastqc_data:
            bs = self.fastqc_data[s_name]['basic_statistics']
            data[s_name] = {
                'percent_gc': bs['%GC'],
                'avg_sequence_length': bs['avg_sequence_length'],
                'total_sequences': bs['Total Sequences'],
            }
            try:
                data[s_name]['percent_duplicates'] = 100 - bs['total_deduplicated_percentage']
            except KeyError:
                pass # Older versions of FastQC don't have this
            # Add count of fail statuses
            num_statuses = 0
            num_fails = 0
            for s in self.fastqc_data[s_name]['statuses'].values():
                num_statuses += 1
                if s == 'fail':
                    num_fails += 1
            data[s_name]['percent_fails'] = (float(num_fails)/float(num_statuses))*100.0

        # Are sequence lengths interesting?
        seq_lengths = [x['avg_sequence_length'] for x in data.values()]
        hide_seq_length = False if max(seq_lengths) - min(seq_lengths) > 10 else True

        headers = OrderedDict()
        headers['percent_duplicates'] = {
            'title': '% Dups',
            'description': '% Duplicate Reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['percent_gc'] = {
            'title': '% GC',
            'description': 'Average % GC Content',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Set1',
            'format': '{:,.0f}'
        }
        headers['avg_sequence_length'] = {
            'title': 'Length',
            'description': 'Average Sequence Length (bp)',
            'min': 0,
            'suffix': ' bp',
            'scale': 'RdYlGn',
            'format': '{:,.0f}',
            'hidden': hide_seq_length
        }
        headers['percent_fails'] = {
            'title': '% Failed',
            'description': 'Percentage of modules failed in FastQC report (includes those not plotted here)',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
            'format': '{:,.0f}',
            'hidden': True
        }
        headers['total_sequences'] = {
            'title': '{} Seqs'.format(config.read_count_prefix),
            'description': 'Total Sequences ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(data, headers)


    def sequence_quality_plot (self):
        """ Create the HTML for the phred quality score plot """

        data = dict()
        for s_name in self.fastqc_data:
            try:
                data[s_name] = {self.avg_bp_from_range(d['base']): d['mean'] for d in self.fastqc_data[s_name]['per_base_sequence_quality']}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('sequence_quality not found in FastQC reports')
            return None

        pconfig = {
            'id': 'fastqc_per_base_sequence_quality_plot',
            'title': 'Mean Quality Scores',
            'ylab': 'Phred Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            'colors': self.get_status_cols('per_base_sequence_quality'),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.add_section (
            name = 'Sequence Quality Histograms',
            anchor = 'fastqc_per_base_sequence_quality',
            description = 'The mean quality value across each base position in the read. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html" target="_bkank">FastQC help</a>.',
            plot = linegraph.plot(data, pconfig)
        )


    def per_seq_quality_plot (self):
        """ Create the HTML for the per sequence quality score plot """

        data = dict()
        for s_name in self.fastqc_data:
            try:
                data[s_name] = {d['quality']: d['count'] for d in self.fastqc_data[s_name]['per_sequence_quality_scores']}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('per_seq_quality not found in FastQC reports')
            return None

        pconfig = {
            'id': 'fastqc_per_sequence_quality_scores_plot',
            'title': 'Per Sequence Quality Scores',
            'ylab': 'Count',
            'xlab': 'Mean Sequence Quality (Phred Score)',
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': self.get_status_cols('per_sequence_quality_scores'),
            'tt_label': '<b>Phred {point.x}</b>: {point.y} reads',
            'xPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.add_section (
            name = 'Per Sequence Quality Scores',
            anchor = 'fastqc_per_sequence_quality_scores',
            description = 'The number of reads with average quality scores. Shows if a subset of reads has poor quality. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html" target="_bkank">FastQC help</a>.',
            plot = linegraph.plot(data, pconfig)
        )


    def sequence_content_plot (self):
        """ Create the epic HTML for the FastQC sequence content heatmap """

        html =  '<p>The proportion of each base position for which each of the four normal DNA bases has been called. \
                    See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html" target="_bkank">FastQC help</a>.</p> \
                 <p class="text-primary"><span class="glyphicon glyphicon-info-sign"></span> Click a heatmap row to see a line plot for that dataset.</p>'

        # Prep the data
        data = OrderedDict()
        for s_name in sorted(self.fastqc_data.keys()):
            try:
                data[s_name] = {self.avg_bp_from_range(d['base']): d for d in self.fastqc_data[s_name]['per_base_sequence_content']}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('sequence_content not found in FastQC reports')
            return None

        html += '<div id="fastqc_per_base_sequence_content_plot"> \n\
            <h5><span class="s_name"><em class="text-muted">rollover for sample name</em></span></h5> \n\
            <div class="fastqc_seq_heatmap_key">\n\
                Position: <span id="fastqc_seq_heatmap_key_pos">-</span>\n\
                <div><span id="fastqc_seq_heatmap_key_t"> %T: <span>-</span></span></div>\n\
                <div><span id="fastqc_seq_heatmap_key_c"> %C: <span>-</span></span></div>\n\
                <div><span id="fastqc_seq_heatmap_key_a"> %A: <span>-</span></span></div>\n\
                <div><span id="fastqc_seq_heatmap_key_g"> %G: <span>-</span></span></div>\n\
            </div>\n\
            <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">\n\
                <div id="fastqc_seq" class="hc-plot"> \n\
                    <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas> \n\
                </div> \n\
            </div> \n\
            <div class="clearfix"></div> \n\
        </div> \n\
        <script type="text/javascript"> \n\
            fastqc_seq_content_data = {d};\n\
            $(function () {{ fastqc_seq_content_heatmap(); }}); \n\
        </script>'.format(d=json.dumps(data))

        self.add_section (
            name = 'Per Base Sequence Content',
            anchor = 'fastqc_per_base_sequence_content',
            content = html
        )


    def gc_content_plot (self):
        """ Create the HTML for the FastQC GC content plot """

        data = dict()
        data_norm = dict()
        for s_name in self.fastqc_data:
            try:
                data[s_name] = {d['gc_content']: d['count'] for d in self.fastqc_data[s_name]['per_sequence_gc_content']}
            except KeyError:
                pass
            else:
                data_norm[s_name] = dict()
                total = sum( [ c for c in data[s_name].values() ] )
                for gc, count in data[s_name].items():
                    data_norm[s_name][gc] = (count / total) * 100
        if len(data) == 0:
            log.debug('per_sequence_gc_content not found in FastQC reports')
            return None

        pconfig = {
            'id': 'fastqc_per_sequence_gc_content_plot',
            'title': 'Per Sequence GC Content',
            'ylab': 'Count',
            'xlab': '%GC',
            'ymin': 0,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            'colors': self.get_status_cols('per_sequence_gc_content'),
            'data_labels': [
                {'name': 'Percentages', 'ylab': 'Percentage'},
                {'name': 'Counts', 'ylab': 'Count'}
            ]
        }

        # Try to find and plot a theoretical GC line
        theoretical_gc = None
        theoretical_gc_raw = None
        theoretical_gc_name = None
        for f in self.find_log_files('fastqc/theoretical_gc'):
            if theoretical_gc_raw is not None:
                log.warn("Multiple FastQC Theoretical GC Content files found, now using {}".format(f['fn']))
            theoretical_gc_raw = f['f']
            theoretical_gc_name = f['fn']
        if theoretical_gc_raw is None:
            tgc = getattr(config, 'fastqc_config', {}).get('fastqc_theoretical_gc', None)
            if tgc is not None:
                theoretical_gc_name = os.path.basename(tgc)
                tgc_fn = 'fastqc_theoretical_gc_{}.txt'.format(tgc)
                tgc_path = os.path.join(os.path.dirname(__file__), 'fastqc_theoretical_gc', tgc_fn)
                if not os.path.isfile(tgc_path):
                    tgc_path = tgc
                try:
                    with io.open (tgc_path, "r", encoding='utf-8') as f:
                        theoretical_gc_raw = f.read()
                except IOError:
                    log.warn("Couldn't open FastQC Theoretical GC Content file {}".format(tgc_path))
                    theoretical_gc_raw = None
        if theoretical_gc_raw is not None:
            theoretical_gc = list()
            for l in theoretical_gc_raw.splitlines():
                if '# FastQC theoretical GC content curve:' in l:
                    theoretical_gc_name = l[39:]
                elif not l.startswith('#'):
                    s = l.split()
                    try:
                        theoretical_gc.append([float(s[0]), float(s[1])])
                    except (TypeError, IndexError):
                        pass

        desc = '''The average GC content of reads. Normal random library typically have a roughly normal distribution of GC content.
                See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html" target="_bkank">FastQC help</a>.</p>'''
        if theoretical_gc is not None:
            # Calculate the count version of the theoretical data based on the largest data store
            max_total = max([sum (d.values()) for d in data.values() ])
            esconfig = {
                'name': 'Theoretical GC Content',
                'dashStyle': 'Dash',
                'lineWidth': 2,
                'color': '#000000',
                'marker': { 'enabled': False },
                'enableMouseTracking': False,
                'showInLegend': False,
            }
            pconfig['extra_series'] = [ [dict(esconfig)], [dict(esconfig)] ]
            pconfig['extra_series'][0][0]['data'] = theoretical_gc
            pconfig['extra_series'][1][0]['data'] = [ [ d[0], (d[1]/100.0)*max_total ] for d in theoretical_gc ]
            desc = '</p><p>The dashed black line shows theoretical GC content: {}.'.format(theoretical_gc_name)

        self.add_section (
            name = 'Per Sequence GC Content',
            anchor = 'fastqc_per_sequence_gc_content',
            description = desc,
            plot = linegraph.plot([data_norm, data], pconfig)
        )


    def n_content_plot (self):
        """ Create the HTML for the per base N content plot """

        data = dict()
        for s_name in self.fastqc_data:
            try:
                data[s_name] = {self.avg_bp_from_range(d['base']): d['n-count'] for d in self.fastqc_data[s_name]['per_base_n_content']}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('per_base_n_content not found in FastQC reports')
            return None

        pconfig = {
            'id': 'fastqc_per_base_n_content_plot',
            'title': 'Per Base N Content',
            'ylab': 'Percentage N-Count',
            'xlab': 'Position in Read (bp)',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': self.get_status_cols('per_base_n_content'),
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ]
        }

        self.add_section (
            name = 'Per Base N Content',
            anchor = 'fastqc_per_base_n_content',
            description = 'The percentage of base calls at each position for which an N was called. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html" target="_bkank">FastQC help</a>.',
            plot = linegraph.plot(data, pconfig)
        )


    def seq_length_dist_plot (self):
        """ Create the HTML for the Sequence Length Distribution plot """

        data = dict()
        seq_lengths = set()
        multiple_lenths = False
        for s_name in self.fastqc_data:
            try:
                data[s_name] = {self.avg_bp_from_range(d['length']): d['count'] for d in self.fastqc_data[s_name]['sequence_length_distribution']}
                seq_lengths.update(data[s_name].keys())
                if len(set(data[s_name].keys())) > 1:
                    multiple_lenths = True
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('sequence_length_distribution not found in FastQC reports')
            return None

        if not multiple_lenths:
            lengths = 'bp , '.join([str(l) for l in list(seq_lengths)])
            desc = 'All samples have sequences of a single length ({}bp).'.format(lengths)
            if len(seq_lengths) > 1:
                desc += ' See the <a href="#general_stats">General Statistics Table</a>.'
            self.add_section (
                name = 'Sequence Length Distribution',
                anchor = 'fastqc_sequence_length_distribution',
                description = '<div class="alert alert-info">{}</div>'.format(desc)
            )
        else:
            pconfig = {
                'id': 'fastqc_sequence_length_distribution_plot',
                'title': 'Sequence Length Distribution',
                'ylab': 'Read Count',
                'xlab': 'Sequence Length (bp)',
                'ymin': 0,
                'yMinTickInterval': 0.1,
                'xDecimals': False,
                'colors': self.get_status_cols('sequence_length_distribution'),
                'tt_label': '<b>{point.x} bp</b>: {point.y}',
            }
            desc =  'The distribution of fragment sizes (read lengths) found. \
                See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html" target="_bkank">FastQC help</a>.'
            self.add_section (
                name = 'Sequence Length Distribution',
                anchor = 'fastqc_sequence_length_distribution',
                description = desc,
                plot = linegraph.plot(data, pconfig)
            )


    def seq_dup_levels_plot (self):
        """ Create the HTML for the Sequence Duplication Levels plot """

        data = dict()
        for s_name in self.fastqc_data:
            try:
                d = {d['duplication_level']: d['percentage_of_total'] for d in self.fastqc_data[s_name]['sequence_duplication_levels']}
                data[s_name] = OrderedDict()
                for k in self.dup_keys:
                    try:
                        data[s_name][k] = d[k]
                    except KeyError:
                        pass
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('sequence_length_distribution not found in FastQC reports')
            return None

        pconfig = {
            'id': 'fastqc_sequence_duplication_levels_plot',
            'title': 'Sequence Duplication Levels',
            'categories': True,
            'ylab': '% of Library',
            'xlab': 'Sequence Duplication Level',
            'ymax': 100,
            'ymin': 0,
            'yMinTickInterval': 0.1,
            'colors': self.get_status_cols('sequence_duplication_levels'),
            'tt_label': '<b>{point.x}</b>: {point.y:.1f}%',
        }

        self.add_section (
            name = 'Sequence Duplication Levels',
            anchor = 'fastqc_sequence_duplication_levels',
            description = 'The relative level of duplication found for every sequence. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html" target="_bkank">FastQC help</a>.',
            plot = linegraph.plot(data, pconfig)
        )

    def overrepresented_sequences (self):
        """Sum the percentages of overrepresented sequences and display them in a bar plot"""

        data = dict()
        for s_name in self.fastqc_data:
            data[s_name] = dict()
            try:
                max_pcnt   = max( [ float(d['percentage']) for d in self.fastqc_data[s_name]['overrepresented_sequences']] )
                total_pcnt = sum( [ float(d['percentage']) for d in self.fastqc_data[s_name]['overrepresented_sequences']] )
                data[s_name]['total_overrepresented'] = total_pcnt
                data[s_name]['top_overrepresented'] = max_pcnt
                data[s_name]['remaining_overrepresented'] = total_pcnt - max_pcnt
            except KeyError:
                if self.fastqc_data[s_name]['statuses']['overrepresented_sequences'] == 'pass':
                    data[s_name]['total_overrepresented'] = 0
                    data[s_name]['top_overrepresented'] = 0
                    data[s_name]['remaining_overrepresented'] = 0
                else:
                    log.debug("Couldn't find data for {}, invalid Key".format(s_name))

        cats = OrderedDict()
        cats['top_overrepresented'] = { 'name': 'Top over-represented sequence' }
        cats['remaining_overrepresented'] = { 'name': 'Sum of remaining over-represented sequences' }

        # Config for the plot
        pconfig = {
            'id': 'fastqc_overrepresented_sequencesi_plot',
            'title': 'Overrepresented sequences',
            'ymin': 0,
            'yCeiling': 100,
            'yMinRange': 20,
            'tt_decimals': 2,
            'tt_suffix': '%',
            'tt_percentages': False,
            'ylab_format': '{value}%',
            'cpswitch': False,
            'ylab': 'Percentage of Total Sequences'
        }

        # Check if any samples have more than 1% overrepresented sequences, else don't make plot.
        if max([ x['total_overrepresented'] for x in data.values()]) < 1:
            plot_html = '<div class="alert alert-info">{} samples had less than 1% of reads made up of overrepresented sequences</div>'.format(len(data))
        else:
            plot_html = bargraph.plot(data, cats, pconfig)

        self.add_section (
            name = 'Overrepresented sequences',
            anchor = 'fastqc_overrepresented_sequences',
            description = 'The total amount of overrepresented sequences found in each library. ' +
                    'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html" target="_bkank">FastQC help for further information</a>.',
            plot  = plot_html
        )


    def adapter_content_plot (self):
        """ Create the HTML for the FastQC adapter plot """

        data = dict()
        for s_name in self.fastqc_data:
            try:
                for d in self.fastqc_data[s_name]['adapter_content']:
                    pos = self.avg_bp_from_range(d['position'])
                    for r in self.fastqc_data[s_name]['adapter_content']:
                        pos = self.avg_bp_from_range(r['position'])
                        for a in r.keys():
                            k = "{} - {}".format(s_name, a)
                            if a != 'position':
                                try:
                                    data[k][pos] = r[a]
                                except KeyError:
                                    data[k] = {pos: r[a]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('adapter_content not found in FastQC reports')
            return None

        # Lots of these datasets will be all zeros.
        # Only take datasets with > 0.1% adapter contamination
        data = {k:d for k,d in data.items() if max(data[k].values()) >= 0.1 }

        pconfig = {
            'id': 'fastqc_adapter_content_plot',
            'title': 'Adapter Content',
            'ylab': '% of Sequences',
            'xlab': 'Position',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'hide_empty': True,
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ],
        }

        if len(data) > 0:
            plot_html = linegraph.plot(data, pconfig)
        else:
            plot_html = '<div class="alert alert-info">No samples found with any adapter contamination > 0.1%</div>'

        # Note - colours are messy as we've added adapter names here. Not
        # possible to break down pass / warn / fail for each adapter, which
        # could lead to misleading labelling (fails on adapter types with
        # little or no contamination)

        self.add_section (
            name = 'Adapter Content',
            anchor = 'fastqc_adapter_content',
            description = 'The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html" target="_bkank">FastQC help</a>. ' +
                        'Only samples with &ge; 0.1% adapter contamination are shown.',
            plot = plot_html
        )


    def avg_bp_from_range(self, bp):
        """ Helper function - FastQC often gives base pair ranges (eg. 10-15)
        which are not helpful when plotting. This returns the average from such
        ranges as an int, which is helpful. If not a range, just returns the int """

        try:
            if '-' in bp:
                maxlen = float(bp.split("-",1)[1])
                minlen = float(bp.split("-",1)[0])
                bp = ((maxlen - minlen)/2) + minlen
        except TypeError:
            pass
        return(int(bp))

    def get_status_cols(self, section):
        """ Helper function - returns a list of colours according to the FastQC
        status of this module for each sample. """
        colours = dict()
        for s_name in self.fastqc_data:
            status = self.fastqc_data[s_name]['statuses'].get(section, 'default')
            colours[s_name] = self.status_colours[status]
        return colours
