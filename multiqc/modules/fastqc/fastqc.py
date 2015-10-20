#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

######################################################
#### LOOKING FOR AN EXAMPLE OF HOW MULTIQC WORKS? ####
######################################################
#### Stop! This module is huge and complicated.   ####
#### Have a look at Bowtie or STAR for a simpler  ####
#### example. CONTRIBUTING.md has documentation.  ####
######################################################

from __future__ import print_function
from collections import defaultdict, OrderedDict
import io
import json
import logging
import os
import re
import shutil
import zipfile

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='FastQC', anchor='fastqc', 
        href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/", 
        info="is a quality control tool for high throughput sequence data,"\
        " written by Simon Andrews at the Babraham Institute in Cambridge.")

        self.fastqc_data = defaultdict(lambda: defaultdict(lambda: defaultdict( lambda: defaultdict() ) ) )
        self.fastqc_stats = dict()
        self.fastqc_statuses = defaultdict(lambda: defaultdict())
        
        # Find and parse unzipped FastQC reports
        for f in self.find_log_files('fastqc_data.txt'):
            s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
            self.parse_fastqc_report(f['f'], s_name)
        
        # Find and parse zipped FastQC reportrs
        for f in self.find_log_files('_fastqc.zip', filecontents=False):
            s_name = f['fn'].rstrip('_fastqc.zip')
            try:
                fqc_zip = zipfile.ZipFile(os.path.join(f['root'], f['fn']))
            except BadZipfile:
                continue
            # FastQC zip files should have just one directory inside, containing report
            d_name = fqc_zip.namelist()[0]
            try:
                with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as f:
                    r_data = f.read().decode('utf8')
                    self.parse_fastqc_report(r_data, s_name)
            except KeyError:
                log.warning("Error - can't find fastqc_raw_data.txt in {}".format(f))

        if len(self.fastqc_stats) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.fastqc_stats)))
        
        # Add full paths to self.css and self.js to be included in template
        self.css = [ os.path.realpath(os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_fastqc.css')) ]
        self.js = [ os.path.realpath(os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_fastqc.js')) ]
        
        # Colours to be used for plotting lines
        self.status_colours = {
            'pass': '#5cb85c',
            'warn': '#f0ad4e',
            'fail': '#d9534f',
            'default': '#999'
        }
        
        # Add to the general statistics table
        self.fastqc_stats_table()
        
        # Write the basic stats table data to a file
        self.write_csv_file(self.fastqc_stats, 'multiqc_fastqc.txt')
        
        # Add the statuses to the intro for multiqc_fastqc.js JavaScript to pick up
        self.intro += '<script type="text/javascript">fastqc_passfails = {};</script>'.format(json.dumps(self.fastqc_statuses))
        
        # Start the sections
        self.sections = list()
        
        # Now add each section in order by calling its def
        self.sequence_quality_plot()
        self.per_seq_quality_plot()
        self.sequence_content_plot()
        self.gc_content_plot()
        self.n_content_plot()
        self.seq_length_dist_plot()
        self.seq_dup_levels_plot()
        self.adapter_content_plot()
        self.kmer_content_plot()


    def parse_fastqc_report(self, file_contents, s_name=None, root=None):
        """ Takes contetns from a fastq_data.txt file and parses out required
        statistics and data. Returns a dict with keys 'stats' and 'data'.
        Data is for plotting graphs, stats are for top table. """
        
        section_headings = {
            'sequence_quality': r'>>Per base sequence quality\s+(pass|warn|fail)',
            'per_seq_quality':  r'>>Per sequence quality scores\s+(pass|warn|fail)',
            'sequence_content': r'>>Per base sequence content\s+(pass|warn|fail)',
            'gc_content':       r'>>Per sequence GC content\s+(pass|warn|fail)',
            'n_content':        r'>>Per base N content\s+(pass|warn|fail)',
            'seq_length_dist':  r'>>Sequence Length Distribution\s+(pass|warn|fail)',
            'seq_dup_levels':   r'>>Sequence Duplication Levels\s+(pass|warn|fail)',
            'adapter_content':  r'>>Adapter Content\s+(pass|warn|fail)',
            'kmer_content':     r'>>Kmer Content\s+(pass|warn|fail)',
        }
        stats_regexes = {
            'filename':           r"^Filename\s+(.+)$",
            'total_sequences':    r"Total Sequences\s+(\d+)",
            'sequence_length':    r"Sequence length\s+([\d-]+)",
            'percent_gc':         r"%GC\s+(\d+)",
            'percent_dedup':      r"#Total Deduplicated Percentage\s+([\d\.]+)",
        }
        
        s = defaultdict(lambda: dict())
        s['seq_len_bp'] = 0
        s['seq_len_read_count'] = 0
        self.seq_lengths = set()
        adapter_types = []
        in_module = None
        for l in file_contents.splitlines():
            
            # Search for general stats
            for k, r in stats_regexes.items():
                r_search = re.search(r, l)
                if r_search:
                    try:
                        s[k] = float(r_search.group(1))
                    except ValueError:
                        s[k] = r_search.group(1)                        
            
            # Parse modules
            if in_module is not None:
                if l == ">>END_MODULE":
                    in_module = None
                else:
                    
                    if in_module == 'sequence_quality':
                        quals = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                        if quals:
                            bp = self.avg_bp_from_range(quals.group(1))
                            groups = ['base', 'mean', 'median', 'lower_quart', 'upper_quart', '10_percentile', '90_percentile']
                            for idx, g in enumerate(groups):
                                try:
                                    self.fastqc_data['sequence_quality'][g][s_name][bp] = float(quals.group( idx + 1 ))
                                except:
                                    self.fastqc_data['sequence_quality'][g][s_name][bp] = quals.group( idx + 1 )
                                    
                    
                    if in_module == 'per_seq_quality' or in_module == 'n_content':
                        sections  = l.split()
                        try:
                            self.fastqc_data[in_module][s_name][float(sections[0])] = float(sections[1])
                        except ValueError:
                            pass # First line - headers
                    
                    if in_module == 'sequence_content':
                        l.replace('NaN','0')
                        seq_matches = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                        if seq_matches:
                            bp = self.avg_bp_from_range(seq_matches.group(1))
                            groups = ['base', 'G', 'A', 'T', 'C']
                            for idx, g in enumerate(groups):
                                if idx == 0:
                                    self.fastqc_data['sequence_content'][s_name][bp][g] = seq_matches.group( idx + 1 )
                                else:
                                    self.fastqc_data['sequence_content'][s_name][bp][g] = float(seq_matches.group( idx + 1 ))
                    
                    if in_module == 'gc_content':
                        gc_matches = re.search("([\d]+)\s+([\d\.E]+)", l)
                        if gc_matches:
                            self.fastqc_data['gc_content'][s_name][int(gc_matches.group(1))] = float(gc_matches.group(2))
                    
                    if in_module == 'seq_length_dist':
                        len_matches = re.search("([\d-]+)\s+([\d\.E]+)", l)
                        if len_matches:
                            bp = self.avg_bp_from_range(len_matches.group(1))
                            self.fastqc_data['seq_length_dist'][s_name][bp] = float(len_matches.group(2))
                            self.seq_lengths.add(bp)
                            s['seq_len_bp'] += float(len_matches.group(2)) * bp
                            s['seq_len_read_count'] += float(len_matches.group(2))
                    
                    if in_module == 'seq_dup_levels':
                        if l[:1] == '#':
                            # Start of module - replace default dict with an OrderedDict
                            self.fastqc_data['seq_dup_levels'][s_name] = OrderedDict()
                            continue # Skip header line
                        sections  = l.split()
                        self.fastqc_data['seq_dup_levels_dedup'][s_name][sections[0]] = float(sections[1])
                        self.fastqc_data['seq_dup_levels'][s_name][sections[0]] = float(sections[2])
        
                    if in_module == 'adapter_content':
                        if l[:1] == '#':
                            adapter_types = l[1:].split("\t")[1:]
                        else:
                            cols = l.split("\t")
                            pos = int(cols[0].split('-', 1)[0])
                            for idx, val in enumerate(cols[1:]):
                                a = adapter_types[idx]
                                k = "{} - {}".format(s_name, a)
                                self.fastqc_data['adapter_content'][k][pos] = float(val)
                    
                    if in_module == 'kmer_content':
                        # fastqc_data.txt doesn't have the breakdown of this
                        # per-base, which is what we need to replicate the
                        # graph in the FastQC report. Do a bar graph instead.
                        if l[:1] == '#':
                            self.fastqc_data['kmer_content_total'][s_name]['Cumulative Max Obs/Exp Scores'] = 0
                            continue
                        sections  = l.split()
                        self.fastqc_data['kmer_content'][s_name][sections[0]] = float(sections[3])
                        self.fastqc_data['kmer_content_total'][s_name]['Cumulative Max Obs/Exp Scores'] += float(sections[3])
                                            
            else:
                # See if this is the start of a new section
                for k, r in section_headings.items():
                    r_search = re.search(r, l)
                    if r_search:
                        in_module = k
                        # Add to the global statuses dict
                        self.fastqc_statuses[k][s_name] = r_search.group(1)
        
        # Work out the average sequence length
        if s['seq_len_read_count'] > 0:
            s['avg_sequence_length'] = s['seq_len_bp'] / s['seq_len_read_count']
        
        # Get percent duplicates (percent unique given)
        if 'percent_dedup' in s:
            s['percent_duplicates'] = 100 - s['percent_dedup']
        
        # Make the sample name from the input filename if we found it
        if 'filename' in s:
            s_name = self.clean_s_name(s['filename'], root)
        
        # Throw a warning if we already have this sample
        # Unzipped reports means that this can be quite frequent
        if s_name in self.fastqc_stats:
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        
        # Add parsed stats to dicts
        self.fastqc_stats[s_name] = s

    def fastqc_stats_table(self):
        """ Add some single-number stats to the basic statistics
        table at the top of the report """
        
        headers = OrderedDict()
        headers['percent_duplicates'] = {
            'title': '% Dups',
            'description': '% Duplicate Reads',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn-rev',
            'format': '{:.1f}%'
        }
        headers['percent_gc'] = {
            'title': '% GC',
            'description': 'Average % GC Content',
            'max': 80,
            'min': 20,
            'scale': 'PRGn',
            'format': '{:.0f}%'
        }
        headers['avg_sequence_length'] = {
            'title': 'Length',
            'description': 'Average Sequence Length (bp)',
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.0f}'
        }
        headers['total_sequences'] = {
            'title': 'M Seqs',
            'description': 'Total Sequences (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.fastqc_stats, headers)


    def sequence_quality_plot (self):
        """ Create the HTML for the phred quality score plot """
        if 'sequence_quality' not in self.fastqc_data or len(self.fastqc_data['sequence_quality']) == 0:
            log.debug('sequence_quality not found in FastQC reports')
            return None    
        
        pconfig = {
            'id': 'fastqc_sequence_quality_plot',
            'title': 'Mean Quality Scores',
            'ylab': 'Phred Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            'colors': self.get_status_cols('sequence_quality'),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.sections.append({
            'name': 'Sequence Quality Histograms',
            'anchor': 'fastqc_sequence_quality',
            'content': '<p>The mean quality value across each base position in the read. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html" target="_bkank">FastQC help</a>.</p>' + 
                        self.plot_xy_data(self.fastqc_data['sequence_quality']['mean'], pconfig)
        })


    def per_seq_quality_plot (self):
        """ Create the HTML for the per sequence quality score plot """
        if 'per_seq_quality' not in self.fastqc_data or len(self.fastqc_data['per_seq_quality']) == 0:
            log.debug('per_seq_quality not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_per_seq_quality_plot',
            'title': 'Per Sequence Quality Scores',
            'ylab': 'Count',
            'xlab': 'Mean Sequence Quality (Phred Score)',
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': self.get_status_cols('per_seq_quality'),
            'tt_label': '<b>Phred {point.x}</b>: {point.y} reads',
            'xPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.sections.append({
            'name': 'Per Sequence Quality Scores',
            'anchor': 'fastqc_per_seq_quality',
            'content': '<p>The number of reads with average quality scores. Shows if a subset of reads has poor quality. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html" target="_bkank">FastQC help</a>.</p>' +
                        self.plot_xy_data(self.fastqc_data['per_seq_quality'], pconfig)
        })

    
    def sequence_content_plot (self):
        """ Create the epic HTML for the FastQC sequence content heatmap """
        if 'sequence_content' not in self.fastqc_data or len(self.fastqc_data['sequence_content']) == 0:
            log.debug('sequence_content not found in FastQC reports')
            return None
        
        html =  '<p>The proportion of each base position for which each of the four normal DNA bases has been called. \
                    See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html" target="_bkank">FastQC help</a>.</p>'
        
        # Order the data by the sample names
        data = OrderedDict(sorted(self.fastqc_data['sequence_content'].items()))
        
        html += '<div id="fastqc_sequence_content_plot"> \n\
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
        
        self.sections.append({
            'name': 'Per Base Sequence Content',
            'anchor': 'fastqc_sequence_content',
            'content': html
        })
    
    
    def gc_content_plot (self):
        """ Create the HTML for the FastQC GC content plot """        
        if 'gc_content' not in self.fastqc_data or len(self.fastqc_data['gc_content']) == 0:
            log.debug('gc_content not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_gc_content_plot',
            'title': 'Per Sequence GC Content',
            'ylab': 'Count',
            'xlab': '%GC',
            'ymin': 0,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            'colors': self.get_status_cols('gc_content'),
        }
        self.sections.append({
            'name': 'Per Sequence GC Content',
            'anchor': 'fastqc_gc_content',
            'content': '<p>The average GC content of reads. Normal random library typically have a roughly normal distribution of GC content. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html" target="_bkank">FastQC help</a>.</p>' +
                        self.plot_xy_data(self.fastqc_data['gc_content'], pconfig)
        })
    
    
    def n_content_plot (self):
        """ Create the HTML for the per base N content plot """
        if 'n_content' not in self.fastqc_data or len(self.fastqc_data['n_content']) == 0:
            log.debug('n_content not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_n_content_plot',
            'title': 'Per Base N Content',
            'ylab': 'Percentage N-Count',
            'xlab': 'Position in Read (bp)',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': self.get_status_cols('n_content'),
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ]
        }
    
        self.sections.append({
            'name': 'Per Base N Content',
            'anchor': 'fastqc_n_content',
            'content': '<p>The percentage of base calls at each position for which an N was called. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html" target="_bkank">FastQC help</a>.</p>' +
                        self.plot_xy_data(self.fastqc_data['n_content'], pconfig)
        })
    
    
    def seq_length_dist_plot (self):
        """ Create the HTML for the Sequence Length Distribution plot """
        if 'seq_length_dist' not in self.fastqc_data or len(self.fastqc_data['seq_length_dist']) == 0:
            log.debug('seq_length_dist not found in FastQC reports')
            return None
        
        if len(self.seq_lengths) < 2:
            html = '<p>All samples have sequences of exactly {} bp in length.</p>'.format(list(self.seq_lengths)[0])
        else:
            pconfig = {
                'id': 'fastqc_seq_length_dist_plot',
                'title': 'Sequence Length Distribution',
                'ylab': 'Read Count',
                'xlab': 'Sequence Length (bp)',
                'ymin': 0,
                'yMinTickInterval': 0.1,
                'xDecimals': False,
                'colors': self.get_status_cols('seq_length_dist'),
                'tt_label': '<b>{point.x} bp</b>: {point.y}',
            }
            html =  '<p>The distribution of fragment sizes (read lengths) found. \
                See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html" target="_bkank">FastQC help</a>.</p>'
            html += self.plot_xy_data(self.fastqc_data['seq_length_dist'], pconfig)
    
        self.sections.append({
            'name': 'Sequence Length Distribution',
            'anchor': 'fastqc_seq_length_dist',
            'content': html
        })
    
    
    def seq_dup_levels_plot (self):
        """ Create the HTML for the Sequence Duplication Levels plot """
        if 'seq_dup_levels' not in self.fastqc_data or len(self.fastqc_data['seq_dup_levels']) == 0:
            log.debug('seq_dup_levels not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_seq_dup_levels_plot',
            'title': 'Sequence Duplication Levels',
            'categories': True,
            'ylab': '% of Library',
            'xlab': 'Sequence Duplication Level',
            'ymax': 100,
            'ymin': 0,
            'yMinTickInterval': 0.1,
            'colors': self.get_status_cols('seq_dup_levels'),
            'tt_label': '<b>{point.x}</b>: {point.y:.1f}%',
        }
    
        self.sections.append({
            'name': 'Sequence Duplication Levels',
            'anchor': 'fastqc_seq_dup_levels',
            'content': '<p>The relative level of duplication found for every sequence. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html" target="_bkank">FastQC help</a>.</p>' +
                        self.plot_xy_data(self.fastqc_data['seq_dup_levels'], pconfig)
        })
        
        
    
    def adapter_content_plot (self):
        """ Create the HTML for the FastQC adapter plot """
        if 'adapter_content' not in self.fastqc_data or len(self.fastqc_data['adapter_content']) == 0:
            log.debug('adapter_content not found in FastQC reports')
            return None
        
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
        # Note - colours are messy as we've added adapter names here. Not
        # possible to break down pass / warn / fail for each adapter, which
        # could lead to misleading labelling (fails on adapter types with
        # little or no contamination)
        
        self.sections.append({
            'name': 'Adapter Content',
            'anchor': 'fastqc_adapter_content',
            'content': '<p>The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html" target="_bkank">FastQC help</a>.</p>' +
                        self.plot_xy_data(self.fastqc_data['adapter_content'], pconfig)
        })
    
    
    def kmer_content_plot (self):
        """ Do the best we can for the FastQC Kmer plot """
        if 'kmer_content' not in self.fastqc_data or len(self.fastqc_data['kmer_content']) == 0:
            log.debug('kmer_content not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_kmer_content_plot',
            'title': 'Kmer Content',
        }
        
        self.sections.append({
            'name': 'Kmer Content',
            'anchor': 'fastqc_kmer_content',
            'content': '<p>The cumulative Observed / Expected scores of over-represented sequences in your library. ' +
                        'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html" target="_bkank">FastQC help</a>.</p>' +
                        self.plot_bargraph(self.fastqc_data['kmer_content_total'], config=pconfig)
        })

    def avg_bp_from_range(self, bp):
        """ Helper function - FastQC often gives base pair ranges (eg. 10-15)
        which are not helpful when plotting. This returns the average from such
        ranges as an int, which is helpful. If not a range, just returns the int """
        
        if '-' in bp:
            maxlen = float(bp.split("-",1)[1])
            minlen = float(bp.split("-",1)[0])
            bp = ((maxlen - minlen)/2) + minlen
        return(int(bp))
    
    def get_status_cols(self, key):
        """ Helper function - returns a list of colours according to the FastQC
        status of this module for each sample. """
        colours = dict()
        for s_name, status in self.fastqc_statuses[key].items():
            colours[s_name] = self.status_colours.get(status, self.status_colours['default'])
        return colours
    