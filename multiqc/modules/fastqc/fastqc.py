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
            fqc_zip = zipfile.ZipFile(os.path.join(f['root'], f['fn']))
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
        self.gc_content_plot()
        self.n_content_plot()
        # self.sequence_content_plot()
        self.adapter_content_plot()


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
            'percent_duplicates': r"#Total Deduplicated Percentage\s+([\d\.]+)",
        }
        
        s = defaultdict(lambda: dict())
        s['seq_len_bp'] = 0
        s['seq_len_read_count'] = 0
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
                            s['seq_len_bp'] += float(len_matches.group(2)) * bp
                            s['seq_len_read_count'] += float(len_matches.group(2))
                    
                    if in_module == 'seq_dup_levels':
                        if l[:1] == '#':
                            continue # headers
                        sections  = l.split()
                        self.fastqc_data['seq_dup_levels'][s_name][sections[0]] = float(sections[1])
                        self.fastqc_data['seq_dup_levels_dedup'][s_name][sections[0]] = float(sections[2])
        
                    if in_module == 'adapter_content':
                        if l[:1] == '#':
                            adapter_types = l[1:].split("\t")[1:]
                        else:
                            cols = l.split("\t")
                            pos = int(cols[0].split('-', 1)[0])
                            for idx, val in enumerate(cols[1:]):
                                a = adapter_types[idx]
                                self.fastqc_data['adapter_content'][s_name][a][pos] = float(val)
                    
                    if in_module == 'kmer_content':
                        # fastqc_data.txt doesn't have the breakdown of this
                        # per-base, which is what we need to replicate the
                        # graph in the FastQC report. Hmmm...
                        pass #TODO
                    
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
            'modify': lambda x: x / 1000000
        }
        self.general_stats_addcols(self.fastqc_stats, headers)


    def sequence_quality_plot (self):
        """ Create the HTML for the phred quality score plot """
        if 'sequence_quality' not in self.fastqc_data or len(self.fastqc_data['sequence_quality']) == 0:
            log.debug('sequence_quality not found in FastQC reports')
            return None
        
        # Setup for multiple datasets
        datasets = [
            {'mean':            'Mean'},
            {'median':          'Median'},
            {'90_percentile':   '90<sup>th</sup> Percentile'},
            {'upper_quart':     'Upper Quartile'},
            {'lower_quart':     'Lower Quartile'},
            {'10_percentile':   '10<sup>th</sup> Percentile'},
        ]
        pdata = [self.fastqc_data['sequence_quality'][dtype] for d in datasets for dtype in d.keys()]
        
        # Find the max phred score, so that we keep the y_max consistent
        max_val = 35
        for dtype in pdata:
            for s_name in dtype:
                for val in dtype[s_name].values():
                    max_val = max(max_val, val)        
        
        pconfig = {
            'data_labels': [{'name': name} for d in datasets for name in d.values()],
            'id': 'fastqc_sequence_quality_plot',
            'title': 'Mean Quality Scores',
            'ylab': 'Phred Score',
            'xlab': 'Position (bp)',
            'ymax': max_val,
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
            'anchor': 'sequence-quality',
            'content': self.plot_xy_data(pdata, pconfig)
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
            'xPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.sections.append({
            'name': 'Per Sequence Quality Scores',
            'anchor': 'per-seq-quality',
            'content': self.plot_xy_data(self.fastqc_data['per_seq_quality'], pconfig)
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
            'colors': self.get_status_cols('per_seq_quality'),
        }
        self.sections.append({
            'name': 'Per Sequence GC Content',
            'anchor': 'gc-content',
            'content': self.plot_xy_data(self.fastqc_data['gc_content'], pconfig)
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
            'ymax': 100,
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': self.get_status_cols('n_content'),
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ]
        }
        self.sections.append({
            'name': 'Per Base N Content',
            'anchor': 'n-content',
            'content': self.plot_xy_data(self.fastqc_data['n_content'], pconfig)
        })
    

    def sequence_content_plot (self):
        """ Create the epic HTML for the FastQC sequence content heatmap """
        if 'sequence_content' not in self.fastqc_data or len(self.fastqc_data['sequence_content']) == 0:
            log.debug('sequence_content not found in FastQC reports')
            return None
        
        self.status_classes = {
            'pass': 'label-success',
            'warn': 'label-warning',
            'fail': 'label-danger',
            'default': 'label-default'
        }
        
        # Get the sample statuses
        data = dict()
        names = list()
        for s in sorted(self.seq_content):
            names.append(s)
            data[s] = self.seq_content[s]

        # Order the table by the sample names
        data = OrderedDict(sorted(data.items()))
        
        images = [{'s_name': s, 'img_path': 'report_data/fastqc/{}_per_base_sequence_content.png'.format(s)}
                    for s in sorted(self.seq_content.keys())]
        statuses = {s: self.statuses['fastqc_seq'][s] for s in self.statuses['fastqc_seq'].keys()}
        
        # Order the table by the sample names
        data = OrderedDict(sorted(data.items()))
        
        if len(names) > 1:
            next_prev_buttons = '<div class="clearfix"><div class="btn-group btn-group-sm"> \n\
                <a href="#{p_n}" class="btn btn-default original_plot_prev_btn" data-target="#fastqc_seq">&laquo; Previous</a> \n\
                <a href="#{n_n}" class="btn btn-default original_plot_nxt_btn" data-target="#fastqc_seq">Next &raquo;</a> \n\
            </div></div>'.format(p_n=names[-1], n_n=names[1])
        else: next_prev_buttons = ''
        
        html = '<div id="fastqc_sequence_content_plot"> \n\
            <h4><span class="s_name">{fn}</span> <span class="label {status_class} s_status">{this_status}</span></h4> \n\
            <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">\n\
                <div id="fastqc_seq" class="hc-plot"> \n\
                    <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas> \n\
                </div> \n\
                <ul id="fastqc_seq_heatmap_key">\n\
                    <li>Position: <span id="fastqc_seq_heatmap_key_pos"></span></li> \n\
                    <li>%T: <span id="fastqc_seq_heatmap_key_t"></span> <span id="fastqc_seq_heatmap_key_colourbar_t" class="heatmap_colourbar"><span></span></span></li>\n\
                    <li>%C: <span id="fastqc_seq_heatmap_key_c"></span> <span id="fastqc_seq_heatmap_key_colourbar_c" class="heatmap_colourbar"><span></span></span></li>\n\
                    <li>%A: <span id="fastqc_seq_heatmap_key_a"></span> <span id="fastqc_seq_heatmap_key_colourbar_a" class="heatmap_colourbar"><span></span></span></li>\n\
                    <li>%G: <span id="fastqc_seq_heatmap_key_g"></span> <span id="fastqc_seq_heatmap_key_colourbar_g" class="heatmap_colourbar"><span></span></span></li>\n\
                    <li><small class="text-muted">Values are approximate</small></li>\n\
                </ul>\n\
            </div> \n\
            <div class="clearfix"></div> \n\
        </div> \n\
        <script type="text/javascript"> \n\
            fastqc_seq_content_data = {d};\n\
            if(typeof fastqc_s_names == "undefined"){{ fastqc_s_names = []; }} \n\
            fastqc_s_names["fastqc_seq"] = {n};\n\
            if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
            fastqc_s_statuses["fastqc_seq"] = {s};\n\
            $(function () {{ \n\
                fastqc_seq_content_heatmap(); \n\
            }}); \n\
        </script>'.format(b=next_prev_buttons, fn=names[0], d=json.dumps(data), n=json.dumps(names), this_status=statuses[names[0]], status_class=self.status_classes.get(statuses[names[0]], 'label-default'), s=json.dumps(statuses))
        
        self.sections.append({
            'name': 'Per Base Sequence Content',
            'anchor': 'sequence-content',
            'content': html
        })

    def adapter_content_plot (self):
        """ Create the HTML for the FastQC adapter plot """
        if 'adapter_content' not in self.fastqc_data or len(self.fastqc_data['adapter_content']) == 0:
            log.debug('adapter_content not found in FastQC reports')
            return None
        
        # Check that there is some adapter contamination in some of the plots
        max_val = 0
        for s in self.fastqc_data['adapter_content'].keys():
            for v in self.fastqc_data['adapter_content'][s].values():
                max_val = max(max_val, v)
        if max_val <= 0.1:
            return '<p>No adapter contamination found in any samples.</p>'
        
        pconfig = {
            'id': 'fastqc_adapter_content_plot',
            'title': 'Adapter Content',
            'ylab': '% of Sequences',
            'xlab': 'Position',
            'ymax': 100,
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'hide_empty': True
        }
        # Note: No point in adding colour by status here. If there's anything to
        # show, it's usually a fail. So everything is red. Boring.
        
        self.sections.append({
            'name': 'Adapter Content',
            'anchor': 'adapter-content',
            'content': self.plot_xy_data(self.fastqc_data['adapter_content'], pconfig)
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
        colours = list()
        for s_name, status in self.fastqc_statuses[key].items():
            colours.append( self.status_colours.get(status, self.status_colours['default']) )
        return colours
    