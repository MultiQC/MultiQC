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
from collections import OrderedDict
import io
import json
import logging
import os
import re
import shutil
import zipfile

import multiqc
from multiqc import config

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('FastQC'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(log)

        # Static variables
        self.name = "FastQC"
        self.anchor = "fastqc"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">FastQC</a> \
            is a quality control tool for high throughput sequence data, written by Simon Andrews at the Babraham Institute in Cambridge.</p> \
            <p class="text-muted"><span class="glyphicon glyphicon-info-sign" aria-hidden="true"></span> <strong>Tip:</strong> Click a data point to view the original FastQC plot.</p>'
        self.data_dir = os.path.join(config.output_dir, 'report_data', 'fastqc')

        # Find and load any FastQC reports
        fastqc_raw_data = {}
        plot_fns = [
            'per_base_quality.png',
            'per_base_sequence_content.png',
            'per_sequence_gc_content.png',
            'adapter_content.png'
        ]
        for root, dirnames, filenames in os.walk(config.analysis_dir, followlinks=True):
            # Extracted FastQC directory
            if 'fastqc_data.txt' in filenames:
                s_name = os.path.basename(root)
                d_path = os.path.join(root, 'fastqc_data.txt')
                with io.open (d_path, "r", encoding='utf-8') as f:
                    r_data = f.read()

                # Get the sample name from inside the file if possible
                fn_search = re.search(r"^Filename\s+(.+)$", r_data, re.MULTILINE)
                if fn_search:
                    s_name = fn_search.group(1).strip()
                s_name = self.clean_s_name(s_name, root)
                if s_name in fastqc_raw_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                fastqc_raw_data[s_name] = r_data

                # Copy across the raw images
                if not os.path.exists(self.data_dir):
                    os.makedirs(self.data_dir)
                for p in plot_fns:
                    try:
                        shutil.copyfile(os.path.join(root, 'Images', p), os.path.join(self.data_dir, "{}_{}".format(s_name, p)))
                    except IOError:
                        pass

            # Zipped FastQC report
            for f in filenames:
                if f[-11:] == '_fastqc.zip':
                    s_name = f[:-11]
                    fqc_zip = zipfile.ZipFile(os.path.join(root, f))
                    # FastQC zip files should have one directory inside, containing report
                    d_name = fqc_zip.namelist()[0]
                    try:
                        with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as f:
                            r_data = f.read().decode('utf8')

                        # Get the sample name from inside the file if possible
                        fn_search = re.search(r"^Filename\s+(.+)$", r_data, re.MULTILINE)
                        if fn_search:
                            s_name = fn_search.group(1).strip()
                        s_name = self.clean_s_name(s_name, root)
                        if s_name in fastqc_raw_data:
                            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                        fastqc_raw_data[s_name] = r_data

                    except KeyError:
                        log.warning("Error - can't find fastqc_raw_data.txt in {}".format(f))
                    else:
                        # Copy across the raw images
                        if not os.path.exists(self.data_dir):
                            os.makedirs(self.data_dir)
                        for p in plot_fns:
                            try:
                                with fqc_zip.open(os.path.join(d_name, 'Images', p)) as f:
                                    img = f.read()
                                with io.open (os.path.join(self.data_dir, "{}_{}".format(s_name, p)), "wb") as f:
                                    f.write(img)
                            except KeyError:
                                pass

        if len(fastqc_raw_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning


        log.info("Found {} reports".format(len(fastqc_raw_data)))

        self.sections = list()
        
        self.statuses = {
            'fastqc_quals' : {},
            'fastqc_gc' : {},
            'fastqc_seq' : {},
            'fastqc_adapter' : {},
        }
        self.status_colours = {
            'pass': '#5cb85c',
            'warn': '#f0ad4e',
            'fail': '#d9534f',
            'default': '#999'
        }
        self.status_classes = {
            'pass': 'label-success',
            'warn': 'label-warning',
            'fail': 'label-danger',
            'default': 'label-default'
        }

        # Get the section pass and fails
        passfails = self.fastqc_get_passfails(fastqc_raw_data)
        self.intro += '<script type="text/javascript">fastqc_passfails = {};</script>'.format(json.dumps(passfails))

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        parsed_stats = self.fastqc_general_stats(fastqc_raw_data)
        self.fastqc_stats_table(parsed_stats)
        
        # Write the basic stats table data to a file
        with io.open (os.path.join(config.output_dir, 'report_data', 'multiqc_fastqc.txt'), "w", encoding='utf-8') as f:
            print( self.dict_to_csv( parsed_stats ), file=f)


        # Section 1 - Quality Histograms
        self.parse_fastqc_seq_quality(fastqc_raw_data)
        if len(self.sequence_quality) > 0:
            self.sections.append({
                'name': 'Sequence Quality Histograms',
                'anchor': 'sequence-quality',
                'content': self.fastqc_quality_overlay_plot()
            })

        # Section 2 - GC Content
        self.parse_fastqc_gc_content(fastqc_raw_data)
        if len(self.gc_content) > 0:
            self.sections.append({
                'name': 'Per Sequence GC Content',
                'anchor': 'gc-content',
                'content': self.fastqc_gc_overlay_plot()
            })

        # Section 3 - Per-base sequence content
        self.parse_fastqc_seq_content(fastqc_raw_data)
        if len(self.seq_content) > 0:
            self.sections.append({
                'name': 'Per Base Sequence Content',
                'anchor': 'sequence-content',
                'content': self.fastqc_seq_heatmap()
            })

        # Section 4 - Adapter Content
        self.fastqc_adapter_content(fastqc_raw_data)
        if len(self.adapter_content) > 0:
            self.sections.append({
                'name': 'Adapter Content',
                'anchor': 'adapter-content',
                'content': self.fastqc_adapter_overlay_plot()
            })


        # Copy across the required module files (CSS / Javascript etc)
        self.init_modfiles()


    def fastqc_get_passfails(self, fastqc_raw_data):
        """ Go through the headings for each report and count how many
        of each module were passes, fails and warnings. Returns a 2D dict,
        first key is section name, second keys are passes / fails / warnings
        and html (html uses bootstrap labels) """

        parsed_data = {
            'sequence-quality': {'pattern': '>>Per base sequence quality\s+(pass|warn|fail)'},
            'sequence-content': {'pattern': '>>Per base sequence content\s+(pass|warn|fail)'},
            'gc-content': {'pattern': '>>Per sequence GC content\s+(pass|warn|fail)'},
            'adapter-content': {'pattern': '>>Adapter Content\s+(pass|warn|fail)'},
        }
        counts = {'pass': 0, 'warn': 0, 'fail': 0}
        for p in parsed_data.keys():
            parsed_data[p].update(counts)

        for s, data in fastqc_raw_data.items():
            for p in parsed_data.keys():
                match = re.search(parsed_data[p]['pattern'], data)
                if match:
                    parsed_data[p][match.group(1)] += 1

        return parsed_data


    def fastqc_general_stats(self, fastqc_raw_data):
        """ Parse fastqc_data.txt for basic stats.
        Returns a 2D dict with basic stats, sample names as first keys,
        then statistic type as second key. """
        parsed_data = {}
        for s, data in fastqc_raw_data.items():
            parsed_data[s] = {}

            dups = re.search("#Total Deduplicated Percentage\s+([\d\.]+)", data)
            if dups:
                parsed_data[s]['percent_duplicates'] = "{0:.1f}%".format(100 - float(dups.group(1)))

            seqlen = re.search("Sequence length\s+([\d-]+)", data)
            if seqlen:
                parsed_data[s]['sequence_length'] = seqlen.group(1)

            gc = re.search("%GC\s+(\d+)", data)
            if gc:
                parsed_data[s]['percent_gc'] = "{}%".format(gc.group(1))

            numseq = re.search("Total Sequences\s+(\d+)", data)
            if numseq:
                parsed_data[s]['total_sequences_m'] = "{0:.1f}".format(float(numseq.group(1))/1000000)

            # Work out the average sequence length as the range is a bit useless
            if '-' in parsed_data[s]['sequence_length']:
                in_module = False
                total_count = 0
                total_bp = 0
                for l in data.splitlines():
                    if l[:30] == ">>Sequence Length Distribution":
                        in_module = True
                    elif l == ">>END_MODULE":
                        in_module = False
                    elif in_module is True:
                        len_matches = re.search("([\d-]+)\s+([\d\.E]+)", l)
                        try:
                            avg_len = len_matches.group(1)
                            if '-' in avg_len:
                                maxlen = float(avg_len.split("-",1)[1])
                                minlen = float(avg_len.split("-",1)[0])
                                avg_len = ((maxlen - minlen)/2) + minlen
                            avg_len = int(avg_len)
                            total_bp += float(len_matches.group(2)) * avg_len
                            total_count += float(len_matches.group(2))
                        except AttributeError:
                            pass
                if total_count > 0:
                    parsed_data[s]['sequence_length'] = '{:.0f}'.format(total_bp / total_count)

        return parsed_data

    def fastqc_stats_table(self, parsed_stats):
        """ Take the parsed stats from the FastQC report and add them to the
        basic stats table at the top of the report """

        config.general_stats['headers']['percent_duplicates'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="FastQC: %&nbsp;Duplicate Reads">% Dups</span></th>'
        config.general_stats['headers']['percent_gc'] = '<th class="chroma-col"  data-chroma-scale="PRGn" data-chroma-max="80" data-chroma-min="20"><span data-toggle="tooltip" title="FastQC: Average %&nbsp;GC Content">% GC</span></th>'
        config.general_stats['headers']['sequence_length'] = '<th class="chroma-col" data-chroma-scale="RdYlGn" data-chroma-min="0"><span data-toggle="tooltip" title="FastQC: Average Sequence Length (bp)">Length</span></th>'
        config.general_stats['headers']['total_sequences_m'] = '<th class="chroma-col" data-chroma-scale="Blues" data-chroma-min="0"><span data-toggle="tooltip" title="FastQC: Total Sequences (millions)">M Seqs</span></th>'

        rowcounts = { 'total_sequences_m' : 0, 'sequence_length' : 0,
            'percent_gc' : 0, 'percent_duplicates' : 0 }

        for samp, vals in parsed_stats.items():
            for k, v in vals.items():
                config.general_stats['rows'][samp][k] = '<td class="text-right">{}</td>'.format(v)
                rowcounts[k] += 1

        # Remove header if we don't have any filled cells for it (eg. % dups in older FastQC reports)
        for k in rowcounts.keys():
            if rowcounts[k] == 0:
                config.general_stats['headers'].pop(k, None)

    def parse_fastqc_seq_quality(self, fastqc_raw_data):
        """ Parse the 'Per base sequence quality' data from fastqc_data.txt
        Creates self.sequence_quality with a dict - keys positions and values
        mean phred score. Scope for several other values here. """

        self.sequence_quality = {}
        for s, data in fastqc_raw_data.items():
            self.sequence_quality[s] = {}
            in_module = False
            for l in data.splitlines():
                if l[:27] == ">>Per base sequence quality":
                    in_module = True
                    self.statuses['fastqc_quals'][s] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    quals = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                    try:
                        bp = int(quals.group(1).split('-', 1)[0])
                        self.sequence_quality[s][bp] = float(quals.group(2))                        
                        # parsed_data[s]['base'].append(quals.group(1))
                        # parsed_data[s]['mean'].append(float(quals.group(2)))
                        # parsed_data[s]['median'].append(float(quals.group(3)))
                        # parsed_data[s]['lower_quart'].append(float(quals.group(4)))
                        # parsed_data[s]['upper_quart'].append(float(quals.group(5)))
                        # parsed_data[s]['10_percentile'].append(float(quals.group(6)))
                        # parsed_data[s]['90_percentile'].append(float(quals.group(7)))
                    except AttributeError:
                        pass
            if len(self.sequence_quality[s]) == 0:
                self.sequence_quality.pop(s, None)

    def fastqc_quality_overlay_plot (self):
        """ Create the HTML for the phred quality score plot """
        
        # Statuses
        statuses = {}
        s_colours = {}
        for s_name, status in self.statuses['fastqc_quals'].items():
            statuses[s_name] = status
            s_colours[s_name] = self.status_colours.get(status, self.status_colours['default'])
        
        pconfig = {
            'id': 'fastqc_quality_plot',
            'title': 'Mean Quality Scores',
            'ylab': 'Phred Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            'colors': s_colours,
            'plotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        
        # Original images
        images = [{'s_name': s, 'img_path': 'report_data/fastqc/{}_per_base_quality.png'.format(s)}
                    for s in sorted(self.sequence_quality.keys())]
        
        html = self.plot_xy_data(self.sequence_quality, pconfig, images)
        
        html += '<script type="text/javascript"> \n\
                    if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
                    fastqc_s_statuses["fastqc_quality_plot"] = {}; \n\
                </script>'.format(json.dumps(statuses))
        
        return html


    def parse_fastqc_gc_content(self, fastqc_raw_data):
        """ Parse the 'Per sequence GC content' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of with keys
        containing percentage and values containing counts """

        self.gc_content = {}
        for s, data in fastqc_raw_data.items():
            self.gc_content[s] = {}
            in_module = False
            for l in data.splitlines():
                if l[:25] == ">>Per sequence GC content":
                    in_module = True
                    self.statuses['fastqc_gc'][s] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    gc_matches = re.search("([\d]+)\s+([\d\.E]+)", l)
                    try:
                        self.gc_content[s][int(gc_matches.group(1))] = float(gc_matches.group(2))
                    except AttributeError:
                        pass
            if len(self.gc_content[s]) == 0:
                self.gc_content.pop(s, None)


    def fastqc_gc_overlay_plot (self):
        """ Create the HTML for the FastQC GC content plot """
        
        # Statuses
        statuses = {}
        s_colours = {}
        for s_name, status in self.statuses['fastqc_gc'].items():
            statuses[s_name] = status
            s_colours[s_name] = self.status_colours.get(status, self.status_colours['default'])
        
        pconfig = {
            'id': 'fastqc_gcontent_plot',
            'title': 'Per Sequence GC Content',
            'ylab': 'Count',
            'xlab': '%GC',
            'ymin': 0,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            'colors': s_colours
        }
        images = [{'s_name': s, 'img_path': 'report_data/fastqc/{}_per_sequence_gc_content.png'.format(s)}
                    for s in sorted(self.gc_content.keys())]
        
        html = self.plot_xy_data(self.gc_content, pconfig, images)
        
        html += '<script type="text/javascript"> \n\
                    if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
                    fastqc_s_statuses["fastqc_gcontent_plot"] = {}; \n\
                </script>'.format(json.dumps(statuses))
        
        return html



    def parse_fastqc_seq_content(self, fastqc_raw_data):
        """ Parse the 'Per base sequence content' data from fastqc_data.txt
        Returns a 3D dict, sample names as first keys, second key the base
        position and third key with the base ([ACTG]). Values contain percentages """

        self.seq_content = {}
        for s, data in fastqc_raw_data.items():
            self.seq_content[s] = {}
            in_module = False
            for l in data.splitlines():
                if l[:27] == ">>Per base sequence content":
                    in_module = True
                    self.statuses['fastqc_seq'][s] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    l.replace('NaN','0')
                    seq_matches = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                    try:
                        bp = int(seq_matches.group(1).split('-', 1)[0])
                        self.seq_content[s][bp] = {}
                        self.seq_content[s][bp]['base'] = seq_matches.group(1)
                        self.seq_content[s][bp]['G'] = float(seq_matches.group(2))
                        self.seq_content[s][bp]['A'] = float(seq_matches.group(3))
                        self.seq_content[s][bp]['T'] = float(seq_matches.group(4))
                        self.seq_content[s][bp]['C'] = float(seq_matches.group(5))
                    except AttributeError:
                        if l[:1] != '#':
                            log.debug("Couldn't parse a line from sequence content report {}: {}".format(s, l))
            if len(self.seq_content[s]) == 0:
                self.seq_content.pop(s, None)

    def fastqc_seq_heatmap (self):
        """ Create the epic HTML for the FastQC sequence content heatmap """
        
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
        
        html = '<p class="text-muted instr">Click to show original FastQC plot.</p>\n\
        <div id="fastqc_seq"> \n\
            <h4><span class="s_name">{fn}</span> <span class="label {status_class} s_status">{this_status}</span></h4> \n\
            <div class="showhide_orig" style="display:none;"> \n\
                {b} <img data-toggle="tooltip" title="Click to return to overlay plot" class="original-plot" src="report_data/fastqc/{fn}_per_base_sequence_content.png"> \n\
            </div>\n\
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
            var fastqc_seq_orig_plots = {oplots};\n\
            $(function () {{ \n\
                fastqc_seq_content_heatmap(); \n\
            }}); \n\
        </script>'.format(b=next_prev_buttons, fn=names[0], d=json.dumps(data), n=json.dumps(names), this_status=statuses[names[0]], status_class=self.status_classes.get(statuses[names[0]], 'label-default'), s=json.dumps(statuses), oplots=json.dumps(images))
        
        return html


    def fastqc_adapter_content(self, fastqc_raw_data):
        """ Parse the 'Adapter Content' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of with keys
        containing adapter type and values containing counts """

        self.adapter_content = {}
        for s, data in fastqc_raw_data.items():
            in_module = False
            adapter_types = []
            for l in data.splitlines():
                if l[:17] == ">>Adapter Content":
                    in_module = True
                    self.statuses['fastqc_adapter'][s] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    if l[:1] == '#':
                        adapter_types = l[1:].split("\t")[1:]
                        for a in adapter_types:
                            name = '{} - {}'.format(s,a)
                            self.adapter_content[name] = {}
                    else:
                        cols = l.split("\t")
                        pos = int(cols[0].split('-', 1)[0])
                        for idx, val in enumerate(cols[1:]):
                            name = '{} - {}'.format(s, adapter_types[idx])
                            self.adapter_content[name][pos] = float(val)
            for a in adapter_types:
                name = '{} - {}'.format(s,a)
                if len(self.adapter_content[name]) == 0:
                    self.adapter_content.pop(name, None)

    def fastqc_adapter_overlay_plot (self):
        """ Create the HTML for the FastQC adapter plot """
        
        # Check that there is some adapter contamination in some of the plots
        max_val = 0
        for s in self.adapter_content.keys():
            for v in self.adapter_content[s].values():
                max_val = max(max_val, v)
        if max_val <= 0.1:
            # Delete original plots - can't see them withoit anything to click on anyway
            adapter_plots = [ f for f in os.listdir(self.data_dir) if f.endswith("_adapter_content.png") ]
            for f in adapter_plots:
                os.remove(os.path.join(self.data_dir, f))
            return '<p>No adapter contamination found in any samples.</p>'
        
        pconfig = {
            'id': 'fastqc_adapter_plot',
            'title': 'Adapter Content',
            'ylab': '% of Sequences',
            'xlab': 'Position',
            'ymax': 100,
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'hide_empty': True
        }
        # NB: No point in adding colour by status here. If there's anything to
        # show, it's usually a fail. So everything is red. Boring.
        
        images = []
        samps = []
        for s in sorted(self.adapter_content.keys()):
            s_name = s.split(" - ")[0];
            if s_name not in samps:
                samps.append(s_name)
                images.append({
                    's_name': s_name,
                    'img_path': 'report_data/fastqc/{}_adapter_content.png'.format(s_name)
                })
        
        html = self.plot_xy_data(self.adapter_content, pconfig, images)
        
        # Make a JS variable holding the FastQC status for each sample
        statuses = {s: self.statuses['fastqc_adapter'][s] for s in self.statuses['fastqc_adapter'].keys()}
        html += '<script type="text/javascript"> \n\
                    if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
                    fastqc_s_statuses["fastqc_adapter_plot"] = {}; \n\
                </script>'.format(json.dumps(statuses))
        
        return html


    def init_modfiles(self):
        """ Copy the required assets into the output directory.
        Need to do this manually as we don't want to over-write
        the existing directory tree."""

        self.css = [ os.path.join('assets', 'css', 'multiqc_fastqc.css') ]
        self.js = [ os.path.join('assets', 'js', 'multiqc_fastqc.js') ]

        for f in self.css + self.js:
            d = os.path.join(config.output_dir, os.path.dirname(f))
            if not os.path.exists(d):
                os.makedirs(d)
            if not os.path.exists(os.path.join(config.output_dir, f)):
                shutil.copy(os.path.join(os.path.dirname(__file__), f), os.path.join(config.output_dir, f))
