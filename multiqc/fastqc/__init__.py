#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

import collections
import json
import logging
import os
import re
import shutil
import zipfile

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "FastQC"
        self.anchor = "fastqc"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">FastQC</a> \
            is a quality control tool for high throughput sequence data, written by Simon Andrews at the Babraham Institute in Cambridge.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']
        self.data_dir = os.path.join(self.output_dir, 'report_data', 'fastqc')

        # Find and load any FastQC reports
        fastqc_raw_data = {}
        plot_fns = [
            'per_base_quality.png',
            'per_base_sequence_content.png',
            'per_sequence_gc_content.png',
            'adapter_content.png'
        ]
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            # Extracted FastQC directory
            if root[-7:] == '_fastqc' and 'fastqc_data.txt' in filenames:
                s_name = os.path.basename(root)
                s_name = s_name[:-7]
                d_path = os.path.join(root, 'fastqc_data.txt')
                with open (d_path, "r") as f:
                    fastqc_raw_data[s_name] = f.read()
                # Copy across the raw images
                if not os.path.exists(self.data_dir):
                    os.makedirs(self.data_dir)
                for p in plot_fns:
                    shutil.copyfile(os.path.join(root, 'Images', p), os.path.join(self.data_dir, "{}_{}".format(s_name, p)))

            # Zipped FastQC report
            for f in filenames:
                if f[-11:] == '_fastqc.zip':
                    s_name = f[:-11]
                    d_name = f[:-4]
                    fqc_zip = zipfile.ZipFile(os.path.join(root, f))
                    try:
                        with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as f:
                            fastqc_raw_data[s_name] = f.read()
                    except KeyError:
                        pass # Can't find fastqc_raw_data.txt in the zip file
                    else:
                        # Copy across the raw images
                        if not os.path.exists(self.data_dir):
                            os.makedirs(self.data_dir)
                        for p in plot_fns:
                            with fqc_zip.open(os.path.join(d_name, 'Images', p)) as f:
                                img = f.read()
                            with open (os.path.join(self.data_dir, "{}_{}".format(s_name, p)), "w") as f:
                                f.write(img)

        if len(fastqc_raw_data) == 0:
            logging.debug("Could not find any FastQC reports in {}".format(analysis_dir))
            return


        logging.info("Found {} FastQC reports".format(len(fastqc_raw_data)))

        self.sections = list()

        # Get the section pass and fails
        passfails = self.fastqc_get_passfails(fastqc_raw_data)
        self.intro += '<script type="text/javascript">fastqc_passfails = {};</script>'.format(json.dumps(passfails))

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        parsed_stats = self.fastqc_basic_stats(fastqc_raw_data)
        self.fastqc_stats_table(parsed_stats, report)

        # Section 1 - Quality Histograms
        histogram_data = self.fastqc_seq_quality(fastqc_raw_data)
        self.sections.append({
            'name': 'Sequence Quality Histograms',
            'anchor': 'sequence-quality',
            'content': self.fastqc_quality_overlay_plot(histogram_data)
        })

        # Section 2 - GC Content
        gc_data = self.fastqc_gc_content(fastqc_raw_data)
        self.sections.append({
            'name': 'Per Sequence GC Content',
            'anchor': 'gc-content',
            'content': self.fastqc_gc_overlay_plot(gc_data)
        })

        # Section 3 - Per-base sequence content
        seq_content = self.fastqc_seq_content(fastqc_raw_data)
        self.sections.append({
            'name': 'Per Base Sequence Content',
            'anchor': 'sequence-content',
            'content': self.fastqc_seq_heatmap(seq_content)
        })

        # Section 4 - Adapter Content
        adapter_data = self.fastqc_adapter_content(fastqc_raw_data)
        self.sections.append({
            'name': 'Adapter Content',
            'anchor': 'adapter-content',
            'content': self.fastqc_adapter_overlay_plot(adapter_data)
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

        for s, data in fastqc_raw_data.iteritems():
            for p in parsed_data.keys():
                match = re.search(parsed_data[p]['pattern'], data)
                if match:
                    parsed_data[p][match.group(1)] += 1

        return parsed_data


    def fastqc_basic_stats(self, fastqc_raw_data):
        """ Parse fastqc_data.txt for basic stats.
        Returns a 2D dict with basic stats, sample names as first keys,
        then statistic type as second key. """
        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
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
                                maxlen = int(avg_len.split("-",1)[1])
                                minlen = int(avg_len.split("-",1)[0])
                                avg_len = ((maxlen - minlen)/2) + minlen
                            total_bp += float(len_matches.group(2)) * avg_len
                            total_count += float(len_matches.group(2))
                        except AttributeError:
                            pass
                if total_count > 0:
                    parsed_data[s]['sequence_length'] = '{:.0f}'.format(total_bp / total_count)

        return parsed_data

    def fastqc_stats_table(self, parsed_stats, report):
        """ Take the parsed stats from the FastQC report and add them to the
        basic stats table at the top of the report """

        report['basic_stats']['headers']['total_sequences_m'] = '<th class="chroma-col" data-chroma-scale="Blues" data-chroma-min="0"><span data-toggle="tooltip" title="FastQC: Total Sequences (millions)">M Seqs</span></th>'
        report['basic_stats']['headers']['sequence_length'] = '<th class="chroma-col" data-chroma-scale="RdYlGn" data-chroma-min="0"><span data-toggle="tooltip" title="FastQC: Average Sequence Length (bp)">Length</span></th>'
        report['basic_stats']['headers']['percent_gc'] = '<th class="chroma-col"  data-chroma-scale="PRGn" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="FastQC: Average % GC Content">%&nbsp;GC</span></th>'
        report['basic_stats']['headers']['percent_duplicates'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="FastQC: % Duplicate Reads">%&nbsp;Dups</span></th>'

        for samp, vals in parsed_stats.iteritems():
            for k, v in vals.iteritems():
                report['basic_stats']['rows'][samp][k] = '<td class="text-right">{}</td>'.format(v)

    def fastqc_seq_quality(self, fastqc_raw_data):
        """ Parse the 'Per base sequence quality' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of lists
        containing base, mean, median, lower_quart, upper_quart, 10_percentile
        and 90_percentile. """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {}
            parsed_data[s]['status'] = ''
            parsed_data[s]['base'] = list()
            parsed_data[s]['mean'] = list()
            parsed_data[s]['median'] = list()
            parsed_data[s]['lower_quart'] = list()
            parsed_data[s]['upper_quart'] = list()
            parsed_data[s]['10_percentile'] = list()
            parsed_data[s]['90_percentile'] = list()
            in_module = False
            for l in data.splitlines():
                if l[:27] == ">>Per base sequence quality":
                    in_module = True
                    parsed_data[s]['status'] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    quals = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                    try:
                        parsed_data[s]['base'].append(quals.group(1))
                        parsed_data[s]['mean'].append(float(quals.group(2)))
                        parsed_data[s]['median'].append(float(quals.group(3)))
                        parsed_data[s]['lower_quart'].append(float(quals.group(4)))
                        parsed_data[s]['upper_quart'].append(float(quals.group(5)))
                        parsed_data[s]['10_percentile'].append(float(quals.group(6)))
                        parsed_data[s]['90_percentile'].append(float(quals.group(7)))
                    except AttributeError:
                        pass

        return parsed_data

    def fastqc_quality_overlay_plot (self, parsed_data):

        data = list()
        names = list()
        statuses = dict()
        for s in sorted(parsed_data):
            pairs = list()
            for k, p in enumerate(parsed_data[s]['base']):
                pairs.append([int(p.split('-', 1)[0]), parsed_data[s]['mean'][k]])
                statuses[s] = parsed_data[s]['status']
            data.append({
                'name': s,
                'data': pairs
            })
            names.append(s);

        html = '<div id="fastqc_qual_original" class="fastqc_orig"> \n\
            <p class="text-muted instr">Click to show original FastQC plot.</p>\n\
            <div class="showhide_orig" style="display:none;"> \n\
                <h4><span class="s_name">'+names[0]+'</span> <span class="label label-default s_status">status</span></h4> \n\
                <div class="btn-group btn-group-sm"> \n\
                    <a href="#'+names[-1]+'" class="btn btn-default fastqc_prev_btn" data-target="#fastqc_qual_original">&laquo; Previous</a> \n\
                    <a href="#'+names[1]+'" class="btn btn-default fastqc_nxt_btn" data-target="#fastqc_qual_original">Next &raquo;</a> \n\
                </div>\n\
                <p><img class="original-plot" src="report_data/fastqc/'+names[0]+'_per_base_quality.png" data-fnsuffix="_per_base_quality.png"></p> \n\
            </div>\n\
        </div>\n\
        <div id="fastqc_quality_overlay" class="fastqc-overlay-plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            fastqc_overlay_hist_data = {};\n\
            if(typeof fastqc_s_names == "undefined"){{ fastqc_s_names = []; }} \n\
            fastqc_s_names["fastqc_qual_original"] = {};\
            if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
            fastqc_s_statuses["fastqc_qual_original"] = {};\
            var quals_pconfig = {{ \n\
                "title": "Mean Quality Scores",\n\
                "ylab": "Phred Score",\n\
                "xlab": "Position (bp)",\n\
                "ymin": 0,\n\
                "tt_label": "Click to show original plot.<br><b>Base {{point.x}}</b>",\n\
                "use_legend": false,\n\
                "click_func": function () {{ \n\
                    fastqc_chg_original (this.series.name, \'#fastqc_qual_original\'); \n\
                    $("#fastqc_qual_original .showhide_orig").delay(100).slideDown(); \n\
                    $("#fastqc_quality_overlay").delay(100).slideUp(); \n\
                }} \n\
            }}; \n\
            $(function () {{ \
                plot_xy_line_graph("#fastqc_quality_overlay", fastqc_overlay_hist_data, quals_pconfig); \
            }}); \
        </script>'.format(json.dumps(data), json.dumps(names), json.dumps(statuses));

        return html


    def fastqc_gc_content(self, fastqc_raw_data):
        """ Parse the 'Per sequence GC content' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of with keys
        containing percentage and values containing counts """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {'status': '', 'vals': dict()}
            in_module = False
            for l in data.splitlines():
                if l[:25] == ">>Per sequence GC content":
                    in_module = True
                    parsed_data[s]['status'] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    gc_matches = re.search("([\d]+)\s+([\d\.E]+)", l)
                    try:
                        parsed_data[s]['vals'][int(gc_matches.group(1))] = float(gc_matches.group(2))
                    except AttributeError:
                        pass

        return parsed_data

    def fastqc_gc_overlay_plot (self, parsed_data):
        data = list()
        names = list()
        statuses = dict()
        for s in sorted(parsed_data):
            pairs = list()
            statuses[s] = parsed_data[s]['status']
            for k, p in iter(sorted(parsed_data[s]['vals'].iteritems())):
                pairs.append([int(k), p])
            data.append({
                'name': s,
                'data': pairs
            })
            names.append(s)

        html = '<div id="fastqc_gc_original" class="fastqc_orig"> \n\
            <p class="text-muted instr">Click to show original FastQC plot.</p>\n\
            <div class="showhide_orig" style="display:none;"> \n\
                <h4><span class="s_name">'+names[0]+'</span> <span class="label label-default s_status">status</span></h4> \n\
                <div class="btn-group btn-group-sm"> \n\
                    <a href="#'+names[-1]+'" class="btn btn-default fastqc_prev_btn" data-target="#fastqc_gc_original">&laquo; Previous</a> \n\
                    <a href="#'+names[1]+'" class="btn btn-default fastqc_nxt_btn" data-target="#fastqc_gc_original">Next &raquo;</a> \n\
                </div>\n\
                <p><img class="original-plot" src="report_data/fastqc/'+names[0]+'_per_sequence_gc_content.png" data-fnsuffix="_per_sequence_gc_content.png"></p> \n\
            </div>\n\
        </div>\n\
        <div id="fastqc_gc_overlay" class="fastqc-overlay-plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \
            var fastqc_overlay_gc_data = {};\
            if(typeof fastqc_s_names == "undefined"){{ fastqc_s_names = []; }} \n\
            fastqc_s_names["fastqc_gc_original"] = {};\
            if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
            fastqc_s_statuses["fastqc_gc_original"] = {};\
            var gc_pconfig = {{ \n\
                "title": "Per Sequence GC Content",\n\
                "ylab": "Count",\n\
                "xlab": "%GC",\n\
                "ymin": 0,\n\
                "xmax": 100,\n\
                "xmin": 0,\n\
                "tt_label": "Click to show original plot.<br><b>{{point.x}}% GC</b>",\n\
                "use_legend": false,\n\
                "click_func": function () {{ \n\
                    fastqc_chg_original (this.series.name, \'#fastqc_gc_original\'); \n\
                    $("#fastqc_gc_original .showhide_orig").delay(100).slideDown(); \n\
                    $("#fastqc_gc_overlay").delay(100).slideUp(); \n\
                }} \n\
            }}; \n\
            $(function () {{ \
                plot_xy_line_graph("#fastqc_gc_overlay", fastqc_overlay_gc_data, gc_pconfig); \
            }}); \
        </script>'.format(json.dumps(data), json.dumps(names), json.dumps(statuses));

        return html



    def fastqc_seq_content(self, fastqc_raw_data):
        """ Parse the 'Per base sequence content' data from fastqc_data.txt
        Returns a 3D dict, sample names as first keys, second key the base
        position and third key with the base ([ACTG]). Values contain percentages """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {'status': '', 'vals': dict()}
            in_module = False
            for l in data.splitlines():
                if l[:27] == ">>Per base sequence content":
                    in_module = True
                    parsed_data[s]['status'] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    seq_matches = re.search("([\d-]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", l)
                    try:
                        bp = int(seq_matches.group(1).split('-', 1)[0])
                        parsed_data[s]['vals'][bp] = {}
                        parsed_data[s]['vals'][bp]['base'] = seq_matches.group(1)
                        parsed_data[s]['vals'][bp]['G'] = float(seq_matches.group(2))
                        parsed_data[s]['vals'][bp]['A'] = float(seq_matches.group(3))
                        parsed_data[s]['vals'][bp]['T'] = float(seq_matches.group(4))
                        parsed_data[s]['vals'][bp]['C'] = float(seq_matches.group(5))
                    except AttributeError:
                        if l[:1] != '#':
                            raise

        return parsed_data

    def fastqc_seq_heatmap (self, parsed_data):

        # Get theh sample statuses
        data = dict()
        names = list()
        statuses = dict()
        for s in sorted(parsed_data):
            names.append(s)
            statuses[s] = parsed_data[s]['status']
            data[s] = parsed_data[s]['vals']

        # Order the table by the sample names
        data = collections.OrderedDict(sorted(data.items()))

        html = '<div id="fastqc_seq_original" class="fastqc_orig"> \n\
            <p class="text-muted instr">Click to show original FastQC plot.</p>\n\
            <h4><span class="s_name">'+names[0]+'</span> <span class="label label-default s_status">'+statuses[names[0]]+'</span></h4> \n\
            <div class="showhide_orig" style="display:none;"> \n\
                <div class="btn-group btn-group-sm"> \n\
                    <a href="#'+names[-1]+'" class="btn btn-default fastqc_prev_btn" data-target="#fastqc_seq_original">&laquo; Previous</a> \n\
                    <a href="#'+names[1]+'" class="btn btn-default fastqc_nxt_btn" data-target="#fastqc_seq_original">Next &raquo;</a> \n\
                </div>\n\
                <p><img class="original-plot" src="report_data/fastqc/'+names[0]+'_per_base_sequence_content.png" data-fnsuffix="_per_base_sequence_content.png"></p> \n\
            </div>\n\
        </div>\n\
        <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">\n\
            <canvas id="fastqc_seq_heatmap" height="300px" width="800px" style="width:100%;"></canvas> \n\
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
        <script type="text/javascript"> \n\
            fastqc_seq_content_data = {};\n\
            if(typeof fastqc_s_names == "undefined"){{ fastqc_s_names = []; }} \n\
            fastqc_s_names["fastqc_seq_original"] = {};\n\
            if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
            fastqc_s_statuses["fastqc_seq_original"] = {};\
            $(function () {{ \n\
                fastqc_seq_content_heatmap(fastqc_seq_content_data); \n\
            }}); \n\
        </script>'.format(json.dumps(data), json.dumps(names), json.dumps(statuses))

        return html


    def fastqc_adapter_content(self, fastqc_raw_data):
        """ Parse the 'Adapter Content' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of with keys
        containing adapter type and values containing counts """

        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {'status': '', 'vals': dict()}
            in_module = False
            adapter_types = []
            for l in data.splitlines():
                if l[:17] == ">>Adapter Content":
                    in_module = True
                    parsed_data[s]['status'] = l.split()[-1]
                elif l == ">>END_MODULE":
                    in_module = False
                elif in_module is True:
                    if l[:1] == '#':
                        adapter_types = l[1:].split("\t")[1:]
                        for a in adapter_types:
                            parsed_data[s]['vals'][a] = dict()
                    else:
                        cols = l.split("\t")
                        pos = int(cols[0].split('-', 1)[0])
                        for idx, val in enumerate(cols[1:]):
                            parsed_data[s]['vals'][adapter_types[idx]][pos] = float(val)
        return parsed_data

    def fastqc_adapter_overlay_plot (self, parsed_data):

        data = list()
        names = list()
        statuses = dict()
        for s in sorted(parsed_data):
            statuses[s] = parsed_data[s]['status']
            for a, d in parsed_data[s]['vals'].iteritems():
                pairs = list()
                for base, count in iter(sorted(d.iteritems())):
                    pairs.append([base, count])
                data.append({
                    'name': '{} - {}'.format(s, a),
                    'data': pairs
                })
            names.append(s)

        html = '<div id="fastqc_adapter_original" class="fastqc_orig"> \n\
            <p class="text-muted instr">Click to show original FastQC plot.</p>\n\
            <div class="showhide_orig" style="display:none;"> \n\
                <h4><span class="s_name">'+names[0]+'</span> <span class="label label-default s_status">status</span></h4> \n\
                <div class="btn-group btn-group-sm"> \n\
                    <a href="#'+names[-1]+'" class="btn btn-default fastqc_prev_btn" data-target="#fastqc_adapter_original">&laquo; Previous</a> \n\
                    <a href="#'+names[1]+'" class="btn btn-default fastqc_nxt_btn" data-target="#fastqc_adapter_original">Next &raquo;</a> \n\
                </div>\n\
                <p><img class="original-plot" src="report_data/fastqc/'+names[0]+'_adapter_content.png" data-fnsuffix="_adapter_content.png"></p> \n\
            </div>\n\
        </div>\n\
        <div id="fastqc_adapter_overlay" class="fastqc-overlay-plot" style="height:500px;"></div>\n\
        <script type="text/javascript"> \
            fastqc_adapter_data = {};\
            if(typeof fastqc_s_names == "undefined"){{ fastqc_s_names = []; }} \n\
            fastqc_s_names["fastqc_adapter_original"] = {};\n\
            if(typeof fastqc_s_statuses == "undefined"){{ fastqc_s_statuses = []; }} \n\
            fastqc_s_statuses["fastqc_adapter_original"] = {};\
            var adapter_pconfig = {{ \n\
                "title": "Adapter Content",\n\
                "ylab": "% of Sequences",\n\
                "xlab": "Position",\n\
                "ymax": 100,\n\
                "ymin": 0,\n\
                "tt_label": "Click to show original plot.<br><b>Base {{point.x}}</b>",\n\
                "use_legend": false,\n\
                "click_func": function () {{ \n\
                    var snames = this.series.name.split(" - ");\n\
                    fastqc_chg_original (snames[0], \'#fastqc_adapter_original\'); \n\
                    $("#fastqc_adapter_original .showhide_orig").delay(100).slideDown(); \n\
                    $("#fastqc_adapter_overlay").delay(100).slideUp(); \n\
                }} \n\
            }}; \n\
            $(function () {{ \
                plot_xy_line_graph("#fastqc_adapter_overlay", fastqc_adapter_data, adapter_pconfig); \
            }}); \
        </script>'.format(json.dumps(data), json.dumps(names), json.dumps(statuses));

        return html


    def init_modfiles(self):
        """ Copy the required assets into the output directory.
        Need to do this manually as we don't want to over-write
        the existing directory tree."""

        self.css = [ os.path.join('assets', 'css', 'multiqc_fastqc.css') ]
        self.js = [ os.path.join('assets', 'js', 'multiqc_fastqc.js') ]

        for f in self.css + self.js:
            d = os.path.join(self.output_dir, os.path.dirname(f))
            if not os.path.exists(d):
                os.makedirs(d)
            if not os.path.exists(os.path.join(self.output_dir, f)):
                shutil.copy(os.path.join(os.path.dirname(__file__), f), os.path.join(self.output_dir, f))
