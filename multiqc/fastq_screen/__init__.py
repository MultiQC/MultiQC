#!/usr/bin/env python

""" MultiQC module to parse output from FastQ Screen """

from collections import OrderedDict
import json
import logging
import os
import re
import shutil

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "FastQ Screen"
        self.anchor = "fastq_screen"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/" target="_blank">FastQ Screen</a> \
            allows you to screen a library of sequences in FastQ format against a set of sequence databases so you can see if the \
            composition of the library matches with what you expect.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any FastQ Screen reports
        fq_screen_raw_data = {}
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            for fn in filenames:
                if fn.endswith("_screen.txt"):
                    s_name = fn[:-11]
                    with open (os.path.join(root,fn), "r") as f:
                        fq_screen_raw_data[s_name] = f.read()

        if len(fq_screen_raw_data) == 0:
            logging.debug("Could not find any FastQ Screen reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} FastQ Screen reports".format(len(fq_screen_raw_data)))

        self.sections = list()

        # Section 1 - Alignment Profiles
        fq_screen_data = self.parse_fqscreen(fq_screen_raw_data)
        self.intro += self.fqscreen_plot(fq_screen_data)


    def parse_fqscreen(self, fq_screen_raw_data):
        """ Parse the FastQ Screen output into a 3D dict """
        parsed_data = dict()
        for s, data in fq_screen_raw_data.iteritems():
            parsed_data[s] = OrderedDict()
            for l in data.splitlines():
                if l[:18] == "%Hit_no_libraries:":
                    org = 'No hits'
                    parsed_data[s][org] = {'percentages':{}}
                    parsed_data[s][org]['percentages']['one_hit_one_library'] = float(l[19:])
                    parsed_data[s][org]['percentages']['unmapped'] = 0
                    parsed_data[s][org]['percentages']['multiple_hits_one_library'] = 0
                    parsed_data[s][org]['percentages']['one_hit_multiple_libraries'] = 0
                    parsed_data[s][org]['percentages']['multiple_hits_multiple_libraries'] = 0
                else:
                    fqs = re.search(r"^(\w+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)$", l)
                    if fqs:
                        org = fqs.group(1)
                        parsed_data[s][org] = {'percentages':{}, 'counts':{}}
                        parsed_data[s][org]['counts']['reads_processed'] = int(fqs.group(2))
                        parsed_data[s][org]['counts']['unmapped'] = int(fqs.group(3))
                        parsed_data[s][org]['percentages']['unmapped'] = float(fqs.group(4))
                        parsed_data[s][org]['counts']['one_hit_one_library'] = int(fqs.group(5))
                        parsed_data[s][org]['percentages']['one_hit_one_library'] = float(fqs.group(6))
                        parsed_data[s][org]['counts']['multiple_hits_one_library'] = int(fqs.group(7))
                        parsed_data[s][org]['percentages']['multiple_hits_one_library'] = float(fqs.group(8))
                        parsed_data[s][org]['counts']['one_hit_multiple_libraries'] = int(fqs.group(9))
                        parsed_data[s][org]['percentages']['one_hit_multiple_libraries'] = float(fqs.group(10))
                        parsed_data[s][org]['counts']['multiple_hits_multiple_libraries'] = int(fqs.group(11))
                        parsed_data[s][org]['percentages']['multiple_hits_multiple_libraries'] = float(fqs.group(12))
            if len(parsed_data[s]) == 0:
                logging.warning("Could not parse FastQ Screen report {}".format(s))
                parsed_data.pop(s, None)

        return parsed_data

    def fqscreen_plot (self, parsed_data):
        categories = list()
        getCats = True
        data = list()
        p_types = OrderedDict()
        p_types['multiple_hits_multiple_libraries'] = {'col': '#7f0000', 'name': 'Multiple Hits, Multiple Libraries' }
        p_types['one_hit_multiple_libraries'] = {'col': '#ff0000', 'name': 'One Hit, Multiple Libraries' }
        p_types['multiple_hits_one_library'] = {'col': '#00007f', 'name': 'Multiple Hits, One Library' }
        p_types['one_hit_one_library'] = {'col': '#0000ff', 'name': 'One Hit, One Library' }
        for k, t in p_types.iteritems():
            first = True
            for s in sorted(parsed_data):
                thisdata = list()
                if len(categories) > 0:
                    getCats = False
                for org in parsed_data[s]:
                    thisdata.append(parsed_data[s][org]['percentages'][k])
                    if getCats:
                        categories.append(org)
                td = {
                    'name': t['name'],
                    'stack': s,
                    'data': thisdata,
                    'color': t['col']
                }
                if first:
                    first = False
                else:
                    td['linkedTo'] = ':previous'
                data.append(td)

        html = '<div id="fq_screen_plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            fq_screen_data = {};\n\
            fq_screen_categories = {};\n\
            $(function () {{ \n\
                $("#fq_screen_plot").highcharts({{ \n\
                    chart: {{ type: "column" }}, \n\
                    title: {{ text: "FastQ Screen Results" }}, \n\
                    xAxis: {{ categories: fq_screen_categories }}, \n\
                    yAxis: {{ \n\
                        max: 100, \n\
                        min: 0, \n\
                        title: {{ text: "Percentage Aligned" }} \n\
                    }}, \n\
                    tooltip: {{ \n\
                        formatter: function () {{ \n\
                            return "<b>" + this.series.stackKey.replace("column","") + " - " + this.x + "</b><br/>" + \n\
                                this.series.name + ": " + this.y + "%<br/>" + \n\
                                "Total Alignment: " + this.point.stackTotal + "%"; \n\
                        }}, \n\
                    }}, \n\
                    plotOptions: {{ \n\
                        column: {{ \n\
                            pointPadding: 0, \n\
                            groupPadding: 0.02, \n\
                            stacking: "normal" }} \n\
                    }}, \n\
                    series: fq_screen_data \n\
                }}); \n\
            }}); \n\
        </script>'.format(json.dumps(data), json.dumps(categories))

        return html
