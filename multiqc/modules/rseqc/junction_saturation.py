#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC junction_saturation.py
http://rseqc.sourceforge.net/#junction-saturation-py """

from collections import OrderedDict
import logging
import re

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC junction_saturation frequency reports and parse their data """

    # Set up vars
    self.junction_saturation_all = dict()
    self.junction_saturation_known = dict()
    self.junction_saturation_novel = dict()
    self.junction_saturation_all_pct = dict()
    self.junction_saturation_known_pct = dict()
    self.junction_saturation_novel_pct = dict()

    # Go through files and parse data
    for f in self.find_log_files('rseqc/junction_saturation'):
        parsed = dict()
        for l in f['f'].splitlines():
            r = re.search(r"^([xyzw])=c\(([\d,]+)\)$", l)
            if r:
                parsed[r.group(1)] = [float(i) for i in r.group(2).split(',')]
        if len(parsed) == 4:
            if parsed['z'][-1] == 0:
                log.warn("Junction saturation data all zeroes, skipping: '{}'".format(f['s_name']))
            else:
                if f['s_name'] in self.junction_saturation_all:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f, section='junction_saturation')
                self.junction_saturation_all[f['s_name']] = OrderedDict()
                self.junction_saturation_known[f['s_name']] = OrderedDict()
                self.junction_saturation_novel[f['s_name']] = OrderedDict()
                for k, v in enumerate(parsed['x']):
                    self.junction_saturation_all[f['s_name']][v] = parsed['z'][k]
                    self.junction_saturation_known[f['s_name']][v] = parsed['y'][k]
                    self.junction_saturation_novel[f['s_name']][v] = parsed['w'][k]

    # Filter to strip out ignored sample names
    self.junction_saturation_all = self.ignore_samples(self.junction_saturation_all)
    self.junction_saturation_known = self.ignore_samples(self.junction_saturation_known)
    self.junction_saturation_novel = self.ignore_samples(self.junction_saturation_novel)

    if len(self.junction_saturation_all) > 0:

        # Make a normalised percentage version of the data
        for s_name in self.junction_saturation_all:
            self.junction_saturation_all_pct[s_name] = OrderedDict()
            self.junction_saturation_known_pct[s_name] = OrderedDict()
            self.junction_saturation_novel_pct[s_name] = OrderedDict()
            total = list(self.junction_saturation_all[s_name].values())[-1]
            if total == 0: # Shouldn't happen..?
                continue
            for k, v in self.junction_saturation_all[s_name].items():
                self.junction_saturation_all_pct[s_name][k] = (v/total)*100
                self.junction_saturation_known_pct[s_name][k] = (self.junction_saturation_known[s_name][k]/total)*100
                self.junction_saturation_novel_pct[s_name][k] = (self.junction_saturation_novel[s_name][k]/total)*100

        # Add line graph to section
        pconfig = {
            'id': 'rseqc_junction_saturation_plot',
            'title': 'RSeQC: Junction Saturation',
            'ylab': 'Percent of Junctions',
            'ymax': 100,
            'ymin': 0,
            'xlab': "Percent of reads",
            'xmin': 0,
            'xmax': 100,
            'tt_label': "<strong>{point.x}% of reads</strong>: {point.y:.2f}",
            'data_labels': [
                {'name': 'All Junctions'},
                {'name': 'Known Junctions'},
                {'name': 'Novel Junctions'},
            ],
            'cursor': 'pointer',
            'click_func': plot_single()
        }
        self.add_section (
            name = 'Junction Saturation',
            anchor = 'rseqc-junction_saturation',
            description = '<a href="http://rseqc.sourceforge.net/#junction-saturation-py" target="_blank">Junction Saturation</a>' \
                " calculates percentage of known splicing junctions that are observed" \
                " in each dataset. If sequencing depth is sufficient, all (annotated) splice junctions should" \
                " be rediscovered, resulting in a curve that reaches a plateau. Missing low abundance splice" \
                " junctions can affect downstream analysis. Counts are normalised to the total number of" \
                " observed junctions.</p>" \
                "<div class='alert alert-info' id='rseqc-junction_sat_single_hint'>" \
                "<span class='glyphicon glyphicon-hand-up'></span> Click a line" \
                " to see the data side by side (as in the original RSeQC plot).</div><p>",
            plot = linegraph.plot([
                    self.junction_saturation_all_pct,
                    self.junction_saturation_known_pct,
                    self.junction_saturation_novel_pct
                ], pconfig)
        )

    # Return number of samples found
    return len(self.junction_saturation_all)


def plot_single():
    """ Return JS code required for plotting a single sample
    RSeQC plot. Attempt to make it look as much like the original as possible. """

    return """
    function(e){

        // Get the three datasets for this sample
        var data = [
            {'name': 'All Junctions'},
            {'name': 'Known Junctions'},
            {'name': 'Novel Junctions'}
        ];
        for (var i = 0; i < 3; i++) {
            var ds = mqc_plots['rseqc_junction_saturation_plot']['datasets'][i];
            for (var k = 0; k < ds.length; k++){
                if(ds[k]['name'] == this.series.name){
                    data[i]['data'] = ds[k]['data'];
                    break;
                }
            }
        }

        // Create single plot div, and hide overview
        var newplot = '<div id="rseqc_junction_saturation_single">'+
            '<p><button class="btn btn-default" id="rseqc-junction_sat_single_return">'+
            'Return to overview</button></p><div class="hc-plot-wrapper">'+
            '<div class="hc-plot hc-line-plot"><small>loading..</small></div></div></div>';
        var pwrapper = $('#rseqc_junction_saturation_plot').parent().parent();
        $(newplot).insertAfter(pwrapper).hide().slideDown();
        pwrapper.slideUp();
        $('#rseqc-junction_sat_single_hint').slideUp();

        // Listener to return to overview
        $('#rseqc-junction_sat_single_return').click(function(e){
          e.preventDefault();
          $('#rseqc_junction_saturation_single').slideUp(function(){
            $(this).remove();
          });
          pwrapper.slideDown();
          $('#rseqc-junction_sat_single_hint').slideDown();
        });

        // Plot the single data
        $('#rseqc_junction_saturation_single .hc-plot').highcharts({
          chart: {
            type: 'line',
            zoomType: 'x'
          },
          colors: ['blue','red','green'],
          title: {
            text: 'RSeQC Junction Saturation: '+this.series.name,
            x: 30 // fudge to center over plot area rather than whole plot
          },
          xAxis: {
            title: { text: 'Percent of total reads' },
            allowDecimals: false,
          },
          yAxis: {
            title: { text: 'Percent of observed splicing junctions' },
            max: 100,
            min: 0,
          },
          legend: {
            floating: true,
            layout: 'vertical',
            align: 'left',
            verticalAlign: 'top',
            x: 60,
            y: 40
          },
          tooltip: {
            shared: true,
            crosshairs: true,
            headerFormat: '<strong>{point.key}% of reads</strong><br/>',
            valueDecimals: 2
          },
          plotOptions: {
            series: {
              animation: false,
              lineWidth: 1,
              marker: {
                lineColor: null,
                fillColor: 'transparent',
                lineWidth: 1,
                symbol: 'circle'
              },
            }
          },
          exporting: { buttons: { contextButton: {
            menuItems: window.HCDefaults.exporting.buttons.contextButton.menuItems,
            onclick: window.HCDefaults.exporting.buttons.contextButton.onclick
          } } },
          series: data
        });
    }
    """

