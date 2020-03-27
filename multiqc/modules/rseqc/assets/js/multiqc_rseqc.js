// Return JS code required for plotting a single sample
// RSeQC plot. Attempt to make it look as much like the original as possible.
// Note: this code is injected by `eval(str)`, not <script type="text/javascript">


function single_sample_plot(e){
    // In case of repeated modules: #rseqc_junction_saturation_plot, #rseqc_junction_saturation_plot-1, ..
    var rseqc_junction_saturation_plot = $(e.currentTarget).closest('.hc-plot');
    var rseqc_junction_saturation_plot_id = rseqc_junction_saturation_plot.attr('id');
    var junction_sat_single_hint = rseqc_junction_saturation_plot.closest('.mqc-section').find('#junction_sat_single_hint');


    // Get the three datasets for this sample
    var data = [
        {'name': 'All Junctions'},
        {'name': 'Known Junctions'},
        {'name': 'Novel Junctions'}
    ];
    var k = 0;
    for (var i = 0; i < 3; i++) {
        var ds = mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][i];
        for (k = 0; k < ds.length; k++){
            if(ds[k]['name'] == this.series.name){
                data[i]['data'] = JSON.parse(JSON.stringify(ds[k]['data']));
                break;
            }
        }
    }

    // Create single plot div, and hide overview
    var newplot = $('<div id="rseqc_junction_saturation_single"> \
        <div id="rseqc_junction_saturation_single_controls"> \
          <button class="btn btn-primary btn-sm" id="rseqc-junction_sat_single_return"> \
            Return to overview \
          </button> \
          <div class="btn-group btn-group-sm"> \
            <button class="btn btn-default rseqc-junction_sat_single_prevnext" data-action="prev">&laquo; Prev</button> \
            <button class="btn btn-default rseqc-junction_sat_single_prevnext" data-action="next">Next &raquo;</button> \
          </div> \
        </div> \
        <div class="hc-plot-wrapper"> \
          <div class="hc-plot hc-line-plot"> \
            <small>loading..</small> \
          </div> \
        </div> \
      </div>');
    var pwrapper = rseqc_junction_saturation_plot.parent().parent();
    newplot.insertAfter(pwrapper).hide().slideDown();
    pwrapper.slideUp();
    junction_sat_single_hint.slideUp();

    // Listener to return to overview
    newplot.find('#rseqc-junction_sat_single_return').click(function(e){
      e.preventDefault();
      newplot.slideUp(function(){
        $(this).remove();
      });
      pwrapper.slideDown();
      junction_sat_single_hint.slideDown();
    });

    // Listeners for previous / next plot
    newplot.find('.rseqc-junction_sat_single_prevnext').click(function(e){
      e.preventDefault();
      if($(this).data('action') == 'prev'){
        k--;
        if(k < 0){
          k = mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][0].length - 1;
        }
      } else {
        k++;
        if(k >= mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][0].length){
          k = 0;
        }
      }
      var hc = newplot.find('.hc-plot').highcharts();
      for (var i = 0; i < 3; i++) {
          hc.series[i].setData(mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][i][k]['data'], false);
      }
      var ptitle = 'RSeQC Junction Saturation: '+mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][0][k]['name'];
      hc.setTitle({text: ptitle});
      hc.redraw({ duration: 200 });
    });

    // Plot the single data
    newplot.find('.hc-plot').highcharts({
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
        title: { text: 'Number of observed splicing junctions' },
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
        headerFormat: '<strong>{point.key}% of reads</strong><br/>'
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
