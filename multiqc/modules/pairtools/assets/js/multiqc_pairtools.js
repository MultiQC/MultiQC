
function get_current_name(name) {
  // some local vars to make the code readable
  const renameFrom = window.mqc_rename_f_texts;
  const renameTo = window.mqc_rename_t_texts;
  const regexp_mode = window.mqc_rename_regex_mode;
  var current_name = name;
  // apply all renamings to the current_name (if there are any):
  if(renameFrom.length > 0){
    $.each(renameFrom, function(idx, f_text){
      if(regexp_mode){
        var re = new RegExp(f_text,"g");
        current_name = current_name.replace(re, renameTo[idx]);
      } else {
        current_name = current_name.replace(f_text, renameTo[idx]);
      }
    });
  }
  // return after renaming applied (if any)
  return current_name;
}


function single_scaling(e) {

    // In case of repeated modules: #rseqc_junction_saturation_plot, #rseqc_junction_saturation_plot-1, ..
    var rseqc_junction_saturation_plot = $(e.currentTarget).closest('.hc-plot');
    var rseqc_junction_saturation_plot_id = rseqc_junction_saturation_plot.attr('id');
    // var junction_sat_single_hint = rseqc_junction_saturation_plot.closest('.mqc-section').find('#junction_sat_single_hint');

    // Get the 4 datasets for this sample
    var data = [
        {'name': 'FF'},
        {'name': 'RF'},
        {'name': 'FR'},
        {'name': 'RR'},
    ];
    var k = 0;
    for (var i = 1; i < 5; i++) {
        var ds = mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][i];
        for (k = 0; k < ds.length; k++){
            // this of course - is not correct when there are degenerate sample
            // names after renaming:
            let current_name = get_current_name(ds[k]['name']);
            if (current_name == this.series.name) {
                // making a copy of data to be able to mess with it 
                data[i-1]['data'] = JSON.parse(JSON.stringify(ds[k]['data']));
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
    // junction_sat_single_hint.slideUp();


    // Do something when something changed in filters ...
    $(document).on('mqc_highlights mqc_renamesamples mqc_hidesamples', function(e){
      e.preventDefault();
      newplot.slideUp(function(){
        $(this).remove();
      });
      pwrapper.slideDown();
      // junction_sat_single_hint.slideDown();
    });

    // Listener to return to overview
    newplot.find('#rseqc-junction_sat_single_return').click(function(e){
      e.preventDefault();
      newplot.slideUp(function(){
        $(this).remove();
      });
      pwrapper.slideDown();
      // junction_sat_single_hint.slideDown();
    });

    // Listeners for previous / next plot
    newplot.find('.rseqc-junction_sat_single_prevnext').click(function(e){
      e.preventDefault();
      if($(this).data('action') == 'prev'){
        k--;
        if(k < 0){
          k = mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][1].length - 1;
        }
      } else {
        k++;
        if(k >= mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][1].length){
          k = 0;
        }
      }
      var hc = newplot.find('.hc-plot').highcharts();
      for (var i = 1; i < 5; i++) {
          hc.series[i-1].setData(mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][i][k]['data'], false);
      }
      var current_name = get_current_name(mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][1][k]['name']);
      var ptitle = 'frequency of interactions: ' + current_name;
      hc.setTitle({text: ptitle});
      hc.redraw({ duration: 200 });
    });

    // Plot the single data
    newplot.find('.hc-plot').highcharts({
      chart: {
        type: 'line',
        zoomType: 'x'
      },
      colors: ['#e41a1c','#377eb8','#4daf4a','#984ea3'],
      title: {
        text: 'frequency of interaction for '+this.series.name,
        x: 30 // fudge to center over plot area rather than whole plot
      },
      xAxis: {
        type: 'logarithmic',
        title: { text: 'Genomic separation (bp)' },
        allowDecimals: false,
      },
      yAxis: {
        type: 'logarithmic',
        title: { text: 'frequency of interactions' },
        allowDecimal: true,
      },
      legend: {
        floating: true,
        layout: 'vertical',
        align: 'right',
        verticalAlign: 'top',
        x: 60,
        y: 40
      },
      tooltip: {
        shared: true,
        crosshairs: true,
        headerFormat: '<strong>genomic separation: {point.key} bp</strong><br/>'
      },
      plotOptions: {
        series: {
          animation: false,
          lineWidth: 4,
          marker: {
            lineColor: null,
            fillColor: 'transparent',
            lineWidth: 0,
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
