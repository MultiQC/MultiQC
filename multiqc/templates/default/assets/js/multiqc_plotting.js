////////////////////////////////////////////////
// HighCharts Plotting Code
////////////////////////////////////////////////

// Global plot data variable
mqc_plots = {};

// Initialise the toolbox filters
window.mqc_highlight_f_texts = [];
window.mqc_highlight_f_cols = [];
window.mqc_highlight_regex_mode = false;
window.mqc_rename_f_texts = [];
window.mqc_rename_t_texts = [];
window.mqc_rename_regex_mode = false;
window.mqc_hide_mode = 'hide';
window.mqc_hide_f_texts = [];
window.mqc_hide_regex_mode = false;
window.HCDefaults = undefined;

// Execute when page load has finished loading
$(function () {

  // Show loading warning
  $('.mqc_loading_warning').show();

  // Decompress the JSON plot data
  mqc_plots = JSON.parse(LZString.decompressFromBase64(mqc_compressed_plotdata));

  // HighCharts Defaults
  window.HCDefaults = $.extend(true, {}, Highcharts.getOptions(), {});
  Highcharts.setOptions({
    credits: {
      enabled: true,
      text: 'Created with MultiQC',
      href: 'http://multiqc.info'
    },
    lang: {
      decimalPoint: (mqc_config['decimalPoint_format'] == undefined ? '.' : mqc_config['decimalPoint_format']),
      thousandsSep: (mqc_config['thousandsSep_format'] == undefined ? ' ' : mqc_config['thousandsSep_format']),
    },
    exporting: {
      buttons: {
        contextButton: {
          menuItems: null,
          onclick: function () {
            // Tick only this plot in the toolbox and slide out
            $('#mqc_export_selectplots input').prop('checked', false);
            $('#mqc_export_selectplots input[value="'+this.renderTo.id+'"]').prop('checked', true);
            // Special case - Table scatter plots are in a modal, need to close this first
            if(this.renderTo.id == 'tableScatterPlot'){
              $('#tableScatterModal').modal('hide');
            }
            mqc_toolbox_openclose('#mqc_exportplots', true);
          },
          text: '<span style="color:#999999;">Export Plot</span>',
          symbol: 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAAXNSR0IArs4c6QAAAAlwSFlzAAALEwAACxMBAJqcGAAAAVlpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDUuNC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6dGlmZj0iaHR0cDovL25zLmFkb2JlLmNvbS90aWZmLzEuMC8iPgogICAgICAgICA8dGlmZjpPcmllbnRhdGlvbj4xPC90aWZmOk9yaWVudGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KTMInWQAAAXNJREFUOBHNUsuqwkAMPX2g4kJd+wOCuKgL//8btAXXIogvtOhCax9xzkBqveLg8gamHZKck+RMgP9mnquh5XIpaZrC8zx0Oh1EUfQ1P3QRkeR6vcL3fdxuN1cqnERhGIKHREEQOIl8V1RE0DyuXCeRC/g39iFeHMdSlqUV+HK5oCgKeyew1+vZEauqwnQ6fcN+aJTnObbbLdrtttWGL0bjiBT/fr9jMBhYX/PzxsrA4/EQ8+zY7/dotVpgdRoFZ5F+v4/ZbPaBCw+Hg8znc5s8Ho8J9kxV04DAxCwZg6aAJWEO7XQ6yWKxQJZlGI1Gr+fXEZhkls8zCTUZfexkMpmg2+2+dUMci1qNlKS5K0YjC0iSRDgSO1EfiblfxOmpxaaDr3Q8HqWpC1+NFbnhu91OSMKC5/OZ19pqIoq5Xq+xWq3qIAnoZxFdCQ3Sx65o9WisqsYENb0rofr1T3Iexi1qs9mIgjTp1z9JhsPhq/qvwG95Tw3FukJt8JteAAAAAElFTkSuQmCC)',
          symbolX: 23,
          symbolY: 19
        }
      }
    }
  });

  // Render plots on page load
  $('.hc-plot.not_rendered:visible:not(.gt_max_num_ds)').each(function(){
    var target = $(this).attr('id');
    // Only one point per dataset, so multiply limit by arbitrary number.
    var max_num = num_datasets_plot_limit * 50;
    // Deferring each plot call prevents browser from locking up
    setTimeout(function(){
        plot_graph(target, undefined, max_num);
        if($('.hc-plot.not_rendered:visible:not(.gt_max_num_ds)').length == 0){
          $('.mqc_loading_warning').hide();
        }
    }, 50);
  });
  if($('.hc-plot.not_rendered:visible:not(.gt_max_num_ds)').length == 0){
    $('.mqc_loading_warning').hide();
  }

  // Render a plot when clicked
  $('body').on('click', '.render_plot', function(e){
    var target = $(this).parent().attr('id');
    plot_graph(target);
    if($('.hc-plot.not_rendered').length == 0){
      $('#mqc-warning-many-samples').hide();
    }
  });

  // Render all plots from header
  $('#mqc-render-all-plots').click(function(){
    $('.hc-plot.not_rendered').each(function(){
      var target = $(this).attr('id');
      plot_graph(target);
    });
    $('#mqc-warning-many-samples').hide();
  });

  // Replot graphs when something changed in filters
  $(document).on('mqc_highlights mqc_renamesamples mqc_hidesamples', function(){
    // Replot graphs
    $('.hc-plot:not(.not_rendered)').each(function(){
      var target = $(this).attr('id');
      plot_graph(target);
    });
  });

  // Switch a HighCharts axis or data source
  $('.hc_switch_group button').click(function(e){
    e.preventDefault();
    $(this).siblings('button.active').removeClass('active');
    $(this).addClass('active');
    var target = $(this).data('target');
    var action = $(this).data('action');
    // Switch between values and percentages
    if(action == 'set_percent' || action == 'set_numbers'){
      var sym = (action == 'set_percent') ? '%' : '#';
      var stack_type = (action == 'set_percent') ? 'percent' : 'normal';
      mqc_plots[target]['config']['stacking'] = stack_type;
      mqc_plots[target]['config']['ytype'] = 'linear';
      plot_graph(target);
      var ylab = $(this).data('ylab');
      if(ylab != undefined){
        $('#'+target).highcharts().yAxis[0].setTitle({ text: ylab });
      }
      var xlab = $(this).data('xlab');
      if(xlab != undefined){
        $('#'+target).highcharts().xAxis[0].setTitle({ text: xlab });
      }
    }
    // Switch to log10 axis
    if(action == 'set_log'){
      mqc_plots[target]['config']['ytype'] = 'logarithmic';
      plot_graph(target);
    }
    // Switch data source
    if(action == 'set_data'){
      var ds = $(this).data('newdata');
      plot_graph(target, ds);
      var ylab = $(this).data('ylab');
      var xlab = $(this).data('xlab');
      var ymax = $(this).data('ymax');
      if(ylab != undefined){
        $('#'+target).highcharts().yAxis[0].setTitle({ text: ylab });
      }
      if(xlab != undefined){
        $('#'+target).highcharts().xAxis[0].setTitle({ text: xlab });
      }
      if(ymax != undefined){
        $('#'+target).highcharts().yAxis[0].setExtremes(null, ymax);
      }
    }
  });

  // Make HighCharts divs height-draggable
  // http://jsfiddle.net/Lkwb86c8/
  $('.hc-plot:not(.no-handle)').each(function(){
    if(!$(this).parent().hasClass('hc-plot-wrapper')){
      $(this).wrap('<div class="hc-plot-wrapper"></div>');
    }
    if(!$(this).siblings().hasClass('hc-plot-handle')){
      $(this).after('<div class="hc-plot-handle"><span></span><span></span><span></span></div>');
    }
    $(this).css({ height: 'auto', top: 0, bottom: '10px', position: 'absolute' });
  });
  $('.hc-plot-handle').on('mousedown', function(e){
    var wrapper = $(this).parent();
    var handle = $(this);
    var startHeight = wrapper.height();
    var pY = e.pageY;
    $(document).on('mouseup', function(e){
      // Clear listeners now that we've let go
      $(document).off('mousemove');
      $(document).off('mouseup');
      // Fire off a custom jQuery event for other javascript chunks to tie into
      // Bind to the plot div, which should have a custom ID
      $(wrapper.parent().find('.hc-plot, .beeswarm-plot')).trigger('mqc_plotresize');
    });
    $(document).on('mousemove', function(me){
      wrapper.css('height', startHeight + (me.pageY - pY));
    });
  });
  // Trigger HighCharts reflow when a plot is resized
  $('.hc-plot, .beeswarm-plot').on('mqc_plotresize', function(e){
    if($(this).highcharts()) {
      $(this).highcharts().reflow();
    }
  });

  // Switch a y axis limit on or off
  $('.mqc_hcplot_plotgroup').on('click', '.mqc_hcplot_yaxis_limit_toggle .mqc_switch_wrapper', function(){
    var target = $( $(this).data('target') ).highcharts();
    var ymax = $(this).data('ymax');
    var ymin = $(this).data('ymin');
    ymax = ymax == 'undefined' ? null : ymax;
    ymin = ymin == 'undefined' ? null : ymin;
    var mqc_switch = $(this).find('.mqc_switch');
    if(mqc_switch.hasClass('on')){
      target.yAxis[0].update({max: null, min:null});
      mqc_switch.removeClass('on').addClass('off').text('off');
    } else {
      target.yAxis[0].update({max: ymax, min: ymin});
      mqc_switch.removeClass('off').addClass('on').text('on');
    }
  });

  // Sort a heatmap by highlighted names
  $('.mqc_heatmap_sortHighlight').click(function(e){
    e.preventDefault();
    var target = $(this).data('target').substr(1);
    if(mqc_plots[target]['config']['sortHighlights'] == true){
      mqc_plots[target]['config']['sortHighlights'] = false;
      $(this).removeClass('active');
    } else {
      mqc_plots[target]['config']['sortHighlights'] = true;
      $(this).addClass('active');
    }
    $(this).blur();
    plot_heatmap(target);
  });

});

// Call to render any plot
function plot_graph(target, ds, max_num){
  if(mqc_plots[target] === undefined){ return false; }
  else {
    // XY Line charts
    if(mqc_plots[target]['plot_type'] == 'xy_line'){
      if(max_num === undefined || mqc_plots[target]['datasets'][0].length < max_num){
        plot_xy_line_graph(target, ds);
        $('#'+target).removeClass('not_rendered');
      } else {
        $('#'+target).addClass('not_rendered gt_max_num_ds').html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
      }
    }
    // Bar graphs
    else if(mqc_plots[target]['plot_type'] == 'bar_graph'){
      if(max_num === undefined || mqc_plots[target]['samples'][0].length < max_num){
        plot_stacked_bar_graph(target, ds);
        $('#'+target).removeClass('not_rendered');
      } else {
        $('#'+target).addClass('not_rendered gt_max_num_ds').html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
      }
    }
    // Scatter plots
    else if(mqc_plots[target]['plot_type'] == 'scatter'){
      if(max_num === undefined || Object.keys(mqc_plots[target]['datasets'][0]).length < max_num){
        plot_scatter_plot(target, ds);
        $('#'+target).removeClass('not_rendered');
      } else {
        $('#'+target).addClass('not_rendered gt_max_num_ds').html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
      }
    }
    // Beeswarm graphs
    else if(mqc_plots[target]['plot_type'] == 'beeswarm'){
      if(max_num === undefined || mqc_plots[target]['samples'][0].length < max_num){
        plot_beeswarm_graph(target, ds);
        $('#'+target).removeClass('not_rendered');
      } else {
        $('#'+target).addClass('not_rendered gt_max_num_ds').html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
      }
    }
    // Heatmap plots
    else if(mqc_plots[target]['plot_type'] == 'heatmap'){
      if(max_num === undefined || mqc_plots[target]['xcats'][0].length < max_num){
        plot_heatmap(target, ds);
        $('#'+target).removeClass('not_rendered');
      } else {
        $('#'+target).addClass('not_rendered gt_max_num_ds').html('<button class="btn btn-default btn-lg render_plot">Show plot</button>');
      }
    }
    // Not recognised
    else { console.log('Did not recognise plot type: '+mqc_plots[target]['plot_type']); }
  }
}

// Basic Line Graph
function plot_xy_line_graph(target, ds){
  if(mqc_plots[target] === undefined || mqc_plots[target]['plot_type'] !== 'xy_line'){
    return false;
  }
  var config = mqc_plots[target]['config'];
  var data = mqc_plots[target]['datasets'];
  if(ds === undefined){ ds = 0; }

  if(config['tt_label'] === undefined){ config['tt_label'] = '{point.x}: {point.y:.2f}'; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  else {
    config['click_func'] = eval("("+config['click_func']+")");
    if(config['cursor'] === undefined){ config['cursor'] = 'pointer'; }
  }
  if (config['xDecimals'] === undefined){ config['xDecimals'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  if (config['pointFormat'] === undefined){
    config['pointFormat'] = '<div style="background-color:{series.color}; display:inline-block; height: 10px; width: 10px; border:1px solid #333;"></div> <span style="text-decoration:underline; font-weight:bold;">{series.name}</span><br>'+config['tt_label'];
  }

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data in tact
  var data = JSON.parse(JSON.stringify(mqc_plots[target]['datasets'][ds]));

  // Rename samples
  if(window.mqc_rename_f_texts.length > 0){
    $.each(data, function(j, s){
      $.each(window.mqc_rename_f_texts, function(idx, f_text){
        if(window.mqc_rename_regex_mode){
          var re = new RegExp(f_text,"g");
          data[j]['name'] = data[j]['name'].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          data[j]['name'] = data[j]['name'].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    });
  }

  // Highlight samples
  if(window.mqc_highlight_f_texts.length > 0){
    $.each(data, function(j, s){
      $.each(window.mqc_highlight_f_texts, function(idx, f_text){
        if((window.mqc_highlight_regex_mode && data[j]['name'].match(f_text)) || (!window.mqc_highlight_regex_mode && data[j]['name'].indexOf(f_text) > -1)){
          data[j]['color'] = window.mqc_highlight_f_cols[idx];
        }
      });
    });
  }

  // Hide samples
  $('#'+target).closest('.mqc_hcplot_plotgroup').parent().find('.samples-hidden-warning').remove();
  $('#'+target).closest('.mqc_hcplot_plotgroup').show();
  if(window.mqc_hide_f_texts.length > 0){
    var num_hidden = 0;
    var num_total = data.length;
    var j = data.length;
    while (j--) {
      var match = false;
      for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
        var f_text = window.mqc_hide_f_texts[i];
        if(window.mqc_hide_regex_mode){
          if(data[j]['name'].match(f_text)){ match = true; }
        } else {
          if(data[j]['name'].indexOf(f_text) > -1){ match = true; }
        }
      }
      if(window.mqc_hide_mode == 'show'){
        match = !match;
      }
      if(match){
        data.splice(j,1);
        num_hidden += 1;
      }
    };
    // Some series hidden. Show a warning text string.
    if(num_hidden > 0) {
      var alert = '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+num_hidden+' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      $('#'+target).closest('.mqc_hcplot_plotgroup').before(alert);
    }
    // All series hidden. Hide the graph.
    if(num_hidden == num_total){
      $('#'+target).closest('.mqc_hcplot_plotgroup').hide();
      return false;
    }
  }

  // Toggle buttons for y-axis limis
  // Handler for this is at top, so doesn't get created multiple times
  if(config['ymax'] != undefined || config['ymin'] != undefined ){
    var pgroup = $('#'+target).closest('.mqc_hcplot_plotgroup');
    var wrapper = $('<div class="mqc_hcplot_yaxis_limit_toggle hidden-xs" />').prependTo(pgroup);
    wrapper.append('<span class="mqc_switch_wrapper" data-ymax="'+config['ymax']+'" data-ymin="'+config['ymin']+'" data-target="#'+target+'">Y-Limits: <span class="mqc_switch on">on</span></span>');
    wrapper.after('<div class="clearfix" />');
  }

  // Make the highcharts plot
  Highcharts.chart(target, {
    chart: {
      type: 'line',
      zoomType: 'x'
    },
    title: {
      text: config['title'],
      x: 30 // fudge to center over plot area rather than whole plot
    },
    xAxis: {
      title: {
        text: config['xlab']
      },
      labels: { format: config['xLabelFormat'] ? config['xLabelFormat']  : '{value}' },
      type: config['xLog'] ? 'logarithmic' : 'linear',
      categories: config['categories'],
      ceiling: config['xCeiling'],
      floor: config['xFloor'],
      max: config['xmax'],
      min: config['xmin'],
      minRange: config['xMinRange'],
      allowDecimals: config['xDecimals'],
      plotBands: config['xPlotBands'],
      plotLines: config['xPlotLines']
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      labels: { format: config['yLabelFormat'] ? config['yLabelFormat'] : '{value}' },
      type: config['yLog'] ? 'logarithmic' : 'linear',
      ceiling: config['yCeiling'],
      floor: config['yFloor'],
      max: config['ymax'],
      min: config['ymin'],
      minRange: config['yMinRange'],
      allowDecimals: config['yDecimals'],
      plotBands: config['yPlotBands'],
      plotLines: config['yPlotLines']
    },
    plotOptions: {
      series: {
        marker: { enabled: false },
        cursor: config['cursor'],
        point: {
          events: {
            click: config['click_func']
          }
        }
      }
    },
    legend: {
      enabled: false
    },
    tooltip: {
      headerFormat: '',
			pointFormat: config['pointFormat'],
			useHTML: true
    },
    series: data
  });
}

// Stacked Bar Graph
function plot_stacked_bar_graph(target, ds){
  if(mqc_plots[target] === undefined || mqc_plots[target]['plot_type'] !== 'bar_graph'){
    return false;
  }
  if(ds === undefined){ ds = 0; }

  // Make a clone of the everything, so that we can mess with it,
  // while keeping the original data in tact
  var data = JSON.parse(JSON.stringify(mqc_plots[target]['datasets'][ds]));
  var cats = JSON.parse(JSON.stringify(mqc_plots[target]['samples'][ds]));
  var config = JSON.parse(JSON.stringify(mqc_plots[target]['config']));

  if (config['stacking'] === undefined){ config['stacking'] = 'normal'; }
  if (config['ytype'] === undefined){ config['ytype'] = 'linear'; }
  if (config['reversedStacks'] === undefined){ config['reversedStacks'] = false; }
  if (config['use_legend'] === undefined){ config['use_legend'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  else { if(config['cursor'] === undefined){ config['cursor'] = 'pointer'; } }
  if (config['tt_percentages'] === undefined){ config['tt_percentages'] = true; }
  if (config['borderWidth'] === undefined){ config['borderWidth'] = 0; }

  if (config['ytype'] == 'logarithmic'){
    if(config['ymin'] == 0 || config['ymin'] == undefined){
      config['ymin'] = 1;
    }
    var minTickInt = 'auto';
  } else {
    var minTickInt = undefined;
  }

  // Rename samples
  if(window.mqc_rename_f_texts.length > 0){
    $.each(cats, function(j, s_name){
      $.each(window.mqc_rename_f_texts, function(idx, f_text){
        if(window.mqc_rename_regex_mode){
          var re = new RegExp(f_text,"g");
          cats[j] = cats[j].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          cats[j] = cats[j].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    });
  }

  // Highlight samples
  if(window.mqc_highlight_f_texts.length > 0){
    $.each(cats, function(j, s_name){
      $.each(window.mqc_highlight_f_texts, function(idx, f_text){
        if(f_text == ''){ return true; } // skip blanks
        if((window.mqc_highlight_regex_mode && s_name.match(f_text)) || (!window.mqc_highlight_regex_mode && s_name.indexOf(f_text) > -1)){
          // Make the data point in each series with this index have a border colour
          $.each(data, function(k, d){
            data[k]['data'][j] = {
              'y': data[k]['data'][j],
              'borderColor': window.mqc_highlight_f_cols[idx]
            }
          });
        }
      });
    });
    // Bump the borderWidth to make the highlights more obvious
    if(config['borderWidth'] <= 2){ config['borderWidth'] = 2; }
  }

  // Hide samples
  $('#'+target).closest('.mqc_hcplot_plotgroup').parent().find('.samples-hidden-warning').remove();
  $('#'+target).closest('.mqc_hcplot_plotgroup').show();
  if(window.mqc_hide_f_texts.length > 0){
    var num_hidden = 0;
    var num_total = cats.length;
    var j = cats.length;
    while (j--) {
      var s_name = cats[j];
      var match = false;
      for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
        var f_text = window.mqc_hide_f_texts[i];
        if(window.mqc_hide_regex_mode){
          if(s_name.match(f_text)){ match = true; }
        } else {
          if(s_name.indexOf(f_text) > -1){ match = true; }
        }
      }
      if(window.mqc_hide_mode == 'show'){
        match = !match;
      }
      if(match){
        cats.splice(j, 1);
        $.each(data, function(k, d){
          data[k]['data'].splice(j, 1);
        });
        num_hidden += 1;
      }
    };
    // Some series hidden. Show a warning text string.
    if(num_hidden > 0) {
      var alert = '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+num_hidden+' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      $('#'+target).closest('.mqc_hcplot_plotgroup').before(alert);
    }
    // All series hidden. Hide the graph.
    if(num_hidden == num_total){
      $('#'+target).closest('.mqc_hcplot_plotgroup').hide();
      return false;
    }
  }

  // Make the highcharts plot
  Highcharts.chart(target, {
    chart: {
      type: 'bar',
      zoomType: 'x'
    },
    title: {
      text: config['title'],
    },
    xAxis: {
      categories: cats,
      min: 0,
      title: {
        text: config['xlab']
      },
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      labels: { format: config['yLabelFormat'] ? config['yLabelFormat'] : '{value}' },
      ceiling: config['yCeiling'],
      floor: config['yFloor'],
      minRange: config['yMinRange'],
      max: config['ymax'],
      min: config['ymin'],
      type: config['ytype'],
      labels: {
        format: config['ylab_format']
      },
      allowDecimals: config['yDecimals'],
      reversedStacks: config['reversedStacks'],
      minorTickInterval: minTickInt
    },
    plotOptions: {
      series: {
        stacking: config['stacking'],
        groupPadding: 0.02,
        borderWidth: config['borderWidth']
      },
      cursor: config['cursor'],
      point: {
        events: {
          click: config['click_func']
        }
      }
    },
    legend: {
      enabled: config['use_legend']
    },
    tooltip: {
      formatter: function () {
        var colspan = config['tt_percentages'] ? 3 : 2;
        var s = '<table><tr><th colspan="'+colspan+'" style="font-weight:bold; text-decoration:underline;">' + this.x + '</th></tr>';
        $.each(this.points, function () {
          yval = Highcharts.numberFormat(this.y, (config['tt_decimals'] == undefined ? 0 : config['tt_decimals'])) + ( config['tt_suffix'] || '');
          ypct = Highcharts.numberFormat(this.percentage, 1);
          s += '<tr> \
            <td style="font-weight:bold; color:'+this.series.color+'; border-bottom:1px solid #dedede;">' + this.series.name + ':</td>\
            <td style="text-align:right; border-bottom:1px solid #dedede; padding: 0 15px;">' + yval + '</td>';
          if(config['tt_percentages']){
            s += '<td style="text-align:right; border-bottom:1px solid #dedede;">(' + ypct + '%)</td>';
          }
          s += '</tr>';
        });
        s += '</table>';
        return s;
      },
      shared: true,
      useHTML: true
    },
    series: data
  });
}


// Scatter plot
function plot_scatter_plot (target, ds){
  if(mqc_plots[target] === undefined || mqc_plots[target]['plot_type'] !== 'scatter'){
    return false;
  }
  var config = mqc_plots[target]['config'];
  var data = mqc_plots[target]['datasets'];
  if(ds === undefined){ ds = 0; }

  if(config['marker_colour'] === undefined){ config['marker_colour'] = 'rgba(124, 181, 236, .5)'; }
  if(config['marker_size'] === undefined){ config['marker_size'] = 5; }
  if(config['marker_line_colour'] === undefined){ config['marker_line_colour'] = '#999'; }
  if(config['marker_line_width'] === undefined){ config['marker_line_width'] = 1; }
  if(config['tt_label'] === undefined){ config['tt_label'] = 'X: <strong>{point.x:.2f}</strong><br/>Y: <strong>{point.y:.2f}</strong>'; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){ }; }
  else {
    config['click_func'] = eval("("+config['click_func']+")");
    if(config['cursor'] === undefined){ config['cursor'] = 'pointer'; }
  }
  if (config['xDecimals'] === undefined){ config['xDecimals'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  if (config['pointFormat'] === undefined){
    config['pointFormat'] = '<div style="background-color:{point.color}; display:inline-block; height: 10px; width: 10px; border:1px solid #333;"></div> <span style="text-decoration:underline; font-weight:bold;">{point.name}</span><br>'+config['tt_label'];
  }

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data in tact
  var data = JSON.parse(JSON.stringify(mqc_plots[target]['datasets'][ds]));

  // Rename samples
  if(window.mqc_rename_f_texts.length > 0){
    $.each(data, function(j, s){
      $.each(window.mqc_rename_f_texts, function(idx, f_text){
        if(window.mqc_rename_regex_mode){
          var re = new RegExp(f_text,"g");
          data[j]['name'] = data[j]['name'].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          data[j]['name'] = data[j]['name'].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    });
  }

  // Highlight samples
  if(window.mqc_highlight_f_texts.length > 0){
    $.each(data, function(j, s){
      if ('marker' in data[j]){
        data[j]['marker']['lineWidth'] = 0;
      } else {
        data[j]['marker'] = {'lineWidth': 0};
      }
      var match = false;
      $.each(window.mqc_highlight_f_texts, function(idx, f_text){
        if(f_text == ''){ return true; }
        if((window.mqc_highlight_regex_mode && data[j]['name'].match(f_text)) || (!window.mqc_highlight_regex_mode && data[j]['name'].indexOf(f_text) > -1)){
          data[j]['color'] = window.mqc_highlight_f_cols[idx];
          match = true;
        }
      });
      if(!match) {
        data[j]['color'] = 'rgba(100,100,100,0.2)';
      }
    });
  }

  // Hide samples
  $('#'+target).closest('.mqc_hcplot_plotgroup').parent().find('.samples-hidden-warning').remove();
  $('#'+target).closest('.mqc_hcplot_plotgroup').show();
  if(window.mqc_hide_f_texts.length > 0){
    var num_hidden = 0;
    var num_total = data.length;
    var j = data.length;
    while (j--) {
      var match = false;
      for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
        var f_text = window.mqc_hide_f_texts[i];
        if(window.mqc_hide_regex_mode){
          if(data[j]['name'].match(f_text)){ match = true; }
        } else {
          if(data[j]['name'].indexOf(f_text) > -1){ match = true; }
        }
      }
      if(window.mqc_hide_mode == 'show'){
        match = !match;
      }
      if(match){
        data.splice(j,1);
        num_hidden += 1;
      }
    };
    // Some series hidden. Show a warning text string.
    if(num_hidden > 0) {
      var alert = '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+num_hidden+' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      $('#'+target).closest('.mqc_hcplot_plotgroup').before(alert);
    }
    // All series hidden. Hide the graph.
    if(num_hidden == num_total){
      $('#'+target).closest('.mqc_hcplot_plotgroup').hide();
      return false;
    }
  }

  // Make the highcharts plot
  Highcharts.chart(target, {
    chart: {
      type: 'scatter',
      zoomType: 'xy',
      plotBorderWidth: 1,
      height: config['square'] ? 500 : undefined,
      width: config['square'] ? 500 : undefined
    },
    title: {
      text: config['title'],
      x: 30 // fudge to center over plot area rather than whole plot
    },
    xAxis: {
      title: {
        text: config['xlab']
      },
      type: config['xLog'] ? 'logarithmic' : 'linear',
      gridLineWidth: 1,
      categories: config['categories'],
      ceiling: config['xCeiling'],
      floor: config['xFloor'],
      max: config['xmax'],
      min: config['xmin'],
      minRange: config['xMinRange'],
      allowDecimals: config['xDecimals'],
      plotBands: config['xPlotBands'],
      plotLines: config['xPlotLines']
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      type: config['yLog'] ? 'logarithmic' : 'linear',
      ceiling: config['yCeiling'],
      floor: config['yFloor'],
      max: config['ymax'],
      min: config['ymin'],
      minRange: config['yMinRange'],
      allowDecimals: config['yDecimals'],
      plotBands: config['yPlotBands'],
      plotLines: config['yPlotLines']
    },
    plotOptions: {
      series: {
        animation: false,
        marker: {
          radius: config['marker_size'],
          lineColor: config['marker_line_colour'],
          lineWidth: config['marker_line_width'],
          states: {
            hover: {
              enabled: config['enableHover'] == undefined ? true : config['enableHover'],
              lineColor: 'rgb(100,100,100)'
            }
          }
        },
        turboThreshold: config['turboThreshold'],
        enableMouseTracking: config['enableMouseTracking'],
        cursor: config['cursor'],
        point: {
          events: {
            click: config['click_func']
          }
        }
      }
    },
    legend: {
      enabled: false
    },
    tooltip: {
      headerFormat: '',
			pointFormat: config['pointFormat'],
			useHTML: true,
      formatter: (function() {
        if(!this.point.noTooltip) {
          // Formatter function doesn't do name for some reason
          fstring = config['pointFormat'].replace('{point.name}', this.point.name);
          return Highcharts.Point.prototype.tooltipFormatter.call(this, fstring);
        }
        return false;
      })
    },
    series: [{
      color: config['marker_colour'],
      data: data
    }]
  },
  // Maintain aspect ratio as chart size changes
  function(this_chart){
    if(config['square']){
      var resizeCh = function(chart){
        // Extra width for legend
        var lWidth = chart.options.legend.enabled ? 30 : 0;
        // Work out new chart width, assuming needs to be narrower
        var chHeight = $(chart.renderTo).height();
        var chWidth = $(chart.renderTo).width();
        var nChHeight = chHeight;
        var nChWidth = chHeight + lWidth;
        // Chart is already too narrow, make it less tall
        if(chWidth < nChWidth){
          nChHeight = chWidth - lWidth;
          nChWidth = chWidth;
        }
        chart.setSize(nChWidth, nChHeight);
      }
      // Resize on load
      resizeCh(this_chart);
      // Resize on graph resize
      $(this_chart.renderTo).on('mqc_plotresize', function(e){
        resizeCh(this_chart);
      });
    }
  });
}

// Beeswarm plot
function plot_beeswarm_graph(target, ds){
  if(mqc_plots[target] === undefined || mqc_plots[target]['plot_type'] !== 'beeswarm'){
    return false;
  }
  var config = mqc_plots[target]['config'];
  if(ds === undefined){ ds = 0; }

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data in tact
  var datasets = JSON.parse(JSON.stringify(mqc_plots[target]['datasets']));
  var samples = JSON.parse(JSON.stringify(mqc_plots[target]['samples']));
  var categories = JSON.parse(JSON.stringify(mqc_plots[target]['categories']));

  // Rename samples
  if(window.mqc_rename_f_texts.length > 0){
    for (i=0; i < samples.length; i++) {
      for (j=0; j < samples[i].length; j++) {
        $.each(window.mqc_rename_f_texts, function(idx, f_text){
          if(window.mqc_rename_regex_mode){
            var re = new RegExp(f_text,"g");
            samples[i][j] = samples[i][j].replace(re, window.mqc_rename_t_texts[idx]);
          } else {
            samples[i][j] = samples[i][j].replace(f_text, window.mqc_rename_t_texts[idx]);
          }
        });
      }
    }
  }

  // Highlight samples
  var baseColour = 'rgb(55,126,184)'; // Blue points by default
  var seriesColours = {};
  if(window.mqc_highlight_f_texts.length > 0){
    baseColour = 'rgb(80,80,80)'; // Grey points if no highlight
    for (i=0; i < samples.length; i++) {
      for (j=0; j < samples[i].length; j++) {
        $.each(window.mqc_highlight_f_texts, function(idx, f_text){
          if((window.mqc_highlight_regex_mode && samples[i][j].match(f_text)) || (!window.mqc_highlight_regex_mode && samples[i][j].indexOf(f_text) > -1)){
            seriesColours[samples[i][j]] = window.mqc_highlight_f_cols[idx];
          }
        });
      }
    }
  }

  // Hide samples
  $('#'+target).closest('.hc-plot-wrapper').parent().find('.samples-hidden-warning').remove();
  $('#'+target).closest('.hc-plot-wrapper').show();
  if(window.mqc_hide_f_texts.length > 0){
    var num_hidden = 0;
    var num_total = 0;
    for (i=0; i < samples.length; i++) {
      num_total = Math.max(num_total, samples[i].length);
      var j = samples[i].length;
      var hidden_here = 0;
      while (j--) {
        var s_name = samples[i][j];
        var match = false;
        for (k = 0; k < window.mqc_hide_f_texts.length; k++) {
          var f_text = window.mqc_hide_f_texts[k];
          if(window.mqc_hide_regex_mode){
            if(s_name.match(f_text)){ match = true; }
          } else {
            if(s_name.indexOf(f_text) > -1){ match = true; }
          }
        }
        if(window.mqc_hide_mode == 'show'){
          match = !match;
        }
        if(match){
          samples[i].splice(j, 1);
          datasets[i].splice(j, 1);
          hidden_here += 1;
        }
      };
      num_hidden = Math.max(num_hidden, hidden_here);
    };
    // Some series hidden. Show a warning text string.
    if(num_hidden > 0) {
      var alert = '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+num_hidden+' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      $('#'+target).closest('.hc-plot-wrapper').before(alert);
    }
    // All series hidden. Hide the graph.
    if(num_hidden == num_total){
      $('#'+target).closest('.hc-plot-wrapper').hide();
      return false;
    }
  }

  // Figure out how tall to make each plot
  var ph_min = 40;
  var ph_max = 100;
  var pheight = 600 / categories.length;
  pheight = Math.min(ph_max, Math.max(ph_min, pheight));

  // Clear the loading text and add hover text placeholder
  $('#'+target).html('<div class="beeswarm-hovertext"><em class="placeholder">Hover over a data point for more information</em></div><div class="beeswarm-plots"></div>');
  // Resize the parent draggable div
  $('#'+target).parent().css('height', ((pheight*categories.length)+40)+'px');

  for (var i = 0; i < categories.length; i++) {

    var borderCol = categories[i]['bordercol'];
    if (borderCol == undefined){
      borderCol = '#cccccc';
    }

    var data = datasets[i];
    var s_names = samples[i];
    if (categories[i]['namespace'] == ''){
      var label = categories[i]['title'];
      var label_long = categories[i]['description'];
    } else{
      var label = categories[i]['namespace'] + '<br/>' + categories[i]['title'];
      var label_long = categories[i]['namespace'] + ': ' + categories[i]['description'];
    }
    var ttSuffix = categories[i]['suffix'];
    var decimalPlaces = categories[i]['decimalPlaces'];
    var minx = categories[i]['min'];
    var maxx = categories[i]['max'];

    // Size and spacing options
    var markerRadius = 2.5
  	var yspace = 70;
    var ysep = 10;
    if(data.length > 50){
      markerRadius = 1.8
      yspace = 50;
      ysep = 20;
    }
    if(data.length > 200){
      markerRadius = 1
      yspace = 30;
      ysep = 30;
    }

    if (maxx == undefined){
    	maxx = Math.max.apply(null, data);
    }
    if (minx == undefined){
    	minx = Math.max.apply(null, data);
    }
    var range = maxx-minx;
    var sep = range/yspace;
    // Get an array of indexes from a sorted data array
    // Leaves the data order in tact so we don't lose s_name association
    var indices = new Array(data.length);
    for (var n = 0; n < data.length; n++) { indices[n] = n; }
    indices.sort(function (a, b) {
      return data[a] < data[b] ? -1 : data[a] > data[b] ? 1 : 0;
    });
    var xydata = [];
    var last = undefined;
    var side = 1;
    for (var s_idx = 0; s_idx < indices.length; s_idx++) {
      row = indices[s_idx];
      s_name = s_names[row];
      d = data[row];
      if (Math.floor(d/sep) !== last){
        last = Math.floor(d/sep);
        side = 1;
      } else {
        side += 1;
      }
      multiplier = (side % 2 == 0) ? 1 : -1;
      var y = (Math.floor(side/2) * multiplier)/ysep;
      // Don't let jitter get too big
      while(y > 1 || y < -1){
        var n = Math.floor(Math.abs(y)) + 1;
        y = (Math.floor(side/2) * multiplier)/(ysep*n);
      }
      // Get the point colour
      var thisCol = baseColour;
      if(s_name in seriesColours) {
        thisCol = seriesColours[s_name];
      }
      xydata.push({
        'x':d,
        'y':y,
        'name':s_name,
        'color': thisCol
      });
    }

    $('<div class="beeswarm-plot" />')
      .appendTo('#'+target+' .beeswarm-plots')
      .css({
        'border-left': '2px solid '+borderCol,
        'height': (100/categories.length)+'%'
      })
      .highcharts({
  			chart: {
            type: 'scatter',
            spacingTop: 0,
            marginBottom: 0,
            marginRight: 20,
            marginLeft: 180,
            backgroundColor: 'transparent',
            // Horrible hacky HighCharts reflow problem.
            // TODO: Come back and find a better solution!
            events: {
              load: function(chart) {
                setTimeout(function(){
                  chart.target.reflow();
                }, 200);
              }
            }
        },
        title: {
          text: label,
          align: 'left',
          verticalAlign: 'middle',
          y: 10,
          useHTML: true,
          style: {
              fontSize: '12px'
          }
        },
        yAxis: {
            title: {text: null},
            max: 1,
            min: -1,
            gridLineWidth: 0,
            title: {text: null},
            labels: {enabled: false},
            lineWidth: 0
        },
        xAxis: {
        	lineWidth: 0,
          tickWidth: 0,
          tickPixelInterval: 200,
          labels: {
            reserveSpace: false,
            y: (-1*(pheight/2))+5,
            zIndex: 1,
            style: {
                color: '#999999'
            }
          },
          min: minx,
          max: maxx,
        },
        tooltip: {
          valueSuffix: ttSuffix,
          valueDecimals: decimalPlaces,
          formatter: function(){
            var value = Highcharts.numberFormat(this.point.x, this.series.tooltipOptions.valueDecimals);
            var suff = this.series.tooltipOptions.valueSuffix;
            var ttstring = '<span style="float:right;">'+this.series.name+'</span><samp>'+this.point.name+'</samp>: &nbsp; <strong>'+value+' '+suff+'</strong>';
            $('#'+target+' .beeswarm-hovertext').html(ttstring);
            return false;
          }
        },
        plotOptions: {
          series: {
            name: label_long,
            turboThreshold: 0,
            marker: {
              radius: markerRadius,
              states: {
                hover: {
                  radiusPlus: 4,
                  lineWidthPlus: 2,
                  lineColor: '#333333'
                }
              }
            },
            stickyTracking: false,
            point: {
              events: {
                mouseOver: function (e) {
                  var hovName = this.name;
                  $('#'+target+' .beeswarm-plot').each(function(){
                    var plot = $(this).highcharts();
                    for (i = 0; i < plot.series[0].data.length; ++i) {
                      if(plot.series[0].data[i].name == hovName){
                        plot.series[0].data[i].setState('hover');
                      }
                    }
                  });

                },
                mouseOut: function () {
                  $('#'+target+' .beeswarm-plot').each(function(){
                    var plot = $(this).highcharts();
                    for (i = 0; i < plot.series[0].data.length; ++i) {
                      plot.series[0].data[i].setState();
                    }
                  });
                  $('#'+target+' .beeswarm-hovertext').html('<em class="placeholder">Hover over a data point for more information</em>');
                }
              }
            }
          }
        },
        legend: { enabled: false },
        credits: { enabled: false },
        exporting: { enabled: false },
        series: [{
          data: xydata,
          // Workaround for HighCharts bug. See https://github.com/highcharts/highcharts/issues/1440
          marker: { states: { hover: { fillColor: {} } } }
        }]

    });

  }
}

// Heatmap plot
function plot_heatmap(target, ds){
  if(mqc_plots[target] === undefined || mqc_plots[target]['plot_type'] !== 'heatmap'){
    return false;
  }
  var config = mqc_plots[target]['config'];

  if(config['square'] === undefined){ config['square'] = true; }

  // Make a clone of the data, so that we can mess with it,
  // while keeping the original data in tact
  var data = JSON.parse(JSON.stringify(mqc_plots[target]['data']));
  var xcats = JSON.parse(JSON.stringify(mqc_plots[target]['xcats']));
  var ycats = JSON.parse(JSON.stringify(mqc_plots[target]['ycats']));

  // Rename samples
  if(window.mqc_rename_f_texts.length > 0){
    for (i=0; i < xcats.length; i++) {
      $.each(window.mqc_rename_f_texts, function(idx, f_text){
        if(window.mqc_rename_regex_mode){
          var re = new RegExp(f_text,"g");
          xcats[i] = xcats[i].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          xcats[i] = xcats[i].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    }
    for (i=0; i < ycats.length; i++) {
      $.each(window.mqc_rename_f_texts, function(idx, f_text){
        if(window.mqc_rename_regex_mode){
          var re = new RegExp(f_text,"g");
          ycats[i] = ycats[i].replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          ycats[i] = ycats[i].replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
    }
  }

  // Sort samples by highlight
  $('.mqc_heatmap_sortHighlight').attr('disabled', false);
  if(config['sortHighlights'] == true){
    if(window.mqc_highlight_f_texts.length > 0){
      // Collect the highlighting indices
      var xcat_hl = Array();
      var ycat_hl = Array();
      for (i=0; i < xcats.length; i++) {
        $.each(window.mqc_highlight_f_texts, function(idx, f_text){
          if(f_text == ''){ xcat_hl[i] = 0; }
          else if((window.mqc_highlight_regex_mode && xcats[i].match(f_text)) || (!window.mqc_highlight_regex_mode && xcats[i].indexOf(f_text) > -1)){
            xcat_hl[i] = window.mqc_highlight_f_texts.length - idx;
          }
        });
      }
      for (i=0; i < ycats.length; i++) {
        $.each(window.mqc_highlight_f_texts, function(idx, f_text){
          if(f_text == ''){ ycat_hl[i] = 0; }
          else if((window.mqc_highlight_regex_mode && ycats[i].match(f_text)) || (!window.mqc_highlight_regex_mode && ycats[i].indexOf(f_text) > -1)){
            ycat_hl[i] = window.mqc_highlight_f_texts.length - idx;
          }
        });
      }
      // Reshape the data - needs deepcopy as indexes are updated
      var newdata = JSON.parse(JSON.stringify(mqc_plots[target]['data']));
      var new_xcats = [], new_ycats = [];
      var xidx = 0, yidx = 0;
      for (hl = window.mqc_highlight_f_texts.length; hl >= 0; hl--){
        for (i=0; i < xcats.length; i++) {
          if(xcat_hl[i] == hl){
            new_xcats.push(xcats[i])
            for (j=0; j < data.length; j++) {
              if(data[j][0] == i){ newdata[j][0] = xidx; }
            }
            xidx += 1;
          }
        }
        for (i=0; i < ycats.length; i++) {
          if(ycat_hl[i] == hl){
            new_ycats.push(ycats[i])
            for (j=0; j < data.length; j++) {
              if(data[j][1] == i){ newdata[j][1] = yidx; }
            }
            yidx += 1;
          }
        }
      }
      data = newdata;
      xcats = new_xcats;
      ycats = new_ycats;
    }
  }

  // Hide samples
  var num_total = Math.max(xcats.length, ycats.length);
  $('#'+target).closest('.hc-plot-wrapper').parent().find('.samples-hidden-warning').remove();
  $('#'+target).closest('.hc-plot-wrapper').show();
  if(window.mqc_hide_f_texts.length > 0){
    var remove = Array();
    var i = xcats.length;
    var xhidden = 0;
    while (i--) {
      var match = false;
      for (j = 0; j < window.mqc_hide_f_texts.length; j++) {
        var f_text = window.mqc_hide_f_texts[j];
        if(window.mqc_hide_regex_mode){
          if(xcats[i].match(f_text)){ match = true; }
        } else {
          if(xcats[i].indexOf(f_text) > -1){ match = true; }
        }
      }
      if(window.mqc_hide_mode == 'show'){
        match = !match;
      }
      if(match){
        xcats.splice(i, 1);
        for (n=0; n < data.length; n++) {
          var x = data[n][1];
          if (x == i){ remove.push(n); }
          else if(x > i){ data[n][1] -= 1; }
        }
        xhidden += 1;
      }
    }
    var i = ycats.length;
    var yhidden = 0;
    while (i--) {
      var match = false;
      for (j = 0; j < window.mqc_hide_f_texts.length; j++) {
        var f_text = window.mqc_hide_f_texts[j];
        if(window.mqc_hide_regex_mode){
          if(ycats[i].match(f_text)){ match = true; }
        } else {
          if(ycats[i].indexOf(f_text) > -1){ match = true; }
        }
      }
      if(window.mqc_hide_mode == 'show'){
        match = !match;
      }
      if(match){
        ycats.splice(i, 1);
        for (n=0; n < data.length; n++) {
          var y = data[n][0];
          if (y == i){
            if(remove.indexOf(n) < 0){ remove.push(n); }
          } else if(y > i){ data[n][0] -= 1; }
        }
        yhidden += 1;
      }
    }
    // Remove the data values that matched
    remove = remove.sort(function(a, b){return a-b}); // Sorts alphabetically by default, even with integers
    var r = remove.length;
    while(r--){
      data.splice( remove[r], 1);
    }
    // Report / hide the plot if we're hiding stuff
    var num_hidden = Math.max(xhidden, yhidden);
    // Some series hidden. Show a warning text string.
    if(num_hidden > 0) {
      var alert = '<div class="samples-hidden-warning alert alert-warning"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+num_hidden+' samples hidden. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
      $('#'+target).closest('.hc-plot-wrapper').before(alert);
    }
    // All series hidden. Hide the graph.
    if(num_hidden >= num_total){
      $('#'+target).closest('.hc-plot-wrapper').hide();
      return false;
    }
  }

  // Highlight samples - do this last as we convert numerical arrays to associative
  if(window.mqc_highlight_f_texts.length > 0){
    $('.mqc_heatmap_sortHighlight').attr('disabled', false);
    var highlight_cells = Array();
    for (i=0; i < xcats.length; i++) {
      $.each(window.mqc_highlight_f_texts, function(idx, f_text){
        if(f_text == ''){ return true; }
        if((window.mqc_highlight_regex_mode && xcats[i].match(f_text)) || (!window.mqc_highlight_regex_mode && xcats[i].indexOf(f_text) > -1)){
          for (n=0; n < data.length; n++) {
            highlight_cells[idx] = ( typeof highlight_cells[idx] != 'undefined' && highlight_cells[idx] instanceof Array ) ? highlight_cells[idx] : [];
            if (data[n][1] == i){ highlight_cells[idx].push(n); }
          }
        }
      });
    }
    for (i=0; i < ycats.length; i++) {
      $.each(window.mqc_highlight_f_texts, function(idx, f_text){
        if(f_text == ''){ return true; }
        if((window.mqc_highlight_regex_mode && ycats[i].match(f_text)) || (!window.mqc_highlight_regex_mode && ycats[i].indexOf(f_text) > -1)){
          for (n=0; n < data.length; n++) {
            if (data[n][0] == i){
              highlight_cells[idx] = ( typeof highlight_cells[idx] != 'undefined' && highlight_cells[idx] instanceof Array ) ? highlight_cells[idx] : [];
              if(highlight_cells[idx].indexOf(n) < 0){ highlight_cells[idx].push(n); }
            }
          }
        }
      });
    }
    // Give highlighted cells a border
    for (var idx in highlight_cells){
      var hl = highlight_cells[idx];
      hl = hl.sort(function(a, b){return a-b}); // Sorts alphabetically by default, even with integers
      var h = hl.length;
      while(h--){
        var i = hl[h];
        data[i] = {
          x: data[i][1] === undefined ? data[i]['x'] : data[i][1],
          y: data[i][0] === undefined ? data[i]['y'] : data[i][0],
          value:data[i][2] === undefined ? data[i]['value'] : data[i][2],
          borderWidth:2,
          borderColor: window.mqc_highlight_f_cols[idx]
        }
      }
    }
  } else {
    $('.mqc_heatmap_sortHighlight').attr('disabled', true);
  }

  // We set undefined config vars so that they stay the same when hiding samples
  if(config['min'] === undefined || config['max'] === undefined){
    var dmin = data[0][2];
    var dmax = data[0][2];
    for (n=0; n < data.length; n++) {
      dmin = Math.min(dmin, data[n][2]);
      dmax = Math.max(dmax, data[n][2]);
    }
    if(config['min'] === undefined){ config['min'] = dmin; }
    if(config['max'] === undefined){ config['max'] = dmax; }
  }
  if(config['colstops'] === undefined){
    config['colstops'] = [
      [0, '#313695'],
      [0.1, '#4575b4'],
      [0.2, '#74add1'],
      [0.3, '#abd9e9'],
      [0.4, '#e0f3f8'],
      [0.5, '#ffffbf'],
      [0.6, '#fee090'],
      [0.7, '#fdae61'],
      [0.8, '#f46d43'],
      [0.9, '#d73027'],
      [1, '#a50026'],
    ];
  }
  if(config['reverseColors'] === undefined){ config['reverseColors'] = false; }
  if(config['decimalPlaces'] === undefined){ config['decimalPlaces'] = 2; }
  if(config['legend'] === undefined){ config['legend'] = true; }
  if(config['borderWidth'] === undefined){ config['borderWidth'] = 0; }
  var datalabels = config['datalabels'];
  if(datalabels === undefined){
    if(data.length < 20){ datalabels = true; }
    else { datalabels = false; }
  }
  // Clone the colstops before we mess around with them
  var colstops = JSON.parse(JSON.stringify(config['colstops']));
  // Reverse the colour scale if the axis is reversed
  if(config['reverseColors']){
    for(var i = 0; i < colstops.length; i++){
      colstops[i][0] = 1 - colstops[i][0];
    }
    colstops.reverse();
  }

  // Make the highcharts plot
  Highcharts.chart(target, {
    chart: {
      type: 'heatmap',
      zoomType: 'xy',
      height: config['square'] ? 500 : undefined,
      width: config['square'] ? 530 : undefined,
      marginTop: config['title'] ? 60 : 50
    },
    plotOptions: {
      series: {
        point: {
          events: {
            mouseOver: function() {
              // Stop highcharts making squares blue on hover
              this.pointAttr.hover.fill = this.color;
            }
          }
        },
        states: {
          hover: {
            borderWidth: 2,
            borderColor: 'red'
          }
        }
      }
    },
    title: {
      text: config['title'],
    },
    xAxis: {
      endOnTick: false,
      maxPadding: 0,
      categories: xcats,
      title: { enabled: true, text: config['xTitle'] },
      labels: {
        formatter: function(){
          try { return this.value.substr(0, 20); }
          catch(err) { return this.value; }
        }
      }
    },
    yAxis: {
      endOnTick: false,
      maxPadding: 0,
      categories: ycats,
      reversed: true,
      opposite: true,
      title: config['yTitle'],
      labels: {
        formatter: function(){
          try { return this.value.substr(0, 20); }
          catch(err) { return this.value; }
        }
      }
    },
    colorAxis: {
      reversed: config['reverseColors'],
      stops: colstops,
      min: config['min'],
      max: config['max'],
    },
    legend: {
      align: 'right',
      layout: 'vertical',
      margin: 0,
      verticalAlign: 'top',
      y: 25,
      symbolHeight: 280,
      enabled: config['legend']
    },
    tooltip: {
      useHTML: true,
      formatter: function () {
        return 'X: <span style="font-weight:bold; font-family:monospace;">'+this.series.xAxis.categories[this.point.x] + '</span><br>' +
        'Y: <span style="font-weight:bold; font-family:monospace;">' + this.series.yAxis.categories[this.point.y] + '</span><br>' +
        '<div style="background-color:'+this.point.color+'; display:inline-block; height: 10px; width: 10px; border:1px solid #333;"></div> ' +
        '<span style="font-weight: bold; text-decoration:underline;">' + Highcharts.numberFormat(this.point.value, config['decimalPlaces']) + '</span>'
      }
    },
    series: [{
      turboThreshold: 0,
      borderWidth: config['borderWidth'],
      data: data,
      dataLabels: {
        enabled: datalabels,
        format: '{point.value:.'+config['decimalPlaces']+'f}',
        color: config['datalabel_colour']
      }
    }]
  }, function(this_chart){
    // Maintain aspect ratio as chart size changes
    if(config['square']){
      var resizeCh = function(chart){
        // Extra width for legend
        var lWidth = chart.options.legend.enabled ? 30 : 0;
        // Work out new chart width, assuming needs to be narrower
        var chHeight = $(chart.renderTo).height();
        var chWidth = $(chart.renderTo).width();
        var nChHeight = chHeight;
        var nChWidth = chHeight + lWidth;
        // Chart is already too narrow, make it less tall
        if(chWidth < nChWidth){
          nChHeight = chWidth - lWidth;
          nChWidth = chWidth;
        }
        chart.setSize(nChWidth, nChHeight);
      }
      // Resize on load
      resizeCh(this_chart);
      // Resize on graph resize
      $(this_chart.renderTo).on('mqc_plotresize', function(e){
        try {
          resizeCh(this_chart);
        } catch(e){
          plot_heatmap($(this).attr('id'));
        }
      });
    }
  });

}

// Highlight text with a fadeout background colour highlight
function highlight_fade_text(obj){
  var orig_col = $(obj).css('color');
  obj.css({
    'display'          : 'inline-block',
    'background-color' : '#5bc0de',
    'color'            : '#FFFFFF',
    'WebkitTransition' : 'background-color 0s, color 0s',
    'MozTransition'    : 'background-color 0s, color 0s',
    'MsTransition'     : 'background-color 0s, color 0s',
    'OTransition'      : 'background-color 0s, color 0s',
    'transition'       : 'background-color 0s, color 0s'
  });
  setTimeout(function(){
    obj.css({
      'background-color' : '#FFFFFF',
      'color'            : orig_col,
      'WebkitTransition' : 'background-color 0.5s, color 0.5s',
      'MozTransition'    : 'background-color 0.5s, color 0.5s',
      'MsTransition'     : 'background-color 0.5s, color 0.5s',
      'OTransition'      : 'background-color 0.5s, color 0.5s',
      'transition'       : 'background-color 0.5s, color 0.5s'
    });
  }, 500);
}