////////////////////////////////////////////////
// HighCharts Plotting Code
////////////////////////////////////////////////

// Execute when page load has finished loading
$(function () {
    
    // Switch a HighCharts axis or data source
  $('.switch_group button').click(function(e){
    e.preventDefault();
    $(this).siblings('button.active').removeClass('active');
    $(this).addClass('active');
    var target = $(this).data('target');
    var action = $(this).data('action');
    var plot_options = highcharts_plot_options[target];
    // Switch between values and percentages
    if(action == 'set_percent' || action == 'set_numbers'){
      var sym = (action == 'set_percent') ? '%' : '#';
      var stack_type = (action == 'set_percent') ? 'percent' : 'normal';
      if(plot_options.yAxis.title.text !== undefined){
        plot_options.yAxis.title.text = sym+plot_options.yAxis.title.text.substr(1);
      }
      plot_options.plotOptions.series.stacking = stack_type;
      $(target).highcharts(plot_options);
    }
    // Switch data source
    if(action == 'set_data'){
      var ds = $(this).data('newdata');
      var data = eval(ds);
      if(data === undefined){
        console.log('Error switching plot dataset - '+ds+' not found..');
      } else {
        plot_options.series = data;
        var ylab = $(this).data('ylab');
        if(ylab !== undefined){ plot_options.yAxis.title.text = ylab; }
        var xlab = $(this).data('xlab');
        if(xlab !== undefined){ plot_options.xAxis.title.text = xlab; }
        $(target).highcharts(plot_options);
      }
    }
  });

  // Make HighCharts divs height-draggable
  // http://jsfiddle.net/Lkwb86c8/
  $('.hc-plot').each(function(){
    if(!$(this).parent().hasClass('hc-plot-wrapper')){
      $(this).wrap('<div class="hc-plot-wrapper"></div>');
    }
    if(!$(this).siblings().hasClass('hc-plot-handle')){
      $(this).after('<div class="hc-plot-handle"><span></span><span></span><span></span></div>');
    }
  });
  $('.hc-plot').css({ height: 'auto', top: 0, bottom: '10px', position: 'absolute' });
  $('.hc-plot-handle').on('mousedown', function(e){
    var wrapper = $(this).parent();
    var handle = $(this);
    var startHeight = wrapper.height();
    var pY = e.pageY;
    $(document).on('mouseup', function(e){
      $(document).off('mousemove');
      $(document).resize(); // send resize trigger for replotting
      // Fire off a custom jQuery event for other javascript chunks to tie into
      // Bind to the plot div, which should have a custom ID
      $(wrapper.parent().find('.hc-plot')).trigger('mqc_plotresize');
    });
    $(document).on('mousemove', function(me){
      wrapper.css('height', startHeight + (me.pageY - pY));
    });
  });
  
});



// We store the config options for every graph in an arrays.
// This way, we can go back and change them via button clicks
// eg. Changing an axis from values to percentages
highcharts_plots = [];
highcharts_plot_options = [];

// Basic Line Graph
function plot_xy_line_graph(div, data, config){
  if(config['tt_label'] === undefined){ config['tt_label'] = '{point.x}: {point.y:.2f}'; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  else { if(config['cursor'] === undefined){ config['cursor'] = 'pointer'; } }
  if (config['xDecimals'] === undefined){ config['xDecimals'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  if (config['pointFormat'] === undefined){
    config['pointFormat'] = '<div style="background-color:{series.color}; display:inline-block; height: 10px; width: 10px; border:1px solid #333;"></div> <span style="text-decoration:underline; font-weight:bold;">{series.name}</span><br>'+config['tt_label'];
  }
  // Collect the starting sample names to preserve after renaming
  var s_names = [];
  $.each(data, function(idx, d){
    s_names.push(d['name']);
  });
  // Make the highcharts config object
  highcharts_plot_options[div] = {
    chart: {
      renderTo: div.replace('#',''),
      type: 'line',
      zoomType: 'x',
      backgroundColor: null,
    },
    title: {
      text: config['title'],
      x: 30 // fudge to center over plot area rather than whole plot
    },
    xAxis: {
      title: {
        text: config['xlab']
      },
      categories: config['categories'],
      ceiling: config['xCeiling'],
      max: config['xmax'],
      min: config['xmin'],
      minTickInterval: config['xMinTickInterval'],
      allowDecimals: config['xDecimals'],
      plotBands: config['xPlotBands']
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      ceiling: config['yCeiling'],
      max: config['ymax'],
      min: config['ymin'],
      minTickInterval: config['yMinTickInterval'],
      allowDecimals: config['yDecimals'],
      plotBands: config['yPlotBands']
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
    credits: {
			enabled: true,
      text: 'Created with MultiQC',
      href: 'https://github.com/ewels/MultiQC'
		},
    tooltip: {
      headerFormat: '',
			pointFormat: config['pointFormat'],
			useHTML: true
    },
    series: data,
    s_names: s_names // Not for highcharts - this is to remember original names for renaming
  }
  highcharts_plots[div] = new Highcharts.Chart(highcharts_plot_options[div]);
}

// Stacked Bar Graph
function plot_stacked_bar_graph(div, cats, data, config){
  if (config['stacking'] === undefined){ config['stacking'] = 'normal'; }
  if (config['use_legend'] === undefined){ config['use_legend'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  else { if(config['cursor'] === undefined){ config['cursor'] = 'pointer'; } }
  // Make the highcharts config object
  highcharts_plot_options[div] = {
    chart: {
      renderTo: div.replace('#',''),
      type: 'bar',
      backgroundColor: null,
    },
    title: {
      text: config['title'],
    },
    xAxis: {
      categories: cats,
      min: 0,
      title: {
        text: config['xlab']
      }
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      max: config['ymax'],
      min: config['ymin'],
      allowDecimals: config['yDecimals'],
      reversedStacks: false
    },
    plotOptions: {
      series: {
        stacking: config['stacking'],
        groupPadding: 0.02
      },
      cursor: config['cursor'],
      point: {
        events: {
          click: config['click_func']
        }
      }
    },
    credits: {
			enabled: true,
      text: 'Created with MultiQC',
      href: 'https://github.com/ewels/MultiQC'
		},
    legend: {
      enabled: config['use_legend']
    },
    tooltip: {
      formatter: function () {
        var s = '<table><tr><th colspan="3" style="font-weight:bold; text-decoration:underline;">' + this.x + '</th></tr>';
        $.each(this.points.reverse(), function () {
          s += '<tr><td style="font-weight:bold; color:'+this.series.color+'; padding-right: 15px;">' + this.series.name + ':</td><td style="text-align:right;">' + numberWithCommas(this.y) + '</td><td style="text-align:right;"> (' + this.percentage.toFixed(1) + '%)</td></tr>';
        });
        s += '</table>';
        return s;
      },
      shared: true,
      useHTML: true
    },
    series: data,
    s_names: cats // Not for highcharts - this is to remember original names for renaming
  }
  highcharts_plots[div] = new Highcharts.Chart(highcharts_plot_options[div]);
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