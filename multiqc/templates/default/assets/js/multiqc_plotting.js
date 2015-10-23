////////////////////////////////////////////////
// HighCharts Plotting Code
////////////////////////////////////////////////

mqc_highcharts = [];

// Execute when page load has finished loading
$(function () {
  
  // Render a plot when clicked
  $('body').on('click', '.render_plot', function(e){
    var target = $(this).parent().attr('id');
    plot_graph(target);
  });
  
  // Switch a HighCharts axis or data source
  $('.switch_group button').click(function(e){
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
      plot_graph(target);
    }
    // Switch data source
    if(action == 'set_data'){
      var ds = $(this).data('newdata');
      plot_graph(target, ds);
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
        $('#'+target).addClass('not_rendered');
      }
    }
    // Bar graphs
    else if(mqc_plots[target]['plot_type'] == 'bar_graph'){
      if(max_num === undefined || mqc_plots[target]['samples'][0].length < max_num){
        plot_stacked_bar_graph(target, ds);
        $('#'+target).removeClass('not_rendered');
      } else {
        $('#'+target).addClass('not_rendered');
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
  else { if(config['cursor'] === undefined){ config['cursor'] = 'pointer'; } }
  if (config['xDecimals'] === undefined){ config['xDecimals'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  if (config['pointFormat'] === undefined){
    config['pointFormat'] = '<div style="background-color:{series.color}; display:inline-block; height: 10px; width: 10px; border:1px solid #333;"></div> <span style="text-decoration:underline; font-weight:bold;">{series.name}</span><br>'+config['tt_label'];
  }

  // Make the highcharts plot
  mqc_highcharts[target] = new Highcharts.Chart({
    chart: {
      renderTo: $('#'+target)[0],
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
      floor: config['xFloor'],
      max: config['xmax'],
      min: config['xmin'],
      minRange: config['xMinRange'],
      allowDecimals: config['xDecimals'],
      plotBands: config['xPlotBands']
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      ceiling: config['yCeiling'],
      floor: config['yFloor'],
      max: config['ymax'],
      min: config['ymin'],
      minRange: config['yMinRange'],
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
    series: data[ds]
  });
}

// Stacked Bar Graph
function plot_stacked_bar_graph(target, ds){
  if(mqc_plots[target] === undefined || mqc_plots[target]['plot_type'] !== 'bar_graph'){
    return false;
  }
  var config = mqc_plots[target]['config'];
  var data = mqc_plots[target]['datasets'];
  var cats = mqc_plots[target]['samples'];
  if(ds === undefined){ ds = 0; }
  
  if (config['stacking'] === undefined){ config['stacking'] = 'normal'; }
  if (config['use_legend'] === undefined){ config['use_legend'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  else { if(config['cursor'] === undefined){ config['cursor'] = 'pointer'; } }
  
  console.log('plotting to '+ '#'+target)
  
  // Make the highcharts plot
  mqc_highcharts[target] = new Highcharts.Chart({
    chart: {
      renderTo: $('#'+target)[0],
      type: 'bar',
      backgroundColor: null,
    },
    title: {
      text: config['title'],
    },
    xAxis: {
      categories: cats[ds],
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
        groupPadding: 0.02,
        borderWidth: 2
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
          yval = this.y.toFixed(0)
          s += '<tr><td style="font-weight:bold; color:'+this.series.color+'; padding-right: 15px; border-bottom:1px solid #dedede;">' + this.series.name + ':</td><td style="text-align:right; border-bottom:1px solid #dedede;">' + numberWithCommas(yval) + '</td><td style="text-align:right; border-bottom:1px solid #dedede;">(' + this.percentage.toFixed(1) + '%)</td></tr>';
        });
        s += '</table>';
        return s;
      },
      shared: true,
      useHTML: true
    },
    series: data[ds]
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