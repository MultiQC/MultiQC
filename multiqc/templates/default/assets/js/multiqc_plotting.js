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
  
  // Show the overlay plots again (clicking the original)
  $('.original-plot').click(function(){
    $(this).closest('.showhide_orig').next('.hc-plot').slideDown();
    $(this).closest('.showhide_orig').slideUp();
    $(this).closest('.mqc-section').find('.instr').text('Click to show original FastQC plot.');
    highlight_fade_text($(this).closest('.mqc-section').find('.instr'));
  });

  // prev / next buttons for original images
  $('.original_plot_prev_btn, .original_plot_nxt_btn').click(function(e){
    e.preventDefault();
    var name = $(this).attr('href').substr(1);
    var target = $(this).data('target').substr(1);
    hc_original_chg_source (name, target);
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
  if (config['orig_click_func'] === true){
    config['cursor'] = 'pointer';
    config['click_func'] = function (e) {
      var id = e.toElement.offsetParent.offsetParent.id;
      var p = $('#'+id).parent();
      if(id !== undefined){
        hc_original_chg_source (this.series.name, id);
        p.find('.showhide_orig').delay(100).slideDown();
        p.find('.hc-plot').delay(100).slideUp();
      }
      // Fire off a custom jQuery event for other javascript chunks to tie into
      $('#'+id).trigger('mqc_original_series_click', [this.series.name]);
    }
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
      max: config['xmax'],
      min: config['xmin'],
      allowDecimals: config['xDecimals'],
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      max: config['ymax'],
      min: config['ymin'],
      allowDecimals: config['yDecimals'],
      plotBands: config['plotBands']
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
			pointFormat: '<span style="color:{series.color}; text-decoration:underline;">{series.name}</span><br>'+config['tt_label'],
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

// Update an original plot source
function hc_original_chg_source (name, id) {
  
  // Get the plot image paths
  try {
    var names = eval(id+'_orig_plots');
  } catch(err) {
    console.log("Couldn't parse original plot paths array for "+id+'_orig_plots - '+err);
    return false;
  }
  if(names.length == 0){
    console.log("Couldn't find original plot paths array for "+id+'_orig_plots - '+err);
    return false;
  }
  
  var target = $('#'+id).parent();
  
  // Wipe the src in case it's not found later
  target.find('img.original-plot').attr('src','assets/img/img_not_found.png');
  
  // Find the original sample name if it's been renamed
  s_name = get_orig_name(name, id);
  
  // Find the image source
  var src = undefined;
  var i;
  $.each(names, function(idx, n){
    if(n['s_name'] == s_name){
      src = n['img_path'];
      i = idx;
      return false;
    }
  });
  if(src !== undefined){
    target.find('img.original-plot').attr('src', src);
    target.find('.s_name').text(get_new_name(name, id, s_name));

    var l = names.length;
    var n_i = i+1 < l ? i+1 : 0;
    var p_i = i-1 >= 0 ? i-1 : l - 1;
    // Sample names for next / prev links. Not renamed, but should be ok.
    var n = names[n_i]['s_name'];
    var p = names[p_i]['s_name'];
    target.find('.original_plot_prev_btn').attr('href', '#'+p);
    target.find('.original_plot_nxt_btn').attr('href', '#'+n);
    if(target.closest('.mqc-section').find('.instr').text() != "Click plot to return to overview plot."){
      target.closest('.mqc-section').find('.instr').text("Click plot to return to overview plot.");
      highlight_fade_text(target.closest('.mqc-section').find('.instr'));
    }
  } else {
    console.log("Couldn't find original image path for "+s_name+", id:"+id);
  }
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  target.trigger('mqc_original_chg_source', [name]);
  
}
