/* Javascript for MultiQC Default Template */

var brewer_scales = ['YlOrRd', 'YlOrBr', 'YlGnBu', 'YlGn', 'Reds', 'RdPu',
  'Purples', 'PuRd', 'PuBuGn', 'PuBu', 'OrRd', 'Oranges', 'Greys', 'Greens',
  'GnBu', 'BuPu', 'BuGn', 'Blues', 'Set3', 'Set2', 'Set1', 'Pastel2', 'Pastel1',
  'Paired', 'Dark2', 'Accent', 'Spectral', 'RdYlGn', 'RdYlBu', 'RdGy', 'RdBu',
  'PuOr', 'PRGn', 'PiYG', 'BrBG'];

// Execute when page load has finished
$(function () {

  // Enable the bootstrap tooltip hovers
  // http://getbootstrap.com/javascript/#tooltips
  $('[data-toggle="tooltip"]').tooltip();

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
      plot_options.yAxis.title.text = sym+plot_options.yAxis.title.text.substr(1)
      plot_options.plotOptions.series.stacking = stack_type;
      $(target).highcharts(plot_options);
    }
    // Switch data source
    if(action == 'set_data'){
      var ds = $(this).data('newdata');
      plot_options.series = window[ds];
      var ylab = $(this).data('ylab');
      if(ylab !== undefined){ plot_options.yAxis.title.text = ylab; }
      var xlab = $(this).data('xlab');
      if(xlab !== undefined){ plot_options.xAxis.title.text = xlab; }
      $(target).highcharts(plot_options);
    }
  });

  // Show / Hide suspected duplicates in general stats table
  var hc = [];
  hc['sample_trimmed'] = hc['sample_untrimmed'] = hc['sample_read1'] = hc['sample_read2'] = 0;
  var rc = 0;
  $("#general_stats_table tbody tr").each(function(){
    var sn = $(this).find('th').text();
    rc += 1;
    if(sn.match(/_val_[12]$/) || sn.match(/_trimmed$/)){
      $(this).addClass('sample_trimmed');
      hc['sample_trimmed'] += 1;
    } else {
      $(this).addClass('sample_untrimmed');
      hc['sample_untrimmed'] += 1;
    }
    if(sn.match(/_1$/) || sn.match(/_R1$/i)){
      $(this).addClass('sample_read1');
      hc['sample_read1'] += 1;
    }
    if(sn.match(/_2$/) || sn.match(/_R2$/i)){
      $(this).addClass('sample_read2');
      hc['sample_read2'] += 1;
    }
  });
  if(hc['sample_trimmed'] > 0 && hc['sample_trimmed'] != rc){ $('.genstat_table_showhide[data-target="sample_trimmed"]').show(); }
  if(hc['sample_untrimmed'] > 0 && hc['sample_untrimmed'] != rc){ $('.genstat_table_showhide[data-target="sample_untrimmed"]').show(); }
  if(hc['sample_read1'] > 0 && hc['sample_read1'] != rc){ $('.genstat_table_showhide[data-target="sample_read1"]').show(); }
  if(hc['sample_read2'] > 0 && hc['sample_read2'] != rc){ $('.genstat_table_showhide[data-target="sample_read2"]').show(); }
  $('.genstat_table_showhide').click(function(e){
    e.preventDefault();
    $(this).toggleClass('active');
    showhide_general_stats_rows();
  });
  $('#genstat_table_showhide_custom').keyup(function(){
    showhide_general_stats_rows();
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
    var startHeight = wrapper.height();
    var pY = e.pageY;
    $(document).on('mouseup', function(e){
      $(document).off('mousemove');
    });
    $(document).on('mousemove', function(me){
      var my = (me.pageY - pY);
      wrapper.css('height', startHeight + my);
      $(document).resize();
    });
  });

  // Colour code table cells using chroma.js
  $('table').each(function(){
    var table = $(this);
    table.find('thead th').each(function(index){
      if($(this).hasClass('chroma-col')){

        // Get the colour scheme if set
        var colscheme_rev = false;
        var colscheme = $(this).data('chroma-scale');
        if(colscheme.substr(colscheme.length - 4) == '-rev'){
          colscheme_rev = true;
          colscheme = colscheme.substr(0, colscheme.length - 4);
        }
        if(colscheme === undefined || brewer_scales.indexOf(colscheme) == -1){
          colscheme = 'GnBu';
        }

        // Collect the data series
        var data = [];
        table.find('tr td:nth-of-type('+(index)+')').each(function(){
          var d = parseFloat($(this).text());
          data.push(d);
        });
        if(data.length == 0){ return true; } // Skip to next loop

        // Get the max and min values if not set with data attributes
        var maxval = $(this).data('chroma-max');
        var minval = $(this).data('chroma-min');
        if(maxval === undefined || minval === undefined){
          $.each(data, function(k, v){
            v = parseFloat(v);
            if(!isNaN(v)){
              if(v > maxval || maxval == undefined){ maxval = v; }
              if(v < minval || minval == undefined){ minval = v; }
            }
          });
        }
        if(isNaN(minval) || isNaN(maxval)){
          console.log('Could not calculate max or min value for '+$(this).text()+': ['+[minval, maxval]+']')
          return true; // Skip to next loop
        }

        // Go through table cells again, adding colour
        var i = 0;
        var scale = chroma.scale(colscheme).domain([minval, maxval]);
        if(colscheme_rev){
          scale = chroma.scale(colscheme).domain([maxval, minval]);
        }
        table.find('tr td:nth-of-type('+(index)+')').each(function(){
          var col = scale(parseFloat($(this).text())).css();
          $(this).css('background-color', col);
          // Change the colour of the text if we can get a better contrast ratio
          if(chroma.contrast('#EEEEEE', col) > chroma.contrast($(this).css('color'), col)){
            $(this).css('color', '#EEEEEE').css('font-weight', '200');
          }
        });
      }
    });
  });

  // Update colours of matching samples
  var hc_colours = chroma.brewer.Set1;
  var hc_colours_idx = 0;
  $('#hc_color_form').submit(function(e){
    e.preventDefault();
    var f_text = $('#hc_colour_filter').val().trim();
    var f_col = $('#hc_colour_filter_color').val().trim();
    var error = false;
    if(f_text.length == 0){ f_text = '[ all ]'; }
    $('#hc_col_filters li .f_text').each(function(){
      if($(this).text() == f_text){
        alert('Error - highlight text "'+f_text+'" already exists');
        error = true;
      }
    });
    if(!error){
      $('#hc_col_filters').append('<li style="color:'+f_col+';"><span class="f_text">'+f_text+'</span><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
      apply_hc_highlights();
      $('#hc_colour_filter').val('');
      hc_colours_idx += 1;
      if(hc_colours_idx >= hc_colours.length){ hc_colours_idx = 0; }
      $('#hc_colour_filter_color').val(hc_colours[hc_colours_idx]);
    }
  });
  $('#hc_col_filters').on('click', 'li', function(){
    $(this).remove();
    apply_hc_highlights();
  });
  // Use jQuery UI to make the filters sortable
  $("#hc_col_filters").sortable();
  $("#hc_col_filters").disableSelection();
  $("#hc_col_filters").on("sortstop", function(event, ui){ apply_hc_highlights(); });


})


////////////////////////////////////////////////
// HighCharts Plotting Functions
////////////////////////////////////////////////

// We store the config options for every graph in an arrays.
// This way, we can go back and change them via button clicks
// eg. Changing an axis from values to percentages
highcharts_plots = [];
highcharts_plot_options = [];

// Basic Line Graph
function plot_xy_line_graph(div, data, config){
  if(config['tt_label'] === undefined){ config['tt_label'] = '{point.x}'; }
  if(config['use_legend'] === undefined){ config['use_legend'] = true; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  if (config['xDecimals'] === undefined){ config['xDecimals'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  highcharts_plot_options[div] = {
    chart: {
      renderTo: div.replace('#',''),
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
      plotLines: [{
        value: 0,
        width: 1,
        color: '#808080'
      }]
    },
    plotOptions: {
      series: {
        marker: { enabled: false },
        cursor: 'pointer',
        point: {
          events: {
            click: config['click_func']
          }
        }
      }
    },
    legend: {
      enabled: config['use_legend'],
      layout: 'vertical',
      align: 'right',
      verticalAlign: 'middle',
      borderWidth: 0
    },
    credits: {
			enabled: true,
      text: 'Created with MultiQC',
      href: 'https://github.com/ewels/MultiQC'
		},
    tooltip: {
      headerFormat: '',
			pointFormat: '<span style="color:{series.color}; text-decoration:underline;">{series.name}</span><br>'+config['tt_label']+': {point.y:.2f}',
			useHTML: true
    },
    series: data
  }
  highcharts_plots[div] = new Highcharts.Chart(highcharts_plot_options[div]);
}

// Stacked Bar Graph
function plot_stacked_bar_graph(div, cats, data, config){
  if (config['use_legend'] === undefined){ config['use_legend'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
  highcharts_plot_options[div] = {
    chart: {
      renderTo: div.replace('#',''),
      type: 'bar',
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
      headerFormat: '',
      pointFormat: '<span style="color:{series.color}; text-decoration:underline;">{series.name}</span><br>'+config['tt_label']+': {point.y:.2f}',
      useHTML: true
    },
    tooltip: {
      formatter: function () {
        var s = '<table><tr><th colspan="3" style="font-weight:bold; text-decoration:underline;">' + this.x + '</th></tr>';
        $.each(this.points.reverse(), function () {
          s += '<tr><td style="font-weight:bold; color:'+this.series.color+'; padding-right: 15px;">' + this.series.name + ':</td><td style="text-align:right;">' + numberWithCommas(this.y) + ' reads</td><td style="text-align:right;">(' + this.percentage.toFixed(1) + '%)</td></tr>';
        });
        s += '</table>';
        return s;
      },
      shared: true,
      useHTML: true
    },
    series: data
  }
  highcharts_plots[div] = new Highcharts.Chart(highcharts_plot_options[div]);
}



//////////////////
// Generic helper functions

// From http://stackoverflow.com/questions/2901102/how-to-print-a-number-with-commas-as-thousands-separators-in-javascript
function numberWithCommas(x) {
  return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

// http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
function findPos(obj) {
  var curleft = 0, curtop = 0;
  if (obj.offsetParent) {
    do {
      curleft += obj.offsetLeft;
      curtop += obj.offsetTop;
    } while (obj = obj.offsetParent);
    return { x: curleft, y: curtop };
  }
  return undefined;
}

// Show / Hide rows in general stats table
function showhide_general_stats_rows(){
  var showhideclasses = [];
  var hideclasses = [];
  $('.genstat_table_showhide').each(function(){
    showhideclasses.push($(this).data('target'));
    if($(this).hasClass('active')){
      hideclasses.push($(this).data('target'));
    }
  });
  $.each(showhideclasses, function(k, v){
    $('.'+v).show();
  });
  $.each(hideclasses, function(k, v){
    $('.'+v).hide();
  });
  var custom_text = $('#genstat_table_showhide_custom').val().trim();
  if(custom_text.length > 0){
    $('#general_stats_table tbody tr').each(function(){
      var sn = $(this).find('th').text();
      if(sn.indexOf(custom_text) > -1){
        $(this).hide();
      }
    });
  }
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

// Apply the Highlight highlights to highcharts plots
function apply_hc_highlights(){
  // First - reset colours on all plots
  $.each(Object.keys(highcharts_plots), function(i, id){
    var p = highcharts_plots[id];
    if(p.options.chart.type == 'line'){
      $.each(p.series, function(j, s){
        var col = highcharts_plot_options[id]['series'][j]['color'];
        if (col != undefined){
          s.color = col;
          s.graph.attr({ stroke: col });
          $.each(s.data, function(i, point) { point.pointAttr.hover.fill = col; });
        }
      });
    } else if(p.options.chart.type == 'bar'){
      $(id+' .highcharts-axis-labels text').css('fill', '#606060');
    }
  });
  $('#general_stats_table tbody th').css('color', '#333');

  // Loop through each filter
  $('#hc_col_filters li .f_text').each(function(){
    var f_text = $(this).text();
    var f_col = $(this).css('color');
    if(f_text == '[ all ]'){ f_text = ''; }
    // Loop through each plot
    $.each(Object.keys(highcharts_plots), function(i, id){
      var p = highcharts_plots[id];
      if(p.options.chart.type == 'line'){
        $.each(p.series, function(j, s){
          if(highcharts_plot_options[id]['series'][j]['color'] === undefined){
            highcharts_plot_options[id]['series'][j]['color'] = s.color;
          }
          if(s.name.indexOf(f_text) > -1){
            s.color = f_col;
            s.graph.attr({ stroke: f_col });
            $.each(s.data, function(i, point) { point.pointAttr.hover.fill = f_col; });
          }
        });
      } else if(p.options.chart.type == 'bar'){
        $(id+' .highcharts-axis-labels text:contains('+f_text+')').css('fill', f_col);
      }
    });
    $('#general_stats_table tbody th:contains('+f_text+')').css('color', f_col);
  });
  $.each(Object.keys(highcharts_plots), function(i, id){
    highcharts_plots[id].redraw();
  });
}
