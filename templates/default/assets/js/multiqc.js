/* Javascript for MultiQC Default Template */

var brewer_scales = ['￼YlOrRd', 'YlOrBr', 'YlGnBu', 'YlGn', 'Reds', 'RdPu',
  'Purples', 'PuRd', 'PuBuGn', 'PuBu', 'OrRd', 'Oranges', 'Greys', 'Greens',
  'GnBu', 'BuPu', 'BuGn', 'Blues', '￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼Set3', 'Set2', 'Set1', '￼￼￼Pastel2', 'Pastel1',
  'Paired', 'Dark2', 'Accent', '￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼Spectral', 'RdYlGn', 'RdYlBu', 'RdGy', 'RdBu',
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
  hc['sample_trimmed'] = hc['sample_read1'] = hc['sample_read2'] = 0;
  var rc = 0;
  $("#general_stats_table tbody tr").each(function(){
    var sn = $(this).find('th').text();
    rc += 1;
    if(sn.match(/_val_[12]$/) || sn.match(/_trimmed$/)){
      $(this).addClass('sample_trimmed');
      hc['sample_trimmed'] += 1;
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


})


////////////////////////////////////////////////
// HighCharts Plotting Functions
////////////////////////////////////////////////

// We store the config options for every graph in an arrays.
// This way, we can go back and change them via button clicks
// eg. Changing an axis from values to percentages
highcharts_plots = [];
highcharts_plot_options = [];
highcharts_plot_colors = [];

// Basic Line Graph
function plot_xy_line_graph(div, data, config){
  if(config['tt_label'] === undefined){ config['tt_label'] = '{point.x}'; }
  if(config['use_legend'] === undefined){ config['use_legend'] = true; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  highcharts_plot_options[div] = {
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
      max: config['xmax'],
      min: config['xmin']
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      max: config['ymax'],
      min: config['ymin'],
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
			enabled: false
		},
    tooltip: {
      headerFormat: '',
			pointFormat: '<span style="color:{series.color}; text-decoration:underline;">{series.name}</span><br>'+config['tt_label']+': {point.y:.2f}',
			useHTML: true
    },
    series: data
  }
  highcharts_plots[div] = $(div).highcharts(highcharts_plot_options[div]);
}

// Stacked Bar Graph
function plot_stacked_bar_graph(div, cats, data, config){
  if (config['colors'] === undefined){ config['colors'] = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9', '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']  }
  highcharts_plot_colors[div] = config['colors'];
  if (config['use_legend'] === undefined){ config['use_legend'] = true; }
  highcharts_plot_options[div] = {
    colors: config['colors'],
    chart: {
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
    },
    plotOptions: {
      series: {
        stacking: config['stacking'],
        groupPadding: 0.02
      }
    },
    credits: {
      enabled: false
    },
    legend: {
      enabled: config['use_legend'],
      reversed: true
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
  highcharts_plots[div] = $(div).highcharts(highcharts_plot_options[div]);
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
