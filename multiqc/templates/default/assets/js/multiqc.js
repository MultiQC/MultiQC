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
      if(plot_options.yAxis.title.text !== undefined){
        plot_options.yAxis.title.text = sym+plot_options.yAxis.title.text.substr(1);
      }
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
  
  // Make rows in general stats table sortable
  $('#general_stats_table tbody').sortable({ handle: '.sorthandle' });
  
  // Colour code table cells using chroma.js
  $('table').each(function(){
    var table = $(this);
    table.find('thead th').each(function(idx){
      var index = idx - 1;
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

  // Highlight colour filters
  var mqc_colours = chroma.brewer.Set1;
  var mqc_colours_idx = 0;
  $('#mqc_color_form').submit(function(e){
    e.preventDefault();
    var f_text = $('#mqc_colour_filter').val().trim();
    var f_col = $('#mqc_colour_filter_color').val().trim();
    var error = false;
    if(f_text.length == 0){ f_text = '[ all ]'; }
    $('#mqc_col_filters li .f_text').each(function(){
      if($(this).val() == f_text){
        alert('Error - highlight text "'+f_text+'" already exists');
        error = true;
      }
    });
    if(error){ return false; }
    $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_colours_idx+100)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    apply_mqc_highlights();
    $('#mqc_colour_filter').val('');
    mqc_colours_idx += 1;
    if(mqc_colours_idx >= mqc_colours.length){ mqc_colours_idx = 0; }
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
  });
  
  // Hide sample filters
  var mqc_hidesamples_idx = 0;
  $('#mqc_hidesamples_form').submit(function(e){
    e.preventDefault();
    var f_text = $('#mqc_hidesamples_filter').val().trim();
    if(f_text.length == 0){
      alert('Error - filter text must not be blank.');
      return false;
    }
    var error = false;
    $('#mqc_hidesamples_filters li .f_text').each(function(){
      if($(this).val() == f_text){
        alert('Error - highlight text "'+f_text+'" already exists');
        error = true;
      }
    });
    if(error){ return false; }
    $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_hidesamples_idx+100)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    apply_mqc_hidesamples();
    $('#mqc_hidesamples_filter').val('');
    mqc_hidesamples_idx += 1;
  });
  
  // Filter text is changed
  $('.mqc_filters').on('blur', 'li input', function(){
    var target = $(this).parent().parent().attr('id');
    if(target == 'mqc_col_filters'){
      if($(this).val().length == 0){ $(this).val('[ all ]'); }
      apply_mqc_highlights();
    }
    if(target == 'mqc_hidesamples_filters'){
      apply_mqc_hidesamples();
    }
  });
  // Enter key pressed whilst editing a filter
  $('.mqc_filters').on('keyup', 'li input', function(e){
    if(e.keyCode == 13) { // Pressed enter
      $(this).blur();
      $(this).parent().next('li').find('input').focus().select();
    }
  });
  // Remove filter button
  $('.mqc_filters').on('click', 'li button', function(){
    var target = $(this).parent().parent().attr('id');
    $(this).parent().remove();
    if(target == 'mqc_col_filters'){ apply_mqc_highlights(); }
    if(target == 'mqc_hidesamples_filters'){ apply_mqc_hidesamples(); }
  });
  // Use jQuery UI to make the colour filters sortable
  $("#mqc_col_filters").sortable();
  $("#mqc_col_filters").on("sortstop", function(event, ui){
    apply_mqc_highlights();
  });
  // Regex mode text
  $('.mqc_regex_mode').click(function(){
    if($(this).text() == 'Regex mode off'){
      $(this).html('Regex mode <ins><strong>on</strong></ins>');
    } else {
      $(this).html('Regex mode <strong>off</strong>');
    }
    if($(this).parent().attr('id') == 'mqc_cols'){ apply_mqc_highlights(); }
    if($(this).parent().attr('id') == 'mqc_hidesamples'){ apply_mqc_hidesamples(); }
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

// Basic Line Graph
function plot_xy_line_graph(div, data, config){
  if(config['tt_label'] === undefined){ config['tt_label'] = '{point.x}'; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  if (config['xDecimals'] === undefined){ config['xDecimals'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
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
      enabled: false
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
  if (config['stacking'] === undefined){ config['stacking'] = true; }
  if (config['use_legend'] === undefined){ config['use_legend'] = true; }
  if (config['yDecimals'] === undefined){ config['yDecimals'] = true; }
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
function apply_mqc_highlights(){
  
  // Collect the filters into an array
  var f_texts = [];
  var f_cols = [];
  var regex_mode = true;
  if($('#mqc_cols .mqc_regex_mode').text() == 'Regex mode off'){
    regex_mode = false;
  }
  $('#mqc_col_filters li .f_text').each(function(){
    f_texts.push($(this).val());
    f_cols.push($(this).css('color'));
  });  
  
  // Loop through each plot
  $.each(Object.keys(highcharts_plots), function(i, id){
    // Put in a try / catch so that one plot doesn't break all highlighting
    try {
      var p = highcharts_plots[id];
      // Colour the lines in line charts
      if(p.options.chart.type == 'line'){
        $.each(p.series, function(j, s){
          if(highcharts_plot_options[id]['series'][j]['color'] === undefined){
            highcharts_plot_options[id]['series'][j]['color'] = s.color;
          }
          var thiscol = highcharts_plot_options[id]['series'][j]['color'];
          $.each(f_texts, function(idx, f_text){
            if(f_text == '[ all ]'){ f_text = ''; }
            if((regex_mode && s.name.match(f_text)) || (!regex_mode && s.name.indexOf(f_text) > -1)){
              thiscol = f_cols[idx];
            }
          });
          s.color = thiscol;
          s.graph.attr({ stroke: thiscol });
          $.each(s.data, function(i, point) { point.pointAttr.hover.fill = thiscol; });
        });
      }
      
      // Colour the text labels on bar charts
      else if(p.options.chart.type == 'bar'){
        $(id+' .highcharts-xaxis-labels text').each(function(i, val){
          var labeltext = $(this).text();
          var thiscol = '#606060';
          $.each(f_texts, function(idx, f_text){
            if(f_text == '[ all ]'){ f_text = ''; }
            if((regex_mode && labeltext.match(f_text)) || (!regex_mode && labeltext.indexOf(f_text) > -1)){
              thiscol = f_cols[idx];
            }
          });
          $(this).css('fill', thiscol);
        });
      }
      p.redraw();
    } catch(err) {
      console.log("Error highlighting highcharts plot "+id)
    }
  });
  
  // Colour the heading text on the general stats table
  $('#general_stats_table tbody th').each(function(i){
    var thtext = $(this).text();
    var thiscol = '#333';
    $.each(f_texts, function(idx, f_text){
      if(f_text == '[ all ]'){ f_text = ''; }
      if((regex_mode && thtext.match(f_text)) || (!regex_mode && thtext.indexOf(f_text) > -1)){
        thiscol = f_cols[idx];
      }
    });
    $(this).css('color', thiscol);
  });
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_highlights', [f_texts, f_cols, regex_mode]);
  
}


// Apply the Highlight highlights to highcharts plots
function apply_mqc_hidesamples(){
  
  // Collect the filters into an array
  var f_texts = [];
  var regex_mode = true;
  if($('#mqc_hidesamples .mqc_regex_mode').text() == 'Regex mode off'){
    regex_mode = false;
  }
  $('#mqc_hidesamples_filters li .f_text').each(function(){
    f_texts.push($(this).val());
  });
  
  // Loop through each plot
  $('.hc-plot').each(function(){
    var plotid = $(this).attr('id');
    // Put in a try / catch so that one plot doesn't break all hiding
    try {
      // Line plots
      if($(this).highcharts().options.chart.type == 'line'){
        $.each($(this).highcharts().series, function(j, s){
          var match = false;
          $.each(f_texts, function(idx, f_text){
            if((regex_mode && s.name.match(f_text)) || (!regex_mode && s.name.indexOf(f_text) > -1)){
              match = true;
            }
          });
          if (s.visible && match) { s.hide(); }
          if (!s.visible && !match) { s.show(); }
        });
      }
      // Bar charts
      else if($(this).highcharts().options.chart.type == 'bar'){
        var replot = $.extend(true, [], highcharts_plot_options['#'+plotid]); // make a copy, not reference
        // Remove the categories
        var idx = replot.xAxis.categories.length;
        while (idx--) {
          var val = replot.xAxis.categories[idx];
          $.each(f_texts, function(i, f_text){
            if((regex_mode && val.match(f_text)) || (!regex_mode && val.indexOf(f_text) > -1)){
              replot.xAxis.categories.splice(idx, 1);
              $.each(replot.series, function(j, s){
                replot.series[j].data.splice(idx, 1);
              });
              return true;
            }
          });
        };
        highcharts_plots['#'+plotid] = new Highcharts.Chart(replot);
      }
    } catch(err) {
      console.log('Error hiding samples in '+$(this).attr('id')+' - '+err.message);
    }
  });
  
  // Hide rows in the general stats table
  $("#general_stats_table tbody th").each(function(){
    var match = false;
    var hfilter = $(this).text();
    $.each(f_texts, function(idx, f_text){
      if((regex_mode && hfilter.match(f_text)) || (!regex_mode && hfilter.indexOf(f_text) > -1)){
        match = true;
      }
    });
    if(match){ $(this).parent().hide(); }
    else { $(this).parent().show(); }
  });
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_hidesamples', [f_texts, regex_mode]);
  
}