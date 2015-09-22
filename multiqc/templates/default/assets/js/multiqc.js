/* Javascript for MultiQC Default Template */

var brewer_scales = ['YlOrRd', 'YlOrBr', 'YlGnBu', 'YlGn', 'Reds', 'RdPu',
  'Purples', 'PuRd', 'PuBuGn', 'PuBu', 'OrRd', 'Oranges', 'Greys', 'Greens',
  'GnBu', 'BuPu', 'BuGn', 'Blues', 'Set3', 'Set2', 'Set1', 'Pastel2', 'Pastel1',
  'Paired', 'Dark2', 'Accent', 'Spectral', 'RdYlGn', 'RdYlBu', 'RdGy', 'RdBu',
  'PuOr', 'PRGn', 'PiYG', 'BrBG'];
var mqc_colours_idx = 0;
var mqc_colours = chroma.brewer.Set1;

// Execute when page load has finished
$(function () {

  // Enable the bootstrap tooltip hovers
  // http://getbootstrap.com/javascript/#tooltips
  $('[data-toggle="tooltip"]').tooltip();
  
  // Enable tablesorter on the general statistics table
  $("#general_stats_table").tablesorter({sortInitialOrder: 'desc'}); 
  
  // Introduction tour
  $('#mqc-launch-into-tour').after(' <small><em>('+tour_steps.length+' steps - takes around 1 minute)</em></small>');
  $.each(tour_steps, function(i, step){
    step['title'] += '<span class="pull-right">'+(i+1)+'/'+tour_steps.length+'</span>';
    var percent = parseInt(((i+1) / tour_steps.length) * 100);
    step['content'] = '<div class="pbar_wrapper"><hr class="pbar" style="width:'+percent+'%;"></div>' + step['content'];
  });
  var intro_tour = new Tour({
    backdrop: true,
    storage: false,
    steps: tour_steps
  });
  intro_tour.init();
  $('#mqc-launch-into-tour').click(function(){
    try{ intro_tour.start(); }
    catch(e){ console.log('Tour broke - '+e); }
  });
  
  // Load any saved configuration
  load_mqc_config();
  

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
  
  // Make rows in general stats table sortable
  $('#general_stats_table tbody').sortable({
    handle: '.sorthandle',
    helper: function fixWidthHelper(e, ui) {
      ui.children().each(function() { $(this).width($(this).width()); });
      return ui;
    }
  });
  
  // Colour code table cells using chroma.js
  $('table').each(function(){
    var table = $(this);
    table.find('thead th').each(function(idx){
      var index = idx + 1;
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
          var val = parseFloat($(this).text());
          var col = scale(val).css();
          var percentage = ((val - minval) / (maxval - minval)) * 100;
          $(this).html('<div class="wrapper"><span class="bar" style="width:'+percentage+'%; background-color: '+col+';"></span><span class="val">'+$(this).text()+'</span></div>');
          $(this).addClass('data-coloured');
        });
        
        // Add to the key modal
        var text = $(this).text();
        var h_text = $(this).find('span').attr('title');
        if(String(h_text).length == 0){
          var h_text = $(this).find('span').data('original-title');
        }
        // Text Colours
        var max_font_col = '#333';
        var min_font_col = '#333';
        if(chroma.contrast('#eeeeee', scale(maxval).css()) > chroma.contrast('#333333', scale(maxval).css())){ max_font_col = '#EEE'; }
        if(chroma.contrast('#eeeeee', scale(minval).css()) > chroma.contrast('#333333', scale(minval).css())){ min_font_col = '#EEE'; }
        // Make really complicated gradient string for this colour scheme
        // Note - NOT simple linear gradients!
        var grad_cols = [];
        var wk_grad_cols = [];
        for (i = 0; i <= 100; i+= 10) {
          var val = ((i/100) * (maxval - minval)) + minval;
          grad_cols.push(scale(val).css()+' '+i+'%');
          wk_grad_cols.push('color-stop('+i+'%, '+scale(val).css()+')');
        }
        grad_cols = grad_cols.join();
        wk_grad_cols = wk_grad_cols.join();
        var grad_k = '<div class="general_stats_key_bar" style="color: '+min_font_col+'; ' +
                        'background: -moz-linear-gradient(left, '+grad_cols+'); '+
                        'background: -webkit-gradient(linear, left top, right top, '+wk_grad_cols+'); '+
                        'background: -webkit-linear-gradient(left, '+grad_cols+'); '+
                        'background: -o-linear-gradient(left, '+grad_cols+'); '+
                        'background: -ms-linear-gradient(left, '+grad_cols+'); '+
                        'background: linear-gradient(to right, '+grad_cols+'); '+
                        '"><span style="float:right; color: '+max_font_col+';">'+maxval+'</span>'+minval+'</div>';
        var tr = $('<tr><th>'+text+'</th><td>'+h_text+'</td><td>'+grad_k+'</td></tr>');
        $('#general_stats_key_modal table tbody').append(tr);
      }
    });
  });
  
  // Toolbox buttons
  $('.mqc-toolbox-buttons a').click(function(e){
    e.preventDefault();
    var target = $(this).attr('href');
    mqc_toolbox_openclose(target);
  });
  // Shortcut keys
  $(window).keydown(function (evt) {
    if (evt.target.tagName.toLowerCase() !== 'input' && evt.target.tagName.toLowerCase() !== 'textarea') {
      // console.log(evt.which);
      if(evt.which == 67){ mqc_toolbox_openclose('#mqc_cols'); } // c
      if(evt.which == 82){ mqc_toolbox_openclose('#mqc_renamesamples'); } // r
      if(evt.which == 72){ mqc_toolbox_openclose('#mqc_hidesamples'); } // h
      if(evt.which == 83){ mqc_toolbox_openclose('#mqc_saveconfig'); } // r
    }
  });

  // Highlight colour filters
  $('#mqc_color_form').submit(function(e){
    e.preventDefault();
    var f_text = $('#mqc_colour_filter').val().trim();
    var f_col = $('#mqc_colour_filter_color').val().trim();
    if(f_text.length == 0){
      alert('Error - highlight text must not be blank.');
      return false;
    }
    $('#mqc_col_filters li .f_text').each(function(){
      if($(this).val() == f_text){
        alert('Error - highlight text "'+f_text+'" already exists');
        return false;
      }
    });
    $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_colours_idx)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    apply_mqc_highlights();
    $('#mqc_colour_filter').val('');
    mqc_colours_idx += 1;
    if(mqc_colours_idx >= mqc_colours.length){ mqc_colours_idx = 0; }
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
  });
  
  // Hide sample filters
  var mqc_hidesamples_idx = 200;
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
    $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_hidesamples_idx)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    apply_mqc_hidesamples();
    $('#mqc_hidesamples_filter').val('');
    mqc_hidesamples_idx += 1;
  });
  
  // Rename samples
  var mqc_renamesamples_idx = 300;
  $('#mqc_renamesamples_form').submit(function(e){
    e.preventDefault();
    var from_text = $('#mqc_renamesamples_from').val().trim();
    var to_text = $('#mqc_renamesamples_to').val().trim();
    if(from_text.length == 0){
      alert('Error - "From" text must not be blank.');
      return false;
    }
    var li = '<li><input class="f_text from_text" value="'+from_text+'" tabindex="'+(mqc_renamesamples_idx)+'" />'
    li += '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="'+to_text+'" tabindex="'+(mqc_renamesamples_idx+1)+'" />'
    li += '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>'
    $('#mqc_renamesamples_filters').append(li);
    apply_mqc_renamesamples();
    $('#mqc_renamesamples_from').val('');
    $('#mqc_renamesamples_to').val('');
    mqc_renamesamples_idx += 2;
    $('#mqc_renamesamples_form input:first').focus();
  });
  
  // Bulk rename samples
  $('#mqc_renamesamples_bulk_collapse').on('shown.bs.collapse', function () {
    $('#mqc_renamesamples_bulk_form textarea').focus();
  });
  $('#mqc_renamesamples_bulk_form').submit(function(e){
    e.preventDefault();
    var raw = $(this).find('textarea').val();
    var lines = raw.match(/^.*([\n\r]+|$)/gm);
    $.each(lines, function(i, l){
      var sections = l.split("\t", 2);
      if(sections.length < 2){ return true; }
      var from_text = sections[0].trim();
      var to_text = sections[1].trim();
      if(from_text.length == 0){ return true; }
      var li = '<li><input class="f_text from_text" value="'+from_text+'" tabindex="'+(mqc_renamesamples_idx)+'" />'
      li += '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="'+to_text+'" tabindex="'+(mqc_renamesamples_idx+1)+'" />'
      li += '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>'
      $('#mqc_renamesamples_filters').append(li);
    });
    apply_mqc_renamesamples();
    $(this).find('textarea').val('');
    $('#mqc_renamesamples_bulk_collapse').collapse('hide');
  });
  
  // Save config
  $('.mqc_saveconfig_btn').click(function(e){
    e.preventDefault();
    if($(this).data('target') == 'autosave'){
      if($(this).text() == 'Auto-save is off'){
        $(this).html('Auto-save is <strong>on</strong>');
      } else {
        mqc_save_config(undefined, undefined, true); // clear autosave
        $(this).html('Auto-save is <strong>off</strong>');
      }
    }
    if($(this).data('target') == 'local'){
      mqc_save_config('general', 'local');
      $('<p class="text-success" id="mqc-save-general-success">General config saved.</p>').hide().insertAfter($(this).parent()).slideDown(function(){
        setTimeout(function(){
          $('#mqc-save-general-success').slideUp(function(){ $(this).remove(); });
        }, 5000);
      });
    }
    if($(this).data('target') == 'file'){
      mqc_save_config('general', 'file');
    }
  });
  
  // Filter text is changed
  $('.mqc_filters').on('blur', 'li input', function(){
    var target = $(this).parent().parent().attr('id');
    if(target == 'mqc_col_filters'){
      apply_mqc_highlights();
    }
    if(target == 'mqc_hidesamples_filters'){
      apply_mqc_hidesamples();
    }
    if(target == 'mqc_renamesamples_filters'){
      apply_mqc_renamesamples();
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
    if(target == 'mqc_renamesamples_filters'){ apply_mqc_renamesamples(); }
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

// Toolbox functions
function mqc_toolbox_openclose (target, open){
  var btn = $('.mqc-toolbox-buttons a[href="'+target+'"]');
  if(open === undefined){
    if(btn.hasClass('active')){ open = false; }
    else { open = true; }
  }
  var already_open = $('.mqc-toolbox').hasClass('active');
  if(open){
    $('.mqc-toolbox, .mqc-toolbox-buttons a, .mqc_filter_section').removeClass('active');
    btn.addClass('active');
    $('.mqc-toolbox, .mqc-toolbox-label, '+target).addClass('active');
    $(document).trigger('mqc_toolbox_open');
    var timeout = already_open ? 0 : 510;
    setTimeout(function(){
      if(target == '#mqc_cols'){ $('#mqc_colour_filter').focus(); }
      if(target == '#mqc_renamesamples'){ $('#mqc_renamesamples_from').focus(); }
      if(target == '#mqc_hidesamples'){ $('#mqc_hidesamples_filter').focus(); }
    }, timeout);
  } else {
    btn.removeClass('active');
    $('.mqc-toolbox, .mqc-toolbox-buttons a, .mqc-toolbox-label').removeClass('active');
    $(document).trigger('mqc_toolbox_close');
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
  
  // Apply a 'background' highlight to remove default colouring first
  if(f_texts.length > 0){
    f_texts.unshift('');
    f_cols.unshift('#cccccc');
  }
  
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
      if((regex_mode && thtext.match(f_text)) || (!regex_mode && thtext.indexOf(f_text) > -1)){
        thiscol = f_cols[idx];
      }
    });
    $(this).css('color', thiscol);
  });
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_highlights', [f_texts, f_cols, regex_mode]);
  mqc_autosave();
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
      if($(this).highcharts() === undefined){ return true; }
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
        var matches = 0;
        // Remove the categories
        var idx = replot.xAxis.categories.length;
        while (idx--) {
          var val = replot.xAxis.categories[idx];
          $.each(f_texts, function(i, f_text){
            if((regex_mode && val.match(f_text)) || (!regex_mode && val.indexOf(f_text) > -1)){
              matches += 1;
              replot.xAxis.categories.splice(idx, 1);
              $.each(replot.series, function(j, s){
                replot.series[j].data.splice(idx, 1);
              });
              return true;
            }
          });
        };
        if(matches > 0){ highcharts_plots['#'+plotid] = new Highcharts.Chart(replot); }
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
  mqc_autosave();
}

function apply_mqc_renamesamples(){
  // Loop through each plot
  $('.hc-plot').each(function(){
    // Skip non-standard plot types
    if($(this).highcharts() === undefined){
      return true;
    }
    // Put in a try / catch so that one plot doesn't break all hiding
    try {
      var pid = $(this).attr('id');
      // Line plots
      if($(this).highcharts().options.chart.type == 'line'){
        $.each($(this).highcharts().series, function(j, s){
          var orig_name = undefined;
          try {
            orig_name = highcharts_plot_options['#'+pid]['s_names'][j];
          } catch(err) {}
          s.name = get_new_name(s.name, pid, orig_name);
        });
      }
      // Bar charts
      else if($(this).highcharts().options.chart.type == 'bar'){
        var replot = $.extend(true, [], highcharts_plot_options['#'+pid]); // make a copy, not reference
        // Rename the categories
        $.each(replot.xAxis.categories, function(idx, val){
          var s_name = replot.xAxis.categories[idx];
          var orig_name = undefined;
          try {
            orig_name = replot['s_names'][idx];
          } catch(err) {}
          var n_name = get_new_name(s_name, pid, orig_name);
          if(n_name != s_name){
            replot.xAxis.categories[idx] = n_name;
          }
        });
        highcharts_plots['#'+pid] = new Highcharts.Chart(replot);
      }
    } catch(err) {
      console.log('Error renaming samples in '+$(this).attr('id')+' - '+err.message);
    }
  });
  
  // Rename original plot titles
  $('.hc-plot-wrapper').each(function(){
    var id = $(this).attr('id');
    var t = $(this).find('.s_name');
    t.text(get_new_name(t.text(), id));
  });
  
  // Rename samples in the general stats table
  $("#general_stats_table tbody th").each(function(){
    var s_name = $(this).text();
    var orig_name = $(this).data('original-sample-name');
    if(orig_name === undefined){
      $(this).data('original-sample-name', s_name);
    }
    $(this).text(get_new_name(s_name, '#general_stats_table', orig_name));
  });
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_renamesamples');
  mqc_autosave();
}
// Try to keep track of what samples are called
// This dict works as a backup for non-standard plots
mqc_renamed_samples = {};
function get_new_name(s_name, obj_id, orig_name){
  // Collect the filters into an array
  var f_texts = [];
  var t_texts = [];
  $('#mqc_renamesamples_filters .from_text').each(function(){
    f_texts.push($(this).val());
  });
  $('#mqc_renamesamples_filters .to_text').each(function(){
    t_texts.push($(this).val());
  });
  // Find the original name from what we currently have, if not supplied
  if(orig_name === undefined){
    orig_name = get_orig_name(s_name, obj_id);
  }
  s_name = orig_name;
  // Now rename it
  $.each(f_texts, function(idx, f_text){
    s_name = s_name.replace(f_text, t_texts[idx]);
  });
  // Update list of names for this plot
  if(mqc_renamed_samples[obj_id] === undefined){
    mqc_renamed_samples[obj_id] = {};
  }
  mqc_renamed_samples[obj_id][orig_name] = s_name;
  return s_name;
}
// Try to guess the original sample name, based on the contents
// of mqc_renamed_samples
function get_orig_name(s_name, obj_id){
  if(mqc_renamed_samples[obj_id] === undefined){
    return s_name;
  }
  $.each(mqc_renamed_samples[obj_id], function(s, r){
    if(r == s_name){
      s_name = s;
      return false; // end loop
    }
  });
  return s_name;
}

// Autosave function
function mqc_autosave(){
  if($('#mqc_saveconfig_autosave').text() == 'Auto-save is off'){
    return;
  }
  mqc_save_config();
}

// Save the current configuration setup
function mqc_save_config(target, method, clear){
  if(target === undefined){ target = $('#mqc_output_path').text(); }
  if(method === undefined){ method = 'local'; }
  var config = {};
  // Autosave status
  config['autosave_tools'] = true;
  if($('#mqc_saveconfig_autosave').text() == 'Auto-save is off'){
    config['autosave_tools'] = false;
  }
  // Collect the highlight filters
  config['highlight_regex'] = true;
  if($('#mqc_cols .mqc_regex_mode').text() == 'Regex mode off'){
    config['highlight_regex'] = false;
  }
  config['highlights_f_texts'] = [];
  config['highlights_f_cols'] = [];
  $('#mqc_col_filters li .f_text').each(function(){
    config['highlights_f_texts'].push($(this).val());
    config['highlights_f_cols'].push($(this).css('color'));
  });
  // Collect the hide sample filters
  config['hidesamples_regex'] = true;
  if($('#mqc_hidesamples .mqc_regex_mode').text() == 'Regex mode off'){
    config['hidesamples_regex'] = false;
  }
  config['hidesamples_f_texts'] = [];
  $('#mqc_hidesamples_filters li .f_text').each(function(){
    config['hidesamples_f_texts'].push($(this).val());
  });
  // Collect the rename filters
  config['rename_from_texts'] = [];
  config['rename_to_texts'] = [];
  $('#mqc_renamesamples_filters .from_text').each(function(){
    config['rename_from_texts'].push($(this).val());
  });
  $('#mqc_renamesamples_filters .to_text').each(function(){
    config['rename_to_texts'].push($(this).val());
  });
  
  if(method == 'local'){
    var prev_config = {};
    // Load existing configs (inc. from other reports)
    try {
      prev_config = localStorage.getItem("mqc_config");
      if(prev_config !== null && prev_config !== undefined){
        prev_config = JSON.parse(prev_config);
      } else {
        prev_config  = {};
      }
    } catch(e){ console.log('Error updating localstorage: '+e); }
    // Update config obj with current config
    if(clear == true){
      prev_config[target] = {'autosave_tools': false};
    } else {
      prev_config[target] = config;
      prev_config[target]['last_updated'] = Date();
    }
    localStorage.setItem("mqc_config", JSON.stringify(prev_config));
  }
  if(method == 'file'){
    var f = "// Config file for MultiQC\n// https://github.com/ewels/MultiQC\n";
    f += "// Generated "+Date()+"\n// Original report path: "+$('#mqc_output_path').text()+"\n";
    f += "\nmqc_config_file_cfg = "+JSON.stringify(config, null, '  ');
    var fblob = new Blob([f], {type: "application/javascript;charset=utf-8"});
    saveAs(fblob, "multiqc_config.js");
    setTimeout(function(){
      alert('Now move multiqc_config.js from your downloads folder to the directory where this report is located.');
    }, 500);
  }
}

// Load previously saved config variables
function load_mqc_config(){
  var config = {};
  // Get local configs - general and then this path
  try {
    var local_config = localStorage.getItem("mqc_config");
    if(local_config !== null && local_config !== undefined){
      local_config = JSON.parse(local_config);
      if(local_config['general'] !== undefined){ console.log('Loaded local general config'); }
      for (var attr in local_config['general']) {
        config[attr] = local_config['general'][attr];
      }
      var path = $('#mqc_output_path').text();
      if(local_config[path] !== undefined){ console.log('Loaded local report config'); }
      for (var attr in local_config[path]) {
        config[attr] = local_config[path][attr];
      }
    }
  } catch(e){ console.log('Could not load local config: '+e); }
  
  // Local file
  if(typeof mqc_config_file_cfg != "undefined"){
    for (var attr in mqc_config_file_cfg) {
      config[attr] = mqc_config_file_cfg[attr];
    }
    $('#mqc_report_location').after('<p>MultiQC report toolbox config loaded from file.</p>');
    console.log('Loaded config from a file');
  }
  var update_highlights = false;
  var update_rename = false;
  var update_hide = false;
  // Apply config - highlights
  if(notEmptyObj(config['highlights_f_texts']) && notEmptyObj(config['highlights_f_cols'])){
    $.each(config['highlights_f_texts'], function(idx, f_text){
      var f_col = config['highlights_f_cols'][idx];
      $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
      mqc_colours_idx += 1;
    });
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
    if(config['highlight_regex'] == true){
      $('#mqc_cols .mqc_regex_mode').html('Regex mode <strong>off</strong>');
    }
    update_highlights = true;
  }
  // Rename samples
  if(notEmptyObj(config['rename_from_texts']) && notEmptyObj(config['rename_to_texts'])){
    $.each(config['rename_from_texts'], function(idx, from_text){
      var to_text = config['rename_to_texts'][idx];
      if(from_text.length == 0){ return true; }
      var li = '<li><input class="f_text from_text" value="'+from_text+'" />'
      li += '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="'+to_text+'" />'
      li += '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>'
      $('#mqc_renamesamples_filters').append(li);
    });
    if(config['hidesamples_regex'] == true){
      $('#mqc_hidesamples .mqc_regex_mode').html('Regex mode <strong>off</strong>');
    }
    update_rename = true;
  }
  // Hide samples
  if(notEmptyObj(config['hidesamples_f_texts'])){
    $.each(config['hidesamples_f_texts'], function(idx, f_text){
      if(f_text.length == 0){ return true; }
      $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    });
    update_hide = true;
  }
  // Autosave
  if(config['autosave_tools'] == false){
    $('#mqc_saveconfig_autosave').html('Auto-save is <strong>off</strong>');
  }
  // Wait for the rest of the page to render, then apply changes.
  // This is ugly. Can anyone think of a better way?
  setTimeout(function(){
    if(update_highlights){ apply_mqc_highlights(); }
    if(update_rename){ apply_mqc_renamesamples(); }
    if(update_hide){ apply_mqc_hidesamples(); }
  }, 1000);
  
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

// Helper config - is defined and object length > 0?
function notEmptyObj (obj){
  try{
    if(obj === undefined){ return false; }
    if(obj.length == 0){ return false; }
  } catch(e){ return false; }
  return true;
}




// MultiQC Tour, using bootstrap tour
var tour_steps = [
  {
    orphan: true,
    title: "Welcome to MultiQC!",
    content: 'MultiQC generates reports based on analysis across many samples.<br>Click next to explore the tools availble in this report.',
    onShow: function(tour){
      mqc_toolbox_openclose('#mqc_cols', false);
      $('.mqc-toolbox').css('z-index', 0);
    }
  },
  {
    element: "#general_stats",
    placement: 'top',
    title: 'General Statistics',
    backdropPadding: {'left': 10},
    content: "This table gives an overview of your samples, with data from all modules.",
  },
  {
    element: "#general_stats_table thead tr th:nth-child(3)",
    title: "Sort Columns",
    content: "Click a header to sort by that column, shift-click to sort by multiple."
  },
  {
    element: "#general_stats_table tbody tr:first-child td:first-child",
    title: "Reorder rows",
    content: "Drag the handle on any row to rearrange."
  },
  {
    element: ".hc-line-plot:first",
    placement: 'top',
    title: "Plots",
    content: "Plots are dynamic and interactive - click and drag line charts to zoom."
  },
  {
    element: ".switch_group:first",
    title: "Counts / Percentages",
    backdropPadding: 5,
    content: "Bar plots can switch between counts and percentages."
  },
  {
    element: ".hc-bar-plot:first",
    title: "Show / Hide Categories",
    placement: 'top',
    content: "Click category names in the legend to hide that data category<br><em>(useful when looking at percentages)</em>."
  },

  // These steps won't play if FastQC not present
  {
    element: "#fastqc_quality_plot_wrapper",
    placement: 'top',
    title: "View Originals",
    content: "Clicking a data point in some plots will show the original data. Click the original plot again to get back to the overview. (Have a go!)",
  },

  {
    element: ".hc-plot-handle:first",
    placement: 'top',
    title: "Resize Plots",
    content: "Drag the grey bar below plots to change their height.",
    backdropPadding: 10,
    onHide: function (tour) { $('.mqc-toolbox').css('z-index', 1200); }
  },

  {
    element: ".mqc-toolbox-buttons a[href=#mqc_cols]",
    placement: 'left',
    title: "MultiQC Toolbox",
    content: "Click one of the icons on the right to open the Toolbox",
  },
  {
    element: ".mqc-toolbox-buttons a[href=#mqc_cols]",
    placement: 'left',
    title: "Highlight Samples",
    content: "This tool allows you to highlight samples across plots and tables. Regexes mode allows for powerful pattern matching.",
    onShow: function (tour) { mqc_toolbox_openclose('#mqc_cols', true); },
  },
  {
    element: ".mqc-toolbox-buttons a[href=#mqc_renamesamples]",
    placement: 'left',
    title: "Rename Samples",
    content: "Here, you can rename samples, for example, cleaning up common suffixes. Data from excel can be pasted in for bulk renaming.",
    onShow: function (tour) { mqc_toolbox_openclose('#mqc_renamesamples', true); },
  },
  {
    element: ".mqc-toolbox-buttons a[href=#mqc_hidesamples]",
    placement: 'left',
    title: "Hide Samples",
    content: "This tool allows you to temporarily hide samples in the report",
    onShow: function (tour) { mqc_toolbox_openclose('#mqc_hidesamples', true); },
  },
  {
    element: ".mqc-toolbox-buttons a[href=#mqc_saveconfig]",
    placement: 'left',
    title: "Save Config",
    content: "By default, every toolbox setting is saved for the report. You can also share your set up with others and set defaults for all reports.",
    onShow: function(tour) { mqc_toolbox_openclose('#mqc_saveconfig', true); },
    onHide: function (tour) {
      mqc_toolbox_openclose('#mqc_saveconfig', false);
      $('.mqc-toolbox').css('z-index', 0);
    },
  },
  {
    orphan: true,
    title: "End of Tour",
    content: 'That\'s it for this tour - for more info, see the homepage: <a href="https://github.com/ewels/MultiQC" target="_blank">https://github.com/ewels/MultiQC</a>',
    onHide: function (tour) { $('.mqc-toolbox').css('z-index', 1200); }
  }];