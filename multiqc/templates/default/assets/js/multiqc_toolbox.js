////////////////////////////////////////////////
// MultiQC Report Toolbox Code
////////////////////////////////////////////////

var mqc_colours_idx = 0;
var mqc_colours = chroma.brewer.Set1;

//////////////////////////////////////////////////////
// TOOLBOX LISTENERS
//////////////////////////////////////////////////////
$(function () {
  
  // Load any saved configuration
  load_mqc_config();

  // Toolbox buttons
  $('.mqc-toolbox-buttons a').click(function(e){
    e.preventDefault();
    var target = $(this).attr('href');
    if(target == '#mqc_cleartoolbox'){
      // Empty the toolbox filter lists
      $('#mqc_col_filters').empty();
      $('#mqc_renamesamples_filters').empty();
      $('#mqc_hidesamples_filters').empty();
      // Reset the highlight colours
      mqc_colours_idx = 0;
      $('#mqc_colour_filter_color').val(mqc_colours[0]);
      // Close the toolbox
      mqc_toolbox_openclose ('#mqc_col_filters', false);
      // Apply new set of filters (reset page)
      apply_mqc_highlights();
      apply_mqc_renamesamples();
      apply_mqc_hidesamples();
    } else {
      mqc_toolbox_openclose(target);
    }
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
    var already_exists = false;
    $('#mqc_col_filters li .f_text').each(function(){
      if($(this).val() == f_text){
        alert('Error - highlight text "'+f_text+'" already exists');
        already_exists = true;
        return false;
      }
    });
    if(already_exists) { return false; }
    $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_colours_idx)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    apply_mqc_highlights();
    $('#mqc_colour_filter').val('');
    mqc_colours_idx += 1;
    if(mqc_colours_idx >= mqc_colours.length){ mqc_colours_idx = 0; }
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
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
    if(target == 'mqc_renamesamples_filters'){
      apply_mqc_renamesamples();
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
    if(target == 'mqc_renamesamples_filters'){ apply_mqc_renamesamples(); }
  });
  // Use jQuery UI to make the colour filters sortable
  $("#mqc_col_filters").sortable();
  $("#mqc_col_filters").on("sortstop", function(event, ui){
    apply_mqc_highlights();
  });
  // Regex mode text
  $('.mqc_regex_mode').click(function(){
    var rswitch = $(this).find('span');
    if(rswitch.text() == 'off'){
      rswitch.removeClass('off').addClass('on').text('on');
    } else {
      rswitch.removeClass('on').addClass('off').text('off');
    }
    if($(this).parent().attr('id') == 'mqc_cols'){ apply_mqc_highlights(); }
    if($(this).parent().attr('id') == 'mqc_hidesamples'){ apply_mqc_hidesamples(); }
  });
  
  
});



//////////////////////////////////////////////////////
// GENERAL TOOLBOX FUNCTIONS
//////////////////////////////////////////////////////
function mqc_toolbox_openclose (target, open){
  $('.mqc-toolbox-buttons a').tooltip('hide');
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


//////////////////////////////////////////////////////
// HIGHLIGHT SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_highlights(){
  
  // Collect the filters into an array
  var f_matches = 0;
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
              f_matches += 1;
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
              f_matches += 1;
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
        f_matches += 1;
      }
    });
    $(this).css('color', thiscol);
  });
  
  // If something was renamed, highlight the toolbox icon
  if(f_matches > 0){
    $('.mqc-toolbox-buttons a[href="#mqc_cols"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_cols"]').removeClass('in_use');
  }
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_highlights', [f_texts, f_cols, regex_mode]);
  mqc_autosave();
}

//////////////////////////////////////////////////////
// RENAME SAMPLES
//////////////////////////////////////////////////////

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
  
  // If something was renamed, highlight the toolbox icon
  if($('#mqc_renamesamples_filters .from_text').length > 0){
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').removeClass('in_use');
  }
  
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



//////////////////////////////////////////////////////
// HIDE SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_hidesamples(){
  // Collect the filters into an array
  var f_matches = 0;
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
        var num_hidden = 0;
        $.each($(this).highcharts().series, function(j, s){
          var match = false;
          $.each(f_texts, function(idx, f_text){
            if((regex_mode && s.name.match(f_text)) || (!regex_mode && s.name.indexOf(f_text) > -1)){
              match = true;
              f_matches += 1;
            }
          });
          if (s.visible && match) { s.hide(); }
          if (!s.visible && !match) { s.show(); }
          if (!s.visible){ num_hidden += 1; }
        });
        // Reset warnings and plot hiding.
        $(this).parent().parent().find('.samples-hidden-warning').remove();
        $(this).parent().parent().find('.hc-plot-wrapper, .btn-group').show();
        // Some series hidden. Show a warning text string.
        if(num_hidden > 0) {
          var alert = '<div class="samples-hidden-warning alert alert-danger"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+num_hidden+' samples hidden in toolbox. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
          if($(this).parent().prev().hasClass('btn-group')){
            $(this).parent().prev().before(alert);
          } else {
            $(this).parent().before(alert);
          }
        }
        // All series hidden. Hide the graph.
        if(num_hidden == $(this).highcharts().series.length){
          $(this).parent().parent().find('.hc-plot-wrapper, .btn-group').hide();
        }
        
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
              f_matches += 1;
              replot.xAxis.categories.splice(idx, 1);
              $.each(replot.series, function(j, s){
                replot.series[j].data.splice(idx, 1);
              });
              return false;
            }
          });
        };
        
        // Reset warnings and plot hiding.
        $(this).parent().parent().find('.samples-hidden-warning').remove();
        $(this).parent().parent().find('.hc-plot-wrapper, .btn-group').show();
        
        // Replot graph
        highcharts_plots['#'+plotid] = new Highcharts.Chart(replot);
        
        // Some series hidden. Show a warning text string.
        if(matches > 0) {
          var alert = '<div class="samples-hidden-warning alert alert-danger"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+matches+' samples hidden in toolbox. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>';
          if($(this).parent().prev().hasClass('btn-group')){
            $(this).parent().prev().before(alert);
          } else {
            $(this).parent().before(alert);
          }
        }
        // All series hidden. Hide the graph.
        if(matches == highcharts_plot_options['#'+plotid].xAxis.categories.length){
          $(this).parent().parent().find('.hc-plot-wrapper, .btn-group').hide();
        }
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
        f_matches += 1;
      }
    });
    if(match){ $(this).parent().hide(); }
    else { $(this).parent().show(); }
  });
  // Hide empty columns
  var gsthidx = 0;
  $("#general_stats_table thead th, #general_stats_table tbody tr td").show();
  $("#general_stats_table thead th").each(function(){
    if(gsthidx == 0){ gsthidx += 1; return true; }
    var count = 0;
    var empties = 0;
    $("#general_stats_table tbody tr td:nth-child("+(gsthidx+2)+")").filter(":visible").each(function(){
      count += 1;
      if($(this).text() == ''){ empties += 1; }
    });
    if(count > 0 && count == empties){
      $(this).hide();
      $("#general_stats_table tbody tr td:nth-child("+(gsthidx+2)+")").hide();
    }
    gsthidx += 1;
  });
  
  // If something was renamed, highlight the toolbox icon
  if(f_matches > 0){
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').removeClass('in_use');
  }
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_hidesamples', [f_texts, regex_mode]);
  mqc_autosave();
}



//////////////////////////////////////////////////////
// SAVE TOOLBOX SETTINGS
//////////////////////////////////////////////////////


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
      $('#mqc_cols .mqc_regex_mode span').removeClass('on').addClass('off').text('off');
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
      $('#mqc_cols .mqc_regex_mode span').removeClass('on').addClass('off').text('off');
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



