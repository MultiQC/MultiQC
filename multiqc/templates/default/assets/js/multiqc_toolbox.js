////////////////////////////////////////////////
// MultiQC Report Toolbox Code
////////////////////////////////////////////////

var mqc_colours_idx = 0;
var mqc_colours = chroma.brewer.Set1;

//////////////////////////////////////////////////////
// TOOLBOX LISTENERS
//////////////////////////////////////////////////////
$(function () {
  
  // Listener to plot graphs when config loaded
  $(document).on('mqc_config_loaded', function(e){
    $('.hc-plot').each(function(){
      var target = $(this).attr('id');
      plot_graph(target, undefined, num_datasets_plot_limit);
    });
  });
  
  // Load any saved configuration
  load_mqc_config();

  // Toolbox buttons
  $('.mqc-toolbox-buttons a').click(function(e){
    e.preventDefault();
    var target = $(this).attr('href');
    mqc_toolbox_openclose(target);
  });

  // Highlight colour filters
  $('#mqc_color_form').submit(function(e){
    e.preventDefault();
    var f_text = $('#mqc_colour_filter').val().trim();
    var f_col = $('#mqc_colour_filter_color').val().trim();
    $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_colours_idx)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    apply_mqc_highlights();
    $('#mqc_colour_filter').val('');
    mqc_colours_idx += 1;
    if(mqc_colours_idx >= mqc_colours.length){ mqc_colours_idx = 0; }
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
  });
  
  // Sort general stats by highlight
  $('#mqc_genstat_sort_highlight').click(function(e){
    e.preventDefault();
    // collect highlighted rows
    var hrows = $('#general_stats_table tbody th.highlighted').parent().detach();
    hrows = hrows.sort(function (a, b) {
      return $(a).find('th').data('highlight') - $(b).find('th').data('highlight');
    });
    if($(this).data('direction') == 'desc'){
      hrows = hrows.get().reverse();
      $('#general_stats_table tbody').prepend(hrows);
      $(this).data('direction', 'asc');
    } else {
      $('#general_stats_table tbody').append(hrows);
      $(this).data('direction', 'desc');
    }
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
    $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_hidesamples_idx)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    apply_mqc_hidesamples();
    $('#mqc_hidesamples_filter').val('');
    mqc_hidesamples_idx += 1;
  });
  
  
  // Save config
  $('.mqc_saveconfig_btn').click(function(e){
    e.preventDefault();
    var tgt = $(this).data('target');
    var act = $(this).data('action');
    if(tgt == 'report'){
      if(act == 'save'){  mqc_save_config(undefined, 'local'); }
      if(act == 'clear'){ mqc_save_config(undefined, 'local', true); }
    }
    if(tgt == 'global'){
      if(act == 'save'){  mqc_save_config('global', 'local'); }
      if(act == 'clear'){ mqc_save_config('global', 'local', true); }
    }
    if(tgt == 'file'){
      mqc_save_config('global', 'file');
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
  // 'Enter' key pressed whilst editing a filter
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
    if($(this).parent().attr('id') == 'mqc_renamesamples'){ apply_mqc_renamesamples(); }
    if($(this).parent().attr('id') == 'mqc_hidesamples'){ apply_mqc_hidesamples(); }
  });
    
});



//////////////////////////////////////////////////////
// GENERAL TOOLBOX FUNCTIONS
//////////////////////////////////////////////////////
function mqc_toolbox_openclose (target, open){
  $('.mqc-toolbox-buttons li a').tooltip('hide');
  var btn = $('.mqc-toolbox-buttons li a[href="'+target+'"]');
  if(open === undefined){
    if(btn.hasClass('active')){ open = false; }
    else { open = true; }
  }
  var already_open = $('.mqc-toolbox').hasClass('active');
  if(open){
    $('.mqc-toolbox, .mqc-toolbox-buttons li a, .mqc_filter_section').removeClass('active');
    btn.addClass('active');
    $('.mqc-toolbox, '+target).addClass('active');
    $(document).trigger('mqc_toolbox_open');
    var timeout = already_open ? 0 : 510;
    setTimeout(function(){
      if(target == '#mqc_cols'){ $('#mqc_colour_filter').focus(); }
      if(target == '#mqc_renamesamples'){ $('#mqc_renamesamples_from').focus(); }
      if(target == '#mqc_hidesamples'){ $('#mqc_hidesamples_filter').focus(); }
    }, timeout);
  } else {
    btn.removeClass('active');
    $('.mqc-toolbox, .mqc-toolbox-buttons li a').removeClass('active');
    $(document).trigger('mqc_toolbox_close');
  }
}


//////////////////////////////////////////////////////
// HIGHLIGHT SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_highlights(){
  
  // Collect the filters into an array
  var f_texts = [];
  var f_cols = [];
  var regex_mode = false;
  if($('#mqc_cols .mqc_regex_mode').text() == 'Regex mode on'){
    regex_mode = true;
  }
  $('#mqc_col_filters li .f_text').each(function(){
    f_texts.push($(this).val());
    f_cols.push($(this).css('color'));
  });
  
  // Apply a 'background' highlight to remove default colouring first
  // Also highlight toolbox drawer icon
  if(f_texts.length > 0){
    f_texts.unshift('');
    f_cols.unshift('#cccccc');
    $('.mqc-toolbox-buttons a[href="#mqc_cols"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_cols"]').removeClass('in_use');
  }
  
  window.mqc_highlight_f_texts = f_texts;
  window.mqc_highlight_f_cols = f_cols;
  window.mqc_highlight_regex_mode = regex_mode;
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_highlights', [f_texts, f_cols, regex_mode]);
}

//////////////////////////////////////////////////////
// RENAME SAMPLES
//////////////////////////////////////////////////////

function apply_mqc_renamesamples(){
  
  // Collect filters
  var f_texts = [];
  var t_texts = [];
  var regex_mode = false;
  $('#mqc_renamesamples_filters .from_text').each(function(){ f_texts.push($(this).val()); });
  $('#mqc_renamesamples_filters .to_text').each(function(){ t_texts.push($(this).val()); });
  if($('#mqc_renamesamples .mqc_regex_mode').text() == 'Regex mode on'){ regex_mode = true; }
  
  // If something was renamed, highlight the toolbox icon
  if(f_texts.length > 0){
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').removeClass('in_use');
  }
  
  window.mqc_rename_f_texts = f_texts;
  window.mqc_rename_t_texts = t_texts;
  window.mqc_rename_regex_mode = regex_mode;
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_renamesamples', [f_texts, t_texts, regex_mode]);
}


//////////////////////////////////////////////////////
// HIDE SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_hidesamples(){
  // Collect the filters into an array
  var f_texts = [];
  var regex_mode = false;
  if($('#mqc_hidesamples .mqc_regex_mode').text() == 'Regex mode on'){
    regex_mode = true;
  }
  $('#mqc_hidesamples_filters li .f_text').each(function(){
    f_texts.push($(this).val());
  });
  
  // If something was renamed, highlight the toolbox icon
  if(f_texts.length > 0){
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').removeClass('in_use');
  }
  
  window.mqc_hide_f_texts = f_texts;
  window.mqc_hide_regex_mode = regex_mode;
  
  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_hidesamples', [f_texts, regex_mode]);
}



//////////////////////////////////////////////////////
// SAVE TOOLBOX SETTINGS
//////////////////////////////////////////////////////

// Save the current configuration setup
function mqc_save_config(target, method, clear){
  if(target === undefined){ target = report_id }
  if(method === undefined){ method = 'local'; }
  var config = {};
  
  // Collect the toolbox vars
  config['highlights_f_texts'] =  window.mqc_highlight_f_texts;
  config['highlights_f_cols'] =   window.mqc_highlight_f_cols;
  config['highlight_regex'] =     window.mqc_highlight_regex_mode;
  config['rename_from_texts'] =   window.mqc_rename_f_texts;
  config['rename_to_texts'] =     window.mqc_rename_t_texts;
  config['rename_regex'] =        window.mqc_rename_regex_mode;
  config['hidesamples_f_texts'] = window.mqc_hide_f_texts;
  config['hidesamples_regex'] =   window.mqc_hide_regex_mode;
  
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
      prev_config[target] = {};
    } else {
      prev_config[target] = config;
      prev_config[target]['last_updated'] = Date();
    }
    localStorage.setItem("mqc_config", JSON.stringify(prev_config));
    // Success message
    var insAfter = '.mqc-save-config-report-section p:first-child';
    if(target == 'global'){ insAfter = '.mqc-save-config-global-section p:first-child'; }
    $('<p class="text-success" id="mqc-save-'+target+'-success">Config saved.</p>').hide().insertAfter($(insAfter)).slideDown(function(){
      setTimeout(function(){
        $('#mqc-save-'+target+'-success').slideUp(function(){ $(this).remove(); });
      }, 5000);
    });
  }
  if(method == 'file'){
    var f = "// Config file for MultiQC\n// http://multiqc.info\n";
    f += "// Generated "+Date()+"\n// Original report path: "+$('#mqc_output_path').text()+"\n";
    f += "\nmqc_config_file_cfg = "+JSON.stringify(config, null, '  ');
    var fblob = new Blob([f], {type: "application/javascript;charset=utf-8"});
    saveAs(fblob, "multiqc_config.js");
    setTimeout(function(){
      alert('Now move multiqc_config.js from your downloads folder to the directory where this report is located.');
    }, 500);
  }
}

//////////////////////////////////////////////////////
// LOAD TOOLBOX SETTINGS
//////////////////////////////////////////////////////
function load_mqc_config(){
  var config = {};
  // Get local configs - global and then this path
  try {
    var local_config = localStorage.getItem("mqc_config");
    if(local_config !== null && local_config !== undefined){
      local_config = JSON.parse(local_config);
      if(local_config['global'] !== undefined){ console.log('Loaded local global config'); }
      for (var attr in local_config['global']) {
        config[attr] = local_config['global'][attr];
      }
      if(local_config[report_id] !== undefined){ console.log('Loaded local report config'); }
      for (var attr in local_config[report_id]) {
        config[attr] = local_config[report_id][attr];
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
  
  // Apply config - highlights
  if(notEmptyObj(config['highlight_regex'])){
    if(config['highlight_regex'] == true){
      $('#mqc_cols .mqc_regex_mode span').removeClass('off').addClass('on').text('on');
      window.mqc_highlight_regex_mode = true;
    }
  }
  if(notEmptyObj(config['highlights_f_texts']) && notEmptyObj(config['highlights_f_cols'])){
    $.each(config['highlights_f_texts'], function(idx, f_text){
      var f_col = config['highlights_f_cols'][idx];
      $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
      window.mqc_highlight_f_texts.push(f_text);
      window.mqc_highlight_f_cols.push(f_col);
      mqc_colours_idx += 1;
    });
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
    $(document).trigger('mqc_highlights', [config['highlights_f_texts'], config['highlights_f_cols'], config['highlight_regex']]);
  }
  
  // Rename samples
  if(notEmptyObj(config['rename_regex'])){
    if(config['rename_regex'] == true){
      $('#mqc_renamesamples .mqc_regex_mode span').removeClass('off').addClass('on').text('on');
      window.mqc_rename_regex_mode = true;
    }
  }
  if(notEmptyObj(config['rename_from_texts']) && notEmptyObj(config['rename_to_texts'])){
    $.each(config['rename_from_texts'], function(idx, from_text){
      var to_text = config['rename_to_texts'][idx];
      if(from_text.length == 0){ return true; }
      var li = '<li><input class="f_text from_text" value="'+from_text+'" />'
      li += '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="'+to_text+'" />'
      li += '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>'
      window.mqc_rename_f_texts.push(from_text);
      window.mqc_rename_t_texts.push(to_text);
      $('#mqc_renamesamples_filters').append(li);
    });
    $(document).trigger('mqc_renamesamples', [config['rename_from_texts'], config['rename_from_texts'], config['rename_regex']]);
  }
  
  // Hide samples
  if(notEmptyObj(config['hidesamples_regex'])){
    if(config['hidesamples_regex'] == true){
      $('#mqc_hidesamples .mqc_regex_mode span').removeClass('off').addClass('on').text('on');
      window.mqc_hide_regex_mode = true;
    }
  }
  if(notEmptyObj(config['hidesamples_f_texts'])){
    $.each(config['hidesamples_f_texts'], function(idx, f_text){
      if(f_text.length == 0){ return true; }
      $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
      window.mqc_hide_f_texts.push(f_text);
    });
    $(document).trigger('mqc_hidesamples', [config['hidesamples_f_texts'], config['hidesamples_regex']]);
  }

  // Trigger loaded event to initialise plots
  $(document).trigger('mqc_config_loaded');
  
}



