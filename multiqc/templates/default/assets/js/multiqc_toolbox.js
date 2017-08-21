////////////////////////////////////////////////
// MultiQC Report Toolbox Code
////////////////////////////////////////////////

var mqc_colours_idx = 0;
var mqc_colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999'];

//////////////////////////////////////////////////////
// TOOLBOX LISTENERS
//////////////////////////////////////////////////////
$(function () {

  // Batch sample renaming buttons
  $('.mqc_sname_switches').click(function(e){
    e.preventDefault();
    if($(this).hasClass('active')){
      return false;
    }
    $('#mqc_sname_switches button').removeClass('active');
    $(this).addClass('active');
    // Clear previous bulk renaming entries
    $('.mqc_sname_switches_li').remove();
    // Build new renaming entries and apply
    var j = $(this).data('index');
    if(j == 0){
      apply_mqc_renamesamples();
    } else {
      for(i=0; i<mqc_sample_names_rename.length; i++){
        var ft = mqc_sample_names_rename[i][0];
        var tt = mqc_sample_names_rename[i][j];
        $('#mqc_renamesamples_filters').append('<li class="mqc_sname_switches_li"> \
          <input class="f_text from_text" value="'+ft+'" /><small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="'+tt+'" /> \
          <button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button> \
        </li>');
      }
      apply_mqc_renamesamples();
    }
  });

  // Hide toolbox when clicking outside
  $(document).mouseup(function (e){
    if (!$(".mqc-toolbox").is(e.target) && $(".mqc-toolbox").has(e.target).length === 0){
      if($('.mqc-toolbox').hasClass('active')){
        mqc_toolbox_openclose(undefined, false);
      }
    }
  });

  // Hide toolbox when a modal is shown
  $('.modal').on('show.bs.modal', function(e){
    if($('.mqc-toolbox').hasClass('active')){
      mqc_toolbox_openclose(undefined, false);
    }
  });

  // Listener to re-plot graphs if config loaded
  $(document).on('mqc_config_loaded', function(e){
    $('.hc-plot').each(function(){
      var target = $(this).attr('id');
      plot_graph(target, undefined, num_datasets_plot_limit);
    });
  });

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
    $('#mqc_col_filters').append('<li style="color:'+f_col+';" id="'+hashCode(f_text+f_col)+'"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" tabindex="'+(mqc_colours_idx)+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    $('#mqc_cols_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
    $('#mqc_colour_filter').val('');
    mqc_colours_idx += 1;
    if(mqc_colours_idx >= mqc_colours.length){ mqc_colours_idx = 0; }
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
  });
  $('#mqc_cols_apply').click(function(e){
    if(apply_mqc_highlights()){
      $(this).attr('disabled', true).removeClass('btn-primary').addClass('btn-default');
    }
  });

  // Rename samples
  var mqc_renamesamples_idx = 300;
  $('#mqc_renamesamples_form').submit(function (event) {
    event.preventDefault();

    var from_text = $('#mqc_renamesamples_from').val().trim();
    var to_text = $('#mqc_renamesamples_to').val().trim();

    if (from_text.length == 0) {
      alert('Error - "From" text must not be blank.');
      return false;
    }

    var li = '<li><input class="f_text from_text" value="'+from_text+'" tabindex="'+(mqc_renamesamples_idx)+'" />'
    li += '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="'+to_text+'" tabindex="'+(mqc_renamesamples_idx+1)+'" />'
    li += '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>'
    $('#mqc_renamesamples_filters').append(li);
    $('#mqc_rename_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');

    // Reset form
    $('#mqc_renamesamples_from').val('');
    $('#mqc_renamesamples_to').val('');
    mqc_renamesamples_idx += 2;
    $('#mqc_renamesamples_form input:first').focus();
  });

  $('#mqc_rename_apply').click(function(e){
    if(apply_mqc_renamesamples()){
      $(this).attr('disabled', true).removeClass('btn-primary').addClass('btn-default');
    }
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
    $('#mqc_rename_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
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
    $('#mqc_hide_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
    $('#mqc_hidesamples_filter').val('');
    mqc_hidesamples_idx += 1;
  });
  $('.mqc_hidesamples_showhide').change(function(e){
    $('#mqc_hide_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
  });
  $('#mqc_hide_apply').click(function(e){
    if(apply_mqc_hidesamples()){
      $(this).attr('disabled', true).removeClass('btn-primary').addClass('btn-default');
    }
  });

  // EXPORTING PLOTS
  // Change text on download button
  $('#mqc_exportplots a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
    if($(e.target).attr('href') == '#mqc_data_download'){
      $('#mqc-dl-plot-txt').text('Data');
    } else {
      $('#mqc-dl-plot-txt').text('Images');
    }
  });
  // Load the plot exporter
  if($('.hc-plot').length > 0){
    $('.hc-plot').each(function(){
      var fname = $(this).attr('id');
      $('#mqc_export_selectplots').append('<div class="checkbox"><label><input type="checkbox" value="'+fname+'" checked> '+fname+'</label></div>');
    });
    // Select all / none for checkboxes
    $('#mqc_export_sall').click(function(e){
      e.preventDefault();
      $('#mqc_export_selectplots input').prop('checked', true);
    });
    $('#mqc_export_snone').click(function(e){
      e.preventDefault();
      $('#mqc_export_selectplots input').prop('checked', false);
    });
    // Aspect ratio fixed
    var mqc_exp_aspect_ratio = $('#mqc_exp_width').val() / $('#mqc_exp_height').val();
    $('#mqc_export_aspratio').change(function(){
      if($(this).is(':checked')){
        mqc_exp_aspect_ratio = $('#mqc_exp_width').val() / $('#mqc_exp_height').val();
      }
    });
    $('#mqc_exp_width').keyup(function(){
      if($('#mqc_export_aspratio').is(':checked')){
        $('#mqc_exp_height').val( $(this).val() / mqc_exp_aspect_ratio );
      }
    });
    $('#mqc_exp_height').keyup(function(){
      if($('#mqc_export_aspratio').is(':checked')){
        $('#mqc_exp_width').val( $(this).val() * mqc_exp_aspect_ratio );
      }
    });

    // Export the plots
    $('#mqc_exportplots').submit(function(e){
      e.preventDefault();
      var skipped_plots = 0;
      ////// EXPORT PLOT IMAGES
      //////
      if($('#mqc_image_download').is(':visible')){
        var ft = $('#mqc_export_ft').val();
        var f_scale = parseInt($('#mqc_export_scaling').val());
        var f_width = parseInt($('#mqc_exp_width').val()) / f_scale;
        var f_height = parseInt($('#mqc_exp_height').val()) / f_scale;
        $('#mqc_export_selectplots input:checked').each(function(){
          var fname = $(this).val();
          var hc = $('#'+fname).highcharts();
          var cfg = {
            type: ft,
            filename: fname,
            sourceWidth: f_width,
            sourceHeight: f_height,
            scale: f_scale
          };
          if(hc !== undefined){
            hc.exportChartLocal(cfg);
          } else if($('#'+fname).hasClass('has-custom-export')){
            $('#'+fname).trigger('mqc_plotexport_image', cfg);
          } else {
            skipped_plots += 1;
          }
        });
        if(skipped_plots > 0){
          alert("Warning: "+skipped_plots+" plots skipped.\n\nNote that it is not currently possible to export dot plot images from reports. Data exports do work.");
        }
      }
      ////// EXPORT PLOT DATA
      //////
      else if($('#mqc_data_download').is(':visible')){
        $('#mqc_export_selectplots input:checked').each(function(){
          try {
            var target = $(this).val();
            var ft = $('#mqc_export_data_ft').val();
            var fname = target+'.'+ft;
            var sep = ft == 'tsv' ? "\t" : ',';
            // Custom plot not in mqc_plots
            if (mqc_plots[target] == undefined){
              if($('#'+target).hasClass('has-custom-export')){
                $('#'+target).trigger('mqc_plotexport_data', { 'target': target, 'ft': ft, 'fname': fname, 'sep': sep });
              } else {
                skipped_plots += 1;
              }
            }
            // If JSON then just dump everything
            else if(ft == 'json'){
              json_str = JSON.stringify(mqc_plots[target], null, 2);
              var blob = new Blob([json_str], {type: "text/plain;charset=utf-8"});
              saveAs(blob, fname);
            }
            // Beeswarm plots must be done manually
            else if(mqc_plots[target]['plot_type'] == 'beeswarm'){
              // Header line
              datastring = 'Sample';
              for(var j=0; j<mqc_plots[target]['categories'].length; j++){
                datastring += sep+mqc_plots[target]['categories'][j]['description'];
              }
              datastring += "\n";
              // This assumes that the same samples are in all rows
              // TODO: Check and throw error if this isn't the case
              var rows = Array();
              for(var j=0; j<mqc_plots[target]['samples'][0].length; j++){
                rows[j]=Array(mqc_plots[target]['samples'][0][j]);
              }
              for(var j=0; j<mqc_plots[target]['datasets'].length; j++){
                for(var k=0; k<mqc_plots[target]['datasets'][j].length; k++){
                  rows[k].push(mqc_plots[target]['datasets'][j][k]);
                }
              }
              for(var j=0; j<rows.length; j++){
                datastring += rows[j].join(sep)+"\n";
              }
              var blob = new Blob([datastring], {type: "text/plain;charset=utf-8"});
              saveAs(blob, fname);
            }
            // Normal plot - use HighCharts plugin to get the data from the plot
            else if(ft == 'tsv' || ft == 'csv'){
              var hc = $('#'+target).highcharts();
              if(hc !== undefined){
                hc.update({ exporting: { csv: { itemDelimiter: sep } } });
                var blob = new Blob([hc.getCSV()], {type: "text/plain;charset=utf-8"});
                saveAs(blob, fname);
              } else {
                skipped_plots += 1;
              }
            } else {
              skipped_plots += 1;
            }
          } catch(e){
            console.error(e);
            skipped_plots += 1;
          }
        });
        if(skipped_plots > 0){
          alert("Warning: Could not export data from "+skipped_plots+" plots.");
        }
      } else { alert("Error - don't know what to export!"); }
    });
  } else {
    $('#mqc_exportplots').hide();
    $('.mqc-toolbox-buttons a[href=#mqc_exportplots]').parent().hide();
  }

  /// SAVING STUFF
  // Load the saved setting names
  populate_mqc_saveselect();
  // Save config
  $('#mqc_saveconfig_form').submit(function(e){
    e.preventDefault();
    var name = $(this).find('input').val().trim();
    if(name == ''){
      alert('Error - you must name the saved settings.');
    } else {
      mqc_save_config(name);
    }
  });
  // Load config
  $('#mqc_loadconfig_form').submit(function(e){
    e.preventDefault();
    var name = $(this).find('select').val().trim();
    if(name == ''){
      alert('Error - No saved setting selected.');
    } else {
      load_mqc_config(name);
    }
  });
  // Delete config
  $('.mqc_config_clear').click(function(e){
    e.preventDefault();
    var name = $('#mqc_loadconfig_form select').val().trim();
    if(name == ''){
      alert('Error - no saved settings selected.');
    } else {
      if(confirm("Delete saved settings '"+name+"'?")){
        mqc_save_config(name, true);
      }
    }
  });
  // Set current config as default
  $('.mqc_config_set_default').click(function(e){
    e.preventDefault();
    var name = $('#mqc_loadconfig_form select').val().trim();
    if(name == ''){
      alert('Error - no saved settings selected.');
    } else {
        load_mqc_config(name);
        mqc_save_config(name, false, true);
    }
  });
    // Clear current config default
  $('.mqc_config_clear_default').click(function(e){
    e.preventDefault();
    mqc_clear_default_config();
  });

  // Filter text is changed
  $('.mqc_filters').on('blur', 'li input', function(){
    var target = $(this).parent().parent().attr('id');
    if(target == 'mqc_col_filters'){
      $('#mqc_cols_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
    }
    if(target == 'mqc_renamesamples_filters'){
      $('#mqc_rename_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
    }
    if(target == 'mqc_hidesamples_filters'){
      $('#mqc_hide_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
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
    if(target == 'mqc_col_filters'){ $('#mqc_cols_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary'); }
    if(target == 'mqc_hidesamples_filters'){ $('#mqc_hide_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary'); }
    if(target == 'mqc_renamesamples_filters'){ $('#mqc_rename_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary'); }
  });
  // Clear all filters button
  $('.mqc_toolbox_clear').click(function(){
    var target = $(this).closest('.mqc_filter_section').find('.mqc_filters').attr('id');
    $('#'+target).empty();
    if(target == 'mqc_col_filters'){ $('#mqc_cols_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary'); }
    if(target == 'mqc_hidesamples_filters'){ $('#mqc_hide_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary'); }
    if(target == 'mqc_renamesamples_filters'){ $('#mqc_rename_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary'); }
  });

  // Use jQuery UI to make the colour filters sortable
  $("#mqc_col_filters").sortable();
  $("#mqc_col_filters").on("sortstop", function(event, ui){
    $('#mqc_cols_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
  });
  // Regex mode text
  $('.mqc_regex_mode').click(function(){
    var rswitch = $(this).find('.re_mode');
    if(rswitch.text() == 'off'){
      rswitch.removeClass('off').addClass('on').text('on');
    } else {
      rswitch.removeClass('on').addClass('off').text('off');
    }
    if($(this).data('target') == 'mqc_cols'){
      $('#mqc_cols_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
    }
    if($(this).data('target') == 'mqc_renamesamples'){
      $('#mqc_rename_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
    }
    if($(this).data('target') == 'mqc_hidesamples'){
      $('#mqc_hide_apply').attr('disabled', false).removeClass('btn-default').addClass('btn-primary');
    }
  });

  /////////////////////////
  // REGEX HELP MODAL
  /////////////////////////
  $('.regex_example_buttons button').click(function(e){
    e.preventDefault();
    $('.regex_example_demo input').val( $(this).data('example') );
    regex_example_test();
  });
  $('.regex_example_demo input').keyup(function(e){
    regex_example_test();
  });
  function regex_example_test(){
    var re = $('.regex_example_demo input').val();
    console.log('Testing '+re);
    $('.regex_example_demo pre span').each(function(){
      $(this).removeClass();
      if( $(this).text().match(re) ){
        console.log('Matches '+$(this).text());
        $(this).addClass('mark text-success');
      } else {
        console.log('Matches '+$(this).text());
        $(this).addClass('text-muted');
      }
    });
  }

});

//////////////////////////////////////////////////////
// UTILITY FUNCTIONS
//////////////////////////////////////////////////////
function hashCode(str) {
  var hash = 0;
  if (str.length == 0) return hash;
  for (i = 0; i < str.length; i++) {
      char = str.charCodeAt(i);
      hash = ((hash<<5)-hash)+char;
      hash = hash & hash; // Convert to 32bit integer
  }
  return hash;
}


//////////////////////////////////////////////////////
// GENERAL TOOLBOX FUNCTIONS
//////////////////////////////////////////////////////
function mqc_toolbox_openclose (target, open){
  // Hide any open tooltip so it's not left dangling
  $('.mqc-toolbox-buttons li a').tooltip('hide');
  // Find if what we clicked is already open
  var btn = $('.mqc-toolbox-buttons li a[href="'+target+'"]');
  if(open === undefined){
    if(btn.hasClass('active')){ open = false; }
    else { open = true; }
  }
  var already_open = $('.mqc-toolbox').hasClass('active');
  if(open){
    if(already_open){
      mqc_toolbox_confirmapply();
    }
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
    mqc_toolbox_confirmapply();
    btn.removeClass('active');
    $('.mqc-toolbox, .mqc-toolbox-buttons li a').removeClass('active');
    $(document).trigger('mqc_toolbox_close');
  }
}
function mqc_toolbox_confirmapply(){
  // Check if there's anything waiting to be applied
  if($('#mqc_cols_apply').is(':enabled') && $('#mqc_cols').is(':visible')){
    $.toast({
      heading: 'Highlights Not Applied',
      text: "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      icon: 'warning',
      position: 'bottom-right',
      hideAfter: 5000
    });
  }
  if($('#mqc_rename_apply').is(':enabled') && $('#mqc_renamesamples').is(':visible')){
    $.toast({
      heading: 'Rename Patterns Not Applied',
      text: "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      icon: 'warning',
      position: 'bottom-right',
      hideAfter: 5000
    });
  }
  if($('#mqc_hide_apply').is(':enabled') && $('#mqc_hidesamples').is(':visible')){
    $.toast({
      heading: 'Hide Samples Not Applied',
      text: "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      icon: 'warning',
      position: 'bottom-right',
      hideAfter: 5000
    });
  }
}

function validate_regexp(pattern) {
  try {
    new RegExp(pattern, 'g');
    return true
  } catch (error) {
    $.toast({
      heading: 'Invalid Regular Expression!',
      text: "Apologies, your regular expression pattern is invalid: <code>"+pattern+"</code><br><br>"+
            'For more help and testing, try it out at <a href="https://regex101.com/" target="_blank">regex101.com</a>.',
      icon: 'error',
      position: 'bottom-right',
      hideAfter: 5000
    });
    return false
  }
}

//////////////////////////////////////////////////////
// HIGHLIGHT SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_highlights(){

  // Collect the filters into an array
  var f_texts = [];
  var f_cols = [];
  var regex_mode = $('#mqc_cols .mqc_regex_mode .re_mode').hasClass('on');
  var num_errors = 0;
  $('#mqc_col_filters li').each(function() {
    var inputElement = $(this).find('.f_text');
    var pattern = inputElement.val();
    // Validate RegExp
    $(this).removeClass('bg-danger');
    if (regex_mode && !validate_regexp(pattern)) {
      $(this).addClass('bg-danger');
      num_errors++;
    }

    // Only add pattern if it hasn't already been added
    if (f_texts.indexOf(pattern) < 0) {
      f_texts.push(pattern);
      f_cols.push(inputElement.css('color'));
    } else {
      f_cols[f_texts.indexOf(pattern)] = inputElement.css('color');
    }
  });
  if(num_errors > 0){
    return false;
  }

  // Apply a 'background' highlight to remove default colouring first
  // Also highlight toolbox drawer icon
  if (f_texts.length > 0 && f_texts.indexOf('') < 0) {
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

  return true;
}

//////////////////////////////////////////////////////
// RENAME SAMPLES
//////////////////////////////////////////////////////

function apply_mqc_renamesamples () {
  var valid_from_texts = [];
  var valid_to_texts = [];
  var regex_mode = $('#mqc_renamesamples .mqc_regex_mode .re_mode').hasClass('on');
  var num_errors = 0;
  // Collect filters
  var f_texts = $('#mqc_renamesamples_filters > li').each(function () {
    var from_text = $(this).find(".from_text").val();
    var to_text = $(this).find(".to_text").val();
    // Validate RegExp
    $(this).removeClass('bg-danger')
    if (regex_mode && !validate_regexp(from_text)) {
      $(this).addClass('bg-danger');
      num_errors++;
    }
    valid_from_texts.push(from_text);
    valid_to_texts.push(to_text);
  });
  if(num_errors > 0){
    return false;
  }

  // If something was renamed, highlight the toolbox icon
  if (valid_from_texts.length > 0) {
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').removeClass('in_use');
  }

  window.mqc_rename_f_texts = valid_from_texts;
  window.mqc_rename_t_texts = valid_to_texts;
  window.mqc_rename_regex_mode = regex_mode;

  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_renamesamples', [
    window.mqc_rename_f_texts,
    window.mqc_rename_t_texts,
    regex_mode
  ]);

  return true;
}

//////////////////////////////////////////////////////
// HIDE SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_hidesamples(){
  // Collect the filters into an array
  var mode = $('.mqc_hidesamples_showhide:checked').val() == 'show' ? 'show' : 'hide';
  var regex_mode = $('#mqc_hidesamples .mqc_regex_mode .re_mode').hasClass('on');
  var f_texts = [];
  var num_errors = 0;
  $('#mqc_hidesamples_filters li').each(function() {
    var pattern = $(this).find('.f_text').val();
    // Validate RegExp
    $(this).removeClass('bg-danger')
    if (regex_mode && !validate_regexp(pattern)) {
      $(this).addClass('bg-danger');
      num_errors++;
    }
    f_texts.push(pattern);
  });
  if(num_errors > 0){
    return false;
  }

  // If something was hidden, highlight the toolbox icon
  if(f_texts.length > 0){
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').addClass('in_use');
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').removeClass('in_use');
  }

  window.mqc_hide_mode = mode;
  window.mqc_hide_f_texts = f_texts;
  window.mqc_hide_regex_mode = regex_mode;

  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger('mqc_hidesamples', [f_texts, regex_mode]);

  return true;
}

//////////////////////////////////////////////////////
// SAVE TOOLBOX SETTINGS
//////////////////////////////////////////////////////

// Save the current configuration setup
function mqc_save_config(name, clear, as_default){
  if(name === undefined){ return false; }
  var config = {};

  // Collect the toolbox vars
  config['highlights_f_texts'] =  window.mqc_highlight_f_texts;
  config['highlights_f_cols'] =   window.mqc_highlight_f_cols;
  config['highlight_regex'] =     window.mqc_highlight_regex_mode;
  config['rename_from_texts'] =   window.mqc_rename_f_texts;
  config['rename_to_texts'] =     window.mqc_rename_t_texts;
  config['rename_regex'] =        window.mqc_rename_regex_mode;
  config['hidesamples_mode'] =    window.mqc_hide_mode;
  config['hidesamples_f_texts'] = window.mqc_hide_f_texts;
  config['hidesamples_regex'] =   window.mqc_hide_regex_mode;

  var prev_config = {};
  // Load existing configs (inc. from other reports)
  try {
    try {

      prev_config = localStorage.getItem("mqc_config");
      if(prev_config !== null && prev_config !== undefined){
        prev_config = JSON.parse(prev_config);
      } else {
        prev_config = {};
      }

      // Update config obj with current config
      if(clear == true){
        delete prev_config[name];
      } else {
        prev_config[name] = config;
        prev_config[name]['last_updated'] = Date();
        if (as_default) {
          for (var c in prev_config) {
            if (prev_config.hasOwnProperty(c)) {
              prev_config[c]['default'] = false;
            }
          }
        }
        prev_config[name]['default'] = as_default;
        if (as_default) console.log('Set new default config!');
      }
      localStorage.setItem("mqc_config", JSON.stringify(prev_config));

    } catch(e){
      console.log('Could not access localStorage');
    }

    if(clear == true){
      // Remove from load select box
      $("#mqc_loadconfig_form select option:contains('"+name+"')").remove();
      // Successfully deleted message
      $('<p class="text-danger" id="mqc-cleared-success">Settings deleted.</p>').hide().insertBefore($('#mqc_loadconfig_form .actions')).slideDown(function(){
        setTimeout(function(){
          $('#mqc-cleared-success').slideUp(function(){ $(this).remove(); });
        }, 5000);
      });
    } else {
      // Remove from load select box
      $("#mqc_loadconfig_form select option:contains('"+name+"')").remove();
      // Add new name to load select box and select it
      $('#mqc_loadconfig_form select').prepend('<option>'+name+(as_default?' [default]':'')+'</option>').val(name+(as_default?' [default]':''));
      // Success message
      $('<p class="text-success" id="mqc-save-success">Settings saved.</p>').hide().insertBefore($('#mqc_saveconfig_form')).slideDown(function(){
        setTimeout(function(){
          $('#mqc-save-success').slideUp(function(){ $(this).remove(); });
        }, 5000);
      });
    }
  } catch(e){ console.log('Error updating localstorage: '+e); }
}

// Clear current default configuration
function mqc_clear_default_config() {
  try {
    var config = localStorage.getItem("mqc_config");
    if (!config) {
      return;
    } else {
      config = JSON.parse(config);
    }
    for (var c in config) {
      if (config.hasOwnProperty(c)) {
        config[c]['default'] = false;
      }
    }
    localStorage.setItem("mqc_config", JSON.stringify(config));
    $('<p class="text-danger" id="mqc-cleared-success">Unset default.</p>').hide().insertBefore($('#mqc_loadconfig_form .actions')).slideDown(function () {
      setTimeout(function () {
        $('#mqc-cleared-success').slideUp(function () { $(this).remove(); });
      }, 5000);
      var name = $('#mqc_loadconfig_form select option:contains("default")').text();
      $('#mqc_loadconfig_form select option:contains("default")').remove();
      name = name.replace(' [default]', '');
      $('#mqc_loadconfig_form select').append('<option>'+name+'</option>').val(name);
    });
  } catch (e) {
    console.log('Could not access localStorage');
  }
}

//////////////////////////////////////////////////////
// LOAD TOOLBOX SAVE NAMES
//////////////////////////////////////////////////////
function populate_mqc_saveselect(){
  var default_config = '';
  try {
    var local_config = localStorage.getItem("mqc_config");
    if(local_config !== null && local_config !== undefined){
      local_config = JSON.parse(local_config);
      default_name = false;
      for (var name in local_config){
        if (local_config[name]['default']) {
          console.log('Loaded default config!');
          load_mqc_config(name);
          default_config = name;
          name = name+' [default]';
          default_name = name;
        }
        $('#mqc_loadconfig_form select').append('<option>'+name+'</option>').val(name);
      }
      // Set the selected select option
      if(default_name !== false){
        $('#mqc_loadconfig_form select option:contains("'+default_name+'")').prop('selected',true);
      } else {
        $('#mqc_loadconfig_form select option:first').prop('selected',true);
      }
    }
  } catch(e){
    console.log('Could not load local config: '+e);
    $('#mqc_saveconfig').html('<h4>Error accessing localStorage</h4>'+
      '<p>This feature uses a web browser feature called "localStorage". '+
      "We're not able to access this at the moment, which probably means that "+
      'you have the <em>"Block third-party cookies and site data"</em> setting ticked (Chrome) '+
      'or equivalent in other browsers.</p><p>Please '+
      '<a href="https://www.google.se/search?q=Block+third-party+cookies+and+site+data" target="_blank">change this browser setting</a>'+
      ' to save MultiQC report configs.</p>');
  }
}

//////////////////////////////////////////////////////
// LOAD TOOLBOX SETTINGS
//////////////////////////////////////////////////////
function load_mqc_config(name){
  if(name === undefined){ return false; }
  var config = {};
  try {
    try {
      var local_config = localStorage.getItem("mqc_config");
    } catch(e){ console.log('Could not access localStorage'); }
    if(local_config !== null && local_config !== undefined){
      local_config = JSON.parse(local_config);
      for (var attr in local_config[name]) {
        config[attr] = local_config[name][attr];
      }
    }
  } catch(e){ console.log('Could not load local config: '+e); }

  // Apply config - rename samples
  if(notEmptyObj(config['rename_regex'])){
    if(config['rename_regex'] == true){
      $('#mqc_renamesamples .mqc_regex_mode .re_mode').removeClass('off').addClass('on').text('on');
      window.mqc_rename_regex_mode = true;
    }
  }
  if(notEmptyObj(config['rename_from_texts']) && notEmptyObj(config['rename_to_texts'])){
    window.mqc_rename_f_texts = [];
    window.mqc_rename_t_texts = [];
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
    $(document).trigger('mqc_renamesamples', [window.mqc_rename_f_texts, window.mqc_rename_t_texts, config['rename_regex']]);
  }

  // Apply config - highlights
  if(notEmptyObj(config['highlight_regex'])){
    if(config['highlight_regex'] == true){
      $('#mqc_cols .mqc_regex_mode .re_mode').removeClass('off').addClass('on').text('on');
      window.mqc_highlight_regex_mode = true;
    }
  }
  if(notEmptyObj(config['highlights_f_texts']) && notEmptyObj(config['highlights_f_cols'])){
    window.mqc_highlight_f_texts = [];
    window.mqc_highlight_f_cols = [];
    $.each(config['highlights_f_texts'], function(idx, f_text){
      var f_col = config['highlights_f_cols'][idx];
      $('#'+hashCode(f_text+f_col)).remove();
      $('#mqc_col_filters').append('<li style="color:'+f_col+';" id="'+hashCode(f_text+f_col)+'"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
      window.mqc_highlight_f_texts.push(f_text);
      window.mqc_highlight_f_cols.push(f_col);
      mqc_colours_idx += 1;
    });
    $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
    $(document).trigger('mqc_highlights', [window.mqc_highlight_f_texts, window.mqc_highlight_f_cols, config['highlight_regex']]);
  }

  // Apply config - hide samples
  if(notEmptyObj(config['hidesamples_regex'])){
    if(config['hidesamples_regex'] == true){
      $('#mqc_hidesamples .mqc_regex_mode .re_mode').removeClass('off').addClass('on').text('on');
      window.mqc_hide_regex_mode = true;
    }
  }
  if(notEmptyObj(config['hidesamples_mode'])){
    if(config['hidesamples_mode'] == 'show'){
      $('.mqc_hidesamples_showhide').prop('checked', false);
      $('.mqc_hidesamples_showhide[val=show]').prop('checked', true);
      window.mqc_hide_mode = 'show';
    }
  }
  if(notEmptyObj(config['hidesamples_f_texts'])){
    window.mqc_hide_f_texts = [];
    $.each(config['hidesamples_f_texts'], function(idx, f_text){
      if(f_text.length == 0){ return true; }
      $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
      window.mqc_hide_f_texts.push(f_text);
    });
    $(document).trigger('mqc_hidesamples', [window.mqc_hide_f_texts, config['hidesamples_regex']]);
  }

  // Trigger loaded event to initialise plots
  $(document).trigger('mqc_config_loaded');

}
