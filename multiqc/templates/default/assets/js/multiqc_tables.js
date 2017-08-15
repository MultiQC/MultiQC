////////////////////////////////////////////////
// MultiQC Table code
////////////////////////////////////////////////

var brewer_scales = ['YlOrRd', 'YlOrBr', 'YlGnBu', 'YlGn', 'Reds', 'RdPu',
  'Purples', 'PuRd', 'PuBuGn', 'PuBu', 'OrRd', 'Oranges', 'Greys', 'Greens',
  'GnBu', 'BuPu', 'BuGn', 'Blues', 'Set3', 'Set2', 'Set1', 'Pastel2', 'Pastel1',
  'Paired', 'Dark2', 'Accent', 'Spectral', 'RdYlGn', 'RdYlBu', 'RdGy', 'RdBu',
  'PuOr', 'PRGn', 'PiYG', 'BrBG'];

// Execute when page load has finished loading
$(function () {

  if($('.mqc_table').length > 0){

    // Enable tablesorter on MultiQC tables
    var strip_non_numeric = function(node){
      return node.innerText.replace(/[^\d.-]/g, '');
    }
    $('.mqc_table').tablesorter({sortInitialOrder: 'desc', textExtraction: strip_non_numeric});

    // Update tablesorter if samples renamed
    $(document).on('mqc_renamesamples', function(e, f_texts, t_texts, regex_mode){
      $('.mqc_table').trigger('update');
    });

    // Copy table contents to clipboard
    var clipboard = new Clipboard('.mqc_table_copy_btn');
    clipboard.on('success', function(e) { e.clearSelection(); });
    $('.mqc_table_copy_btn').click(function(){
      var btn = $(this);
      btn.addClass('active').html('<span class="glyphicon glyphicon-copy"></span> Copied!');
      setTimeout(function(){
        btn.removeClass('active').html('<span class="glyphicon glyphicon-copy"></span> Copy table');
      }, 2000);
    });

    // Make table headers fixed when table body scrolls (use CSS transforms)
    // http://stackoverflow.com/a/25902860/713980
    $('.mqc-table-responsive').scroll(function() {
      $(this).find('thead').css('transform', "translate(0,"+$(this).scrollTop()+"px)");
    });

    // Table header-specific bootstrap tooltips
    $('.mqc_table_tooltip').tooltip({ container: 'body' });

    // Expand tables to full height
    $('.mqc-table-expand').click(function(){
      if($(this).find('span').hasClass('glyphicon-chevron-down')){
        $(this).parent().find('.mqc-table-responsive').css('max-height', 'none');
        $(this).find('span').removeClass('glyphicon-chevron-down').addClass('glyphicon-chevron-up');
      } else {
        $(this).parent().find('.mqc-table-responsive').css('max-height', '400px');
        $(this).find('span').removeClass('glyphicon-chevron-down').addClass('glyphicon-chevron-down');
      }
    });

    /////// COLUMN CONFIG
    // show + hide columns
    $('.mqc_table_col_visible').change(function(){
      var target = $(this).data('target');
      mqc_table_col_updateVisible(target);
    });
    // Bulk set visible / hidden
    $('.mqc_configModal_bulkVisible').click(function(e){
      e.preventDefault();
      var target = $(this).data('target');
      var visible = $(this).data('action') == 'showAll';
      $(target+'_configModal_table tbody .mqc_table_col_visible').prop('checked', visible);
      mqc_table_col_updateVisible(target);
    });
    function mqc_table_col_updateVisible(target){
      $(target+'_configModal_table .mqc_table_col_visible').each(function(){
        var cclass = $(this).val();
        if($(this).is(":checked")) {
          $(target+' .'+cclass).removeClass('hidden');
          $(target+'_configModal_table .'+cclass).removeClass('text-muted');
        } else {
          $(target+' .'+cclass).addClass('hidden');
          $(target+'_configModal_table .'+cclass).addClass('text-muted');
        }
      });
      // Hide empty rows
      $(target+' tbody tr').show();
      $(target+' tbody tr').each(function(){
        var hasVal = false;
        $(this).find('td:visible').each(function(){
          if(!$(this).hasClass('sorthandle') && $(this).text() !== ''){
            hasVal = true;
          }
        });
        if(!hasVal){
          $(this).hide();
        }
      });
      // Update counts
      $(target+'_numrows').text( $(target+' tbody tr:visible').length );
      $(target+'_numcols').text( $(target+' thead th:visible').length - 1 );
    }

    // Make rows in MultiQC tables sortable
    $('.mqc_table.mqc_sortable tbody').sortable({
      handle: '.sorthandle',
      helper: function fixWidthHelper(e, ui) {
        ui.children().each(function() { $(this).width($(this).width()); });
        return ui;
      }
    });

    // Change order of columns
    $('.mqc_configModal_table').on('sortstop', function(e, ui){
      change_mqc_table_col_order( $(this) );
    });
    $('.mqc_configModal_table').bind('sortEnd',function() { 
      change_mqc_table_col_order( $(this) );
    });

    // TOOLBOX LISTENERS

    // highlight samples
    $(document).on('mqc_highlights', function(e, f_texts, f_cols, regex_mode){
      $('.mqc_table_sortHighlight').hide();
      $('.mqc_table tbody th').removeClass('highlighted').removeData('highlight');
      $('.mqc_table tbody th').each(function(i){
        var th = $(this);
        var thtext = $(this).text();
        var thiscol = '#333';
        $.each(f_texts, function(idx, f_text){
          if((regex_mode && thtext.match(f_text)) || (!regex_mode && thtext.indexOf(f_text) > -1)){
            thiscol = f_cols[idx];
            th.addClass('highlighted').data('highlight', idx);
            $('.mqc_table_sortHighlight').show();
          }
        });
        $(this).css('color', thiscol);
      });
    });

    // Sort MultiQC tables by highlight
    $('.mqc_table_sortHighlight').click(function(e){
      e.preventDefault();
      var target = $(this).data('target');
      // collect highlighted rows
      var hrows = $(target+' tbody th.highlighted').parent().detach();
      hrows = hrows.sort(function (a, b) {
        return $(a).find('th').data('highlight') - $(b).find('th').data('highlight');
      });
      if($(this).data('direction') == 'desc'){
        hrows = hrows.get().reverse();
        $(target+' tbody').prepend(hrows);
        $(this).data('direction', 'asc');
      } else {
        $(target+' tbody').append(hrows);
        $(this).data('direction', 'desc');
      }
    });

    // Rename samples
    $(document).on('mqc_renamesamples', function(e, f_texts, t_texts, regex_mode){
      $(".mqc_table tbody th").each(function(){
        var s_name = $(this).data('original-sn');
        $.each(f_texts, function(idx, f_text){
          if(regex_mode){
            var re = new RegExp(f_text,"g");
            s_name = s_name.replace(re, t_texts[idx]);
          } else {
            s_name = s_name.replace(f_text, t_texts[idx]);
          }
        });
        $(this).text(s_name);
      });
    });

    // Hide samples
    $(document).on('mqc_hidesamples', function(e, f_texts, regex_mode){

      // Hide rows in MultiQC tables
      $(".mqc_table tbody th").each(function(){
        var match = false;
        var hfilter = $(this).text();
        $.each(f_texts, function(idx, f_text){
          if((regex_mode && hfilter.match(f_text)) || (!regex_mode && hfilter.indexOf(f_text) > -1)){
            match = true;
          }
        });
        if(window.mqc_hide_mode == 'show'){
          match = !match;
        }
        if(match){
          $(this).parent().hide().addClass('hidden');
        } else {
          $(this).parent().show().removeClass('hidden');
        }
      });
      $('.mqc_table_numrows').each(function(){
        var tid = $(this).attr('id').replace('_numrows','');
        $(this).text( $('#'+tid+' tbody tr:visible').length );
      });

      // Hide empty columns
      $('.mqc_table').each(function(){
        var table = $(this);
        var gsthidx = 0;
        table.find("thead th, tbody tr td").show();
        table.find("thead th").each(function(){
          if(gsthidx == 0){ gsthidx += 1; return true; }
          var count = 0;
          var empties = 0;
          table.find("tbody tr td:nth-child("+(gsthidx+2)+")").filter(":visible").each(function(){
            count += 1;
            if($(this).text() == ''){ empties += 1; }
          });
          if(count > 0 && count == empties){
            $(this).hide();
            table.find("tbody tr td:nth-child("+(gsthidx+2)+")").hide();
          }
          gsthidx += 1;
        });
      });
      $('.mqc_table_numcols').each(function(){
        var tid = $(this).attr('id').replace('_numcols','');
        $(this).text( $('#'+tid+' thead th:visible').length - 1 );
      });
    });

  } // End of check for table

  // Table Scatter Modal
  $('#tableScatterForm').submit(function(e){
    e.preventDefault();
  });
  $('.mqc_table_makeScatter').click(function(e){
    // Reset dropdowns
    if($('#tableScatter_tid').val() != $(this).data('table')){
      $('#tableScatter_col1, #tableScatter_col2').html('<option value="">Select Column</option>');
      // Add columns to dropdowns
      $($(this).data('table')+' thead tr th').each(function(e){
        var c_id = $(this).attr('id');
        if(c_id != undefined){
          var c_name = $(this).attr('data-namespace') + ': ' + $(this).text();
          $('#tableScatter_col1, #tableScatter_col2').append('<option value="'+c_id+'">'+c_name+'</select>');
        }
      });
      $('#tableScatter_tid').val($(this).data('table'));
      $('#tableScatterPlot').html('<small>Please select two table columns.</small>').addClass('not_rendered');
    }
  });
  $('#tableScatterForm select').change(function(e){
    var tid = $('#tableScatter_tid').val();
    var col1 = $('#tableScatter_col1').val().replace('header_', '');
    var col2 = $('#tableScatter_col2').val().replace('header_', '');
    var col1_name = $('#tableScatter_col1 option:selected').text();
    var col2_name = $('#tableScatter_col2 option:selected').text();
    var col1_max = parseFloat($(tid+' thead th#header_'+col1).data('dmax'));
    var col1_min = parseFloat($(tid+' thead th#header_'+col1).data('dmin'));
    var col2_max = parseFloat($(tid+' thead th#header_'+col2).data('dmax'));
    var col2_min = parseFloat($(tid+' thead th#header_'+col2).data('dmin'));
    if(isNaN(col1_max)){ col1_max = undefined; }
    if(isNaN(col1_min)){ col1_min = undefined; }
    if(isNaN(col2_max)){ col2_max = undefined; }
    if(isNaN(col2_min)){ col2_min = undefined; }
    if(col1 != '' && col2 != ''){
      $('#tableScatterPlot').html('<small>loading..</small>');
      if ($(tid).attr('data-title')) {
        plot_title = $(tid).attr('data-title');
      } else {
        plot_title = tid.replace(/^#/, '').replace(/_/g, ' ');
      }
      // Get the data values
      mqc_plots['tableScatterPlot'] = {
        'plot_type': 'scatter',
        'config': {
          'id': 'tableScatter_'+tid,
          'title': plot_title,
          'xlab': col1_name,
          'ylab': col2_name,
          'xmin': col1_min,
          'xmax': col1_max,
          'ymin': col2_min,
          'ymax': col2_max,
        },
        'datasets': [[]]
      };
      $(tid+' tbody tr').each(function(e){
        var s_name = $(this).children('th.rowheader').text();
        var val_1 = $(this).children('td.'+col1).text().replace(/[^\d\.]/g,'');
        var val_2 = $(this).children('td.'+col2).text().replace(/[^\d\.]/g,'');
        if(!isNaN(parseFloat(val_1)) && isFinite(val_1) && !isNaN(parseFloat(val_2)) && isFinite(val_2)){
          mqc_plots['tableScatterPlot']['datasets'][0].push({
            'name': s_name,
            'x': parseFloat(val_1),
            'y': parseFloat(val_2)
          });
        }
      });
      if(Object.keys(mqc_plots['tableScatterPlot']['datasets'][0]).length > 0){
        if(plot_scatter_plot('tableScatterPlot') == false){
          $('#tableScatterPlot').html('<small>Error: Something went wrong when plotting the scatter plot.</small>');
          $('#tableScatterPlot').addClass('not_rendered');
        } else {
          $('#tableScatterPlot').removeClass('not_rendered');
        }
      } else {
        $('#tableScatterPlot').html('<small>Error: No data pairs found for these columns.</small>');
        $('#tableScatterPlot').addClass('not_rendered');
      }
    } else {
      $('#tableScatterPlot').html('<small>Please select two table columns.</small>');
      $('#tableScatterPlot').addClass('not_rendered');
    }
  });

});

// Reorder columns in MultiQC tables.
// Note: Don't have to worry about floating headers, as 'Configure Columns'
// button is only visible when this is hidden. Ace!
function change_mqc_table_col_order(table){

  // Find the targets of this sorting
  var tid = table.attr('id');
  var target = tid.replace('_configModal_table','');

  // Collect the desired order of columns
  var classes = [];
  $('#'+tid+' tbody tr').each(function(){
    classes.push($(this).attr('class'));
  });
  // Go through each row
  $('#'+target+' tr').each(function(){
    var cols = {};
    var row = $(this);
    // Detach any cell that matches a known class from above
    row.find('td, th').each(function(){
      var cell = $(this);
      $.each(classes, function(idx, c){
        if(cell.hasClass(c)){
          cols[c] = cell.detach();
        }
      });
    });
    // Insert detached cells back in the order given in the sorted table
    for (var idx in classes){
      var c = classes[idx];
      if(cols[c] !== undefined){
        row.append(cols[c]);
      }
    }
  });
}
