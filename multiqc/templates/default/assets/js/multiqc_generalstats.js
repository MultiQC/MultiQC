////////////////////////////////////////////////
// MultiQC General Statistics Table code
////////////////////////////////////////////////

var brewer_scales = ['YlOrRd', 'YlOrBr', 'YlGnBu', 'YlGn', 'Reds', 'RdPu',
  'Purples', 'PuRd', 'PuBuGn', 'PuBu', 'OrRd', 'Oranges', 'Greys', 'Greens',
  'GnBu', 'BuPu', 'BuGn', 'Blues', 'Set3', 'Set2', 'Set1', 'Pastel2', 'Pastel1',
  'Paired', 'Dark2', 'Accent', 'Spectral', 'RdYlGn', 'RdYlBu', 'RdGy', 'RdBu',
  'PuOr', 'PRGn', 'PiYG', 'BrBG'];

// Execute when page load has finished loading
$(function () {
    
  // Enable the bootstrap tooltip hovers
  // http://getbootstrap.com/javascript/#tooltips
  $('[data-toggle="tooltip"]').tooltip();
  
  // Enable tablesorter on the general statistics table
  $("#general_stats_table").tablesorter({sortInitialOrder: 'desc'}); 
  
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
            if(percentage > 100){ percentage = 100; }
            if(percentage < 0){ percentage = 0; }
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
    
});
