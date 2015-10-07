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
          
        }
      });
    });
    
});
