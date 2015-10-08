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
  
  // Freeze the top header when scrolling
  var gsTab = $('#general_stats_table');
  var gsTabDiv = $('#general_stats_table_container .table-responsive');
  var gsHeight = gsTab.height();
  
  $(window).scroll(function(){
    var wTop = $(window).scrollTop();
    var gsOffset = gsTab.offset().top;
    var gsTop = gsOffset - wTop;
    var gsVisible = gsTop + gsHeight;
    var gsHeadHeight = 30;
    if(gsTop < 0 && gsVisible > 0 ){
      ctw = $("#gsCloneWrapper");
      if(ctw.length == 0){
        // Make a copy of the general stats table
        ct = gsTab.clone();
        ct.attr('id', 'gsClone').width(gsTab.width());
        // Hide everything except the header. Scroll it sideways if needed.
        ct.css({visibility:'hidden', 'margin-left': -gsTabDiv.scrollLeft()});
        ct.find('thead').css({visibility:'visible'});
        // Wrap it and add to the container with position: fixed
        ctw = $('<div id="gsCloneWrapper" />').append(ct);
        ctw.css({'position':'fixed', 'top':0, 'height': gsHeadHeight});
        $("#general_stats_table_container").append(ctw);
      }
      // Nicely scroll out of the way instead of dissapearing
      if(gsVisible < gsHeadHeight * 2){
        $("#gsCloneWrapper").css('top', gsVisible - (gsHeadHeight * 2) );
      } else {
        $("#gsCloneWrapper").css('top', 0);
      }
    } else {
      // Not needed - remove it (avoids printing errors etc)
      $("#gsCloneWrapper").remove();
    }
  });
  // Scroll left and right in the responsive container
  gsTabDiv.scroll(function(){
    $("#gsClone").css('margin-left', -$(this).scrollLeft());
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

        // Get the max and min values from data attributes
        var maxval = $(this).data('chroma-max');
        var minval = $(this).data('chroma-min');
        if(isNaN(minval) || isNaN(maxval)){
          console.log('Could not find max or min value for '+$(this).text()+': ['+[minval, maxval]+']')
          return true; // Skip to next loop
        }

        // Go through table cells again, adding colour
        var i = 0;
        var scale = chroma.scale(colscheme).domain([minval, maxval]);
        if(colscheme_rev){
          scale = chroma.scale(colscheme).domain([maxval, minval]);
        }
        table.find('tr td:nth-of-type('+(idx)+')').each(function(){
          var val = parseFloat($(this).text());
          var col = scale(val).css();
          $(this).find('.wrapper .bar').css('background-color', col);
        });
        
      }
    });
  });
    
});
