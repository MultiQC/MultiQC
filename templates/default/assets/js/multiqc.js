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
        if(maxval === undefined){ maxval = Math.max.apply(Math, data); }
        var minval = $(this).data('chroma-min');
        if(minval === undefined){ minval = Math.min.apply(Math, data); }

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

// Basic Line Graph
function plot_xy_line_graph(div, data, title, ylab, xlab, ymax, ymin, xmax, xmin, tt_label){
  if(tt_label === undefined){ tt_label = '{point.x}'; }
  $(div).highcharts({
    chart: {
      type: 'line',
      zoomType: 'x'
    },
    title: {
      text: title,
      x: -20 //center
    },
    xAxis: {
      title: {
        text: xlab
      },
      max: xmax,
      min: xmin
    },
    yAxis: {
      title: {
        text: ylab
      },
      max: ymax,
      min: ymin,
      plotLines: [{
        value: 0,
        width: 1,
        color: '#808080'
      }]
    },
    plotOptions: {
      series: {
        marker: {
          enabled: false
        }
      }
    },
    legend: {
      layout: 'vertical',
      align: 'right',
      verticalAlign: 'middle',
      borderWidth: 0
    },
    credits: {
			enabled: false
		},
    tooltip: {
      headerFormat: '<b>'+tt_label+'</b><table>',
			pointFormat: '<tr><td><span style="color:{series.color};">{series.name}:</span></td><td>{point.y:.2f}</td></tr>',
			footerFormat: '</table>',
			useHTML: true
    },
    series: data
  });
}



//////////////////
// Generic helper functions
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
