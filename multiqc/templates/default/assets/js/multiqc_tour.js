////////////////////////////////////////////////
// MultiQC Tour, using bootstrap tour
////////////////////////////////////////////////

// Execute when page load has finished
$(function () {
 
  // Introduction tour
  $('#mqc-launch-into-tour').after(' <small><em>('+tour_steps.length+' steps - takes around 1 minute)</em></small>');
  $.each(tour_steps, function(i, step){
    step['title'] += '<span class="pull-right">'+(i+1)+'/'+tour_steps.length+'</span>';
    var percent = parseInt(((i+1) / tour_steps.length) * 100);
    step['content'] = '<div class="pbar_wrapper"><hr class="pbar" style="width:'+percent+'%;"></div>' + step['content'];
  });
  orig_z_index = $('.mqc-toolbox').css('z-index'); 
  var intro_tour = new Tour({
    backdrop: true,
    storage: false,
    onStart: function(tour){
      $('.mqc-toolbox').css('z-index', 0);
      mqc_toolbox_openclose('#mqc_cols', false);
    },
    onEnd: function (tour) {
      $('.mqc-toolbox').css('z-index', 1200);
      mqc_toolbox_openclose('#mqc_cols', false);
    },
    steps: tour_steps
  });
  intro_tour.init();
  $('#mqc-launch-into-tour').click(function(ev){
    ev.preventDefault();
    try{ intro_tour.restart(); }
    catch(e){ console.log('Tour broke - '+e); }
  });
  
});

var orig_z_index = 1200;
var tour_steps = [
  {
    orphan: true,
    title: "Welcome to MultiQC!",
    content: 'MultiQC generates reports based on analysis across many samples.<br>Click next to explore the tools availble in this report.',
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
  {
    element: ".hc-plot-handle:first",
    placement: 'top',
    title: "Resize Plots",
    content: "Drag the grey bar below plots to change their height.",
    backdropPadding: 10,
  },
  {
    element: ".highcharts-button:first path",
    title: "Export Plots",
    placement: 'left',
    backdropPadding: {
      'top': 5,
      'left': 5,
      'bottom': 35,
      'right': 35,
    },
    content: "Plots can be exported in a range of formats (including <code>svg</code> and <code>pdf</code>, suitable for publications).",
    onHide: function (tour) { $('.mqc-toolbox').css('z-index', orig_z_index); }
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
    content: "You can save your configuration for this report, or as a default for all reports. You can also share your set up with others via a downloaded file.",
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
  }];