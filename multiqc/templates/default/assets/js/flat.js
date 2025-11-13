////////////////////////////////////////////////
// Flat plots controls
////////////////////////////////////////////////

// On page load
$(function () {
  function switchPlot(button) {
    let plotGroup = $(button).closest(".mqc_mplplot_plotgroup");
    let plotAnchor = plotGroup.data("plot-anchor");
    let activeDs = plotGroup.find(".dataset-switch-group button.active");
    if (activeDs !== undefined && activeDs.length > 0) plotAnchor = activeDs.data("datasetUid");

    let p = plotGroup.find(".percent-switch");
    let l = plotGroup.find(".log10-switch");
    if (p.length > 0 || l.length > 0) {
      let pActive = p.length > 0 && p.hasClass("active");
      let lActive = l.length > 0 && l.hasClass("active");
      if (pActive) plotAnchor += "-pct";
      if (lActive) plotAnchor += "-log";
      if (!pActive && !lActive) plotAnchor += "-cnt";
    }

    plotGroup.find(".mqc_mplplot").hide();
    $("#" + plotAnchor).show();
  }

  // Switch between counts and percentages in a bar plot
  $(".mpl_switch_group.percent-switch,.mpl_switch_group.log10-switch").click(function (e) {
    e.preventDefault();
    $(this).toggleClass("active");
    switchPlot(this);
  });

  // Switch datasets in a bar plot
  $(".mpl_switch_group.dataset-switch-group button").click(function (e) {
    e.preventDefault();
    if ($(this).hasClass("active")) return;
    $(this).siblings("button.active").removeClass("active");
    $(this).addClass("active");
    switchPlot(this);
  });
});
