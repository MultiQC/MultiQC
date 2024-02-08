////////////////////////////////////////////////
// Flat plots controls
////////////////////////////////////////////////

// On page load
$(function () {
  function switchPlot(button) {
    let plotgroup = $(button).closest(".mqc_mplplot_plotgroup");
    let pid = plotgroup.data("pid");
    let activeDs = plotgroup.find(".dataset-switch-group button.active");
    if (activeDs !== undefined && activeDs.length > 0) pid = activeDs.data("datasetUid");

    let target = "#" + pid;

    let p = plotgroup.find(".percent-switch");
    let l = plotgroup.find(".log10-switch");
    if (p.length > 0 || l.length > 0) {
      let pActive = p.length > 0 && p.hasClass("active");
      let lActive = l.length > 0 && l.hasClass("active");
      if (pActive) target += "-pct";
      if (lActive) target += "-log";
      if (!pActive && !lActive) target += "-cnt";
    }

    plotgroup.find(".mqc_mplplot").hide();
    $(target).show();
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
