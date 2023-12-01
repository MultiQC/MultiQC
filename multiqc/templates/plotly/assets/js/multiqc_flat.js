////////////////////////////////////////////////
// Static MatPlotLib Plots Javascript Code
////////////////////////////////////////////////

// On page load
$(function () {
  // Switch between counts and percentages in a bar plot
  $(".mpl_switch_group.percent-switch,.mpl_switch_group.log10-switch").click(function (e) {
    e.preventDefault();
    $(this).toggleClass("active");

    let plotgroup = $(this).closest(".mqc_mplplot_plotgroup");
    let percent_switch = plotgroup.children(".percent-switch");
    let log10_switch = plotgroup.children(".log10-switch");

    let pid = plotgroup[0].id;
    let activeDatasetId = plotgroup.find(".dataset-switch-group button.active").data("datasetUid");
    if (activeDatasetId !== undefined) pid = activeDatasetId;
    let target = "#flat-" + pid;

    if (percent_switch.hasClass("active")) target += "_pc";
    if (log10_switch.hasClass("active")) target += "_log";

    plotgroup.find(".mqc_mplplot").hide();
    $(target).show();
  });

  // Switch datasets in a bar plot
  $(".mpl_switch_group.dataset-switch-group button").click(function (e) {
    e.preventDefault();
    if ($(this).hasClass("active")) return;
    $(this).siblings("button.active").removeClass("active");
    $(this).addClass("active");

    let plotgroup = $(this).closest(".mqc_mplplot_plotgroup");
    let target = "#flat-" + $(this).data("datasetUid");
    if (plotgroup.children(".percent-switch").hasClass("active")) target += "_pc";
    if (plotgroup.children(".log10-switch").hasClass("active")) target += "_log";

    plotgroup.find(".mqc_mplplot").hide();
    $(target).show();
  });
});
