////////////////////////////////////////////////
// Static MatPlotLib Plots Javascript Code
////////////////////////////////////////////////

// On page load
$(function () {
  // Switch between counts and percentages in a bar plot
  $(".mpl_switch_group.percent-switch").click(function (e) {
    e.preventDefault();
    $(this).toggleClass("active");
    let plotgroup = $(this).closest(".mqc_mplplot_plotgroup");
    let current = "#" + plotgroup.find(".mqc_mplplot:visible").attr("id");

    let target;
    if (current.endsWith("_pc_log")) {
      target = current.slice(0, -"_pc_log".length) + "_log";
    } else if (current.endsWith("_pc")) {
      target = current.slice(0, -"_pc".length);
    } else if (current.endsWith("_log")) {
      target = current.slice(0, -"_log".length) + "_pc_log";
    } else {
      target = current + "_pc";
    }
    plotgroup.find(".mqc_mplplot").hide();
    $(target).show();
  });

  // Switch log scale in a bar plot
  $(".mpl_switch_group.log10-switch").click(function (e) {
    e.preventDefault();
    $(this).toggleClass("active");
    let plotgroup = $(this).closest(".mqc_mplplot_plotgroup");
    let current = "#" + plotgroup.find(".mqc_mplplot:visible").attr("id");

    let target;
    if (current.endsWith("_log")) target = current.slice(0, -"_log".length);
    else target = current + "_log";

    plotgroup.find(".mqc_mplplot").hide();
    $(target).show();
  });

  // Switch datasets in a bar plot
  $(".mpl_switch_group.dataset-switch-group button").click(function (e) {
    e.preventDefault();
    if ($(this).hasClass("active")) return;
    $(this).siblings("button.active").removeClass("active");
    $(this).addClass("active");
    let ds_uid = $(this).data("datasetUid");
    let target_id = "#" + ds_uid;
    let plotgroup = $(target_id).closest(".mqc_mplplot_plotgroup");
    if (plotgroup.children(".flat-percent-switch").hasClass("active")) {
      target_id += "_pc";
    } else if (plotgroup.children(".flat-log10-switch").hasClass("active")) {
      target_id += "_log";
    }
    plotgroup.find(".mqc_mplplot").hide();
    $(target_id).show();
  });
});
