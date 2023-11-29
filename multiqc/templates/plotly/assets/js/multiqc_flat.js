////////////////////////////////////////////////
// Static MatPlotLib Plots Javascript Code
////////////////////////////////////////////////

// On page load
$(function () {
  // Switch between counts and percentages in a bar plot
  $(".flat_plot_toggle_percent_log button").click(function (e) {
    e.preventDefault();
    if (!$(this).hasClass("active")) {
      $(this).siblings("button.active").removeClass("active");
      $(this).addClass("active");
      let wrapper = $(this).closest(".mqc_mplplot_plotgroup");
      let current = "#" + wrapper.find(".mqc_mplplot:visible").attr("id");
      let target;
      if (current.substr(current.length - 3) === "_pc") {
        target = current.substr(0, current.length - 3);
      } else if (current.substr(current.length - 4) === "_log") {
        target = current.substr(0, current.length - 4);
      } else {
        target = current + "_pc";
      }
      wrapper.find(".mqc_mplplot").hide();
      $(target).show();
    }
  });

  // Switch datasets in a bar plot
  $(".flat_plot_switch_datasets button").click(function (e) {
    e.preventDefault();
    if (!$(this).hasClass("active")) {
      $(this).siblings("button.active").removeClass("active");
      $(this).addClass("active");
      let target = $(this).data("target");
      let wrapper = $(target).closest(".mqc_mplplot_plotgroup");
      if (wrapper.find(".flat_plot_toggle_percent_log .percents").hasClass("active")) {
        target += "_pc";
      } else if (wrapper.find(".flat_plot_toggle_percent_log .log").hasClass("active")) {
        target += "_log";
      }
      wrapper.find(".mqc_mplplot").hide();
      $(target).show();
    }
  });
});
