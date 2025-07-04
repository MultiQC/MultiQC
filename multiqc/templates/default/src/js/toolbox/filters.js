////////////////////////////////////////////////
// MultiQC Report Common Filter Handling
//
// This file contains common code for handling filter rows across different
// toolbox sections (highlights, rename samples, hide samples). The filter
// functionality includes adding/removing filter rows, handling input events,
// and managing regex mode toggles.
////////////////////////////////////////////////

// Make function available globally
window.initFilters = function () {
  // Filter text is changed
  $(".mqc_filters").on("blur", "li input", function () {
    var target = $(this).parent().parent().attr("id");
    if (target == "mqc_col_filters") {
      $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if (target == "mqc_renamesamples_filters") {
      $("#mqc_rename_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if (target == "mqc_hidesamples_filters") {
      $("#mqc_hide_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
  });

  // 'Enter' key pressed whilst editing a filter
  $(".mqc_filters").on("keyup", "li input", function (e) {
    if (e.keyCode == 13) {
      // Pressed enter
      $(this).blur();
      $(this).parent().next("li").find("input").focus().select();
    }
  });

  // Remove filter button
  $(".mqc_filters").on("click", "li button", function () {
    var target = $(this).parent().parent().attr("id");
    $(this).parent().remove();
    if (target == "mqc_col_filters") {
      $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if (target == "mqc_hidesamples_filters") {
      $("#mqc_hide_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if (target == "mqc_renamesamples_filters") {
      $("#mqc_rename_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
  });

  // Clear all filters button
  $(".mqc_toolbox_clear").click(function () {
    var target = $(this).closest(".mqc_filter_section").find(".mqc_filters").attr("id");
    $("#" + target).empty();
    if (target == "mqc_col_filters") {
      $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if (target == "mqc_hidesamples_filters") {
      $("#mqc_hide_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if (target == "mqc_renamesamples_filters") {
      $("#mqc_rename_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
  });

  // Regex mode text
  $(".mqc_regex_mode").click(function () {
    var rswitch = $(this).find(".re_mode");
    if (rswitch.text() == "off") {
      rswitch.removeClass("off").addClass("on").text("on");
    } else {
      rswitch.removeClass("on").addClass("off").text("off");
    }
    if ($(this).data("target") == "mqc_cols") {
      $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if ($(this).data("target") == "mqc_renamesamples") {
      $("#mqc_rename_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
    if ($(this).data("target") == "mqc_hidesamples") {
      $("#mqc_hide_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    }
  });
};
