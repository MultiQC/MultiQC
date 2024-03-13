// Functions

//////////////////////////////////////////////////
// Button Status Switcher

////////////////////////
// Button Disabler
function btn_disable(mode, regex, hide_list, user_hide_list) {
  // List of example buttons / items, this will be pretty much constant throughout life of HTStream
  var exempt_list = [
    "Base: A",
    "Base: C",
    "Base: G",
    "Base: T",
    "Base: N",
    "PE Reads",
    "SE Reads",
    "PE Bps",
    "SE Bps",
    "Read 1",
    "Read 2",
    "Single End",
  ];

  // Find htstream report section and appropriate buttons
  var parent_div = $("#mqc-module-section-htstream");
  var unfiltered_btn_divs = parent_div.find("[data-action='set_data'], *[id*=htstream_]");
  var btn_divs = unfiltered_btn_divs.filter(":button").filter(":not(*[onClick*=htstream_plot_switch])");

  // Disable buttons that are not in exempt list or specifically listed to hide by MultiQC
  $.each(btn_divs.get().reverse(), function (x, value) {
    var btn_text = $(value).text();

    // Mode, show or hide
    if (mode == "show") {
      // Is regex on
      if (regex) {
        var show = false;

        // Iterate throughlist
        for (i = 0; i < user_hide_list.length; i++) {
          // Show specified buttons (regex)
          if (btn_text.match(user_hide_list[i])) {
            $(value).prop("disabled", false);
            show = true;
            break;
          }
        }

        // Hide buttons not included in exempt and specified lists
        if (show == false && !exempt_list.includes(btn_text)) {
          $(value).prop("disabled", true);
        }
      } else {
        // Show specified buttons, disable nonspecified
        if (hide_list.indexOf(btn_text) == -1 && exempt_list.indexOf(btn_text) == -1) {
          $(value).prop("disabled", true);
        } else {
          $(value).prop("disabled", false);
        }
      }
    } else {
      // Regex?
      if (regex) {
        show = true;

        // Iterate through buttons
        for (i = 0; i < user_hide_list.length; i++) {
          // If button matches pattern (regex), disable.
          if (btn_text.match(user_hide_list[i]) && !exempt_list.includes(btn_text)) {
            $(value).prop("disabled", true);
            show = false;
            break;
          }
        }

        // Show nonmatching ones
        if (show == true) {
          $(value).prop("disabled", false);
        }
      } else {
        // Show specified buttons, disable nonspecified
        if (hide_list.indexOf(btn_text) == -1 || exempt_list.indexOf(btn_text) != -1) {
          $(value).prop("disabled", false);
        } else {
          $(value).prop("disabled", true);
        }
      }
    }
  });
}

////////////////////////
// Button Activator
function btn_activator() {
  // Find htstream report section and find target buttons
  var parent_div = $("#mqc-module-section-htstream");
  var btn_groups = parent_div.find(".hc_switch_group").filter(":not(*[class*=htstream_exempt])");

  // Disable and activate appropriate buttons in the list
  $.each(btn_groups, function (x, value) {
    var btns = $(value).find(".btn");
    var first = true;

    for (i = 0; i < btns.length; i++) {
      var attr = $(btns[i]).attr("disabled");

      if (typeof attr == typeof undefined && first == true) {
        $(btns[i]).addClass("active");
        first = false;
      } else {
        $(btns[i]).removeClass("active");
      }
    }
  });
}

//////////////////////////////////////////////////
// Div and Plot Switches

// ////////////////////////
// // Div Switch
// function htstream_div_switch(ele) {
//   // Get ID of button, get corresponding plot divs.
//   var plot_id = ele.id.split("_btn")[0];
//   var parent_node = $(ele).closest(".htstream_fadein");
//   var plot_div = parent_node.find(".hc-plot");

//   // Change div attributes
//   plot_div.attr("id", plot_id);
//   plot_div.attr("class", "hc-plot not_rendered hc-heatmap");

//   // Plot chart (usually heatmaps)
//   renderPlot(plot_id);
// }

////////////////////////
// Plot Switch
function htstream_plot_switch(ele, target) {
  // Get ID of button and find ID's of corresponding on and off divs.
  var unique_id = ele.id.split("_btn")[0];

  //var btns = ele.parentElement.children;
  $(this).addClass("active");
  $("button[id='" + target + "_btn']").removeClass("active");

  // turn off and on proper divs
  var on = $("#" + unique_id);
  var off = $("#" + target);

  // Hide corresponding divs
  off.css("display", "none");
  on.css("display", "block");

  // Find plot div and plot data
  renderPlot(on.find(".hc-plot").attr("id"));
}

////////////////////////
// Button click action for base by cycle graphs
function hts_btn_click(ele) {
  // change value in button
  const unique_id = ele.id.split("_").pop();
  const plot_id = $("div[id$='" + unique_id + "'][id^='htstream_stats_base_line_'][class*='hc-plot']").attr("id");
  const datasets = mqc_plots[plot_id].datasets;
  var target = -1;

  // find index of correct datset
  for (var i = 0; i < datasets.length; i++) {
    if (datasets[i].label == ele.innerText) {
      target = i;
      break;
    }
  }

  const dropdown_id = "htstream_stats_dropdown_" + unique_id;
  const dropdown = document.getElementById(dropdown_id);
  dropdown.innerHTML = ele.innerText + ' <span class="caret"></span>';

  // get children of buttons
  var children = dropdown.parentElement.nextElementSibling.children;

  // iterate through and click appropriate button
  for (var i = 0; i < children.length; i++) {
    var child = children[i];

    if (child.innerText == ele.innerText) {
      child.click();
    }
  }

  // stolen from the plotting.js from multiqc ;)
  mqc_plots[plot_id].activeDatasetIdx = target;
  renderPlot(plot_id);
}

//////////////////////////////////////////////////
// Hide, Rename, and Highlight Samples Handlers

var global_f_add = ["Base: A", "Base: C", "Base: G", "Base: T", "Base: N", "PE Reads", "SE Reads", "PE Bps", "SE Bps"];

var global_on_colors = [
  "#B62612",
  "#82A7E0",
  "#0B8E0B",
  "#DE7D00",
  "#000000",
  "#1EC2D0",
  "#EA8645",
  "#1EC2D0",
  "#EA8645",
];

var sample_num;

////////////////////////
// Hide
$(document).on("mqc_hidesamples", function (e, f_texts, regex_mode) {
  // Deduce mode
  mode = $(".mqc_hidesamples_showhide:checked").val() == "show" ? "show" : "hide";

  // Initialize hide arrrays
  var hide_list = [];
  var user_hide_list = [];

  // Mode is show, add always on list
  if (mode != "hide") {
    var f_add = global_f_add;
  }

  if (f_texts.length != 0) {
    user_hide_list = f_texts.slice();

    // Iterate through list to find buttons to disable
    for (i = 0; i < f_add.length; i++) {
      f_texts.push(f_add[i]);
    }

    hide_list = f_texts.slice();
  }

  // Call functions to disable certain buttons
  btn_disable(mode, regex_mode, hide_list, user_hide_list);
  btn_activator();
});

////////////////////////
// Highlight
$(document).on("mqc_highlights", function (e, f_texts, f_cols, regex_mode) {
  // Iterate through and add highlights for constant samples and colors
  for (i = 0; i < global_f_add.length; i++) {
    f_texts.push(global_f_add[i]);
    f_cols.push(global_on_colors[i]);
  }
});

////////////////////////
// Rename
$(document).on("mqc_renamesamples", function (e, f_texts, t_texts, regex_mode) {
  // Itereate through and add to list (f_texts)
  for (i = 0; i < global_f_add.length; i++) {
    f_texts.push(global_f_add[i]);
    t_texts.push(global_f_add[i]);
  }
});

//////////////////////////////////////////////////
// Page Load Magic
$("document").ready(function () {
  // parse included div with extra configs
  var data = JSON.parse($("#htstream_config").text());

  if (data["sample_colors"] != null) {
    // Check for added sample coloring configurations
    var samples = Object.keys(data["sample_colors"]);
    var colors = data["sample_colors"];

    sample_num = data["htstream_number_of_samples"];

    // Add html coloring samples through MultiQC recolor infrastructure
    if (samples.length != 0) {
      for (i = 0; i < samples.length; i++) {
        $("#mqc_col_filters").append(
          '<li style="color:' +
            colors[samples[i]] +
            ';" id="' +
            samples[i] +
            '"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="' +
            samples[i] +
            '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
        );
      }

      $("#mqc_cols_apply").click();
    }
  }
});
