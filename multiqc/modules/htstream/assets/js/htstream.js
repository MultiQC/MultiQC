// Functions

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

////////////////////////
// Plot Switch
function htstream_plot_switch(ele, target) {
  // Get ID of button and find ID's of corresponding on and off divs.
  const unique_id = ele.id.split("_btn")[0];
  $("#" + ele.id).addClass("active");
  $("button[id='" + target + "_btn']").removeClass("active"); // remove other button's active status

  // turn off and on proper divs
  let on = $("#" + unique_id);
  let off = $("#" + target);

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
  let target = "";

  // find index of correct datset
  for (var i = 0; i < datasets.length; i++) {
    if (datasets[i].label == ele.innerText) {
      target = i;
      break;
    }
  }

  if (target == null) {
    console.error("Sample " + ele.innerText + " not available.");
  }

  // change text of dropdown menu
  const dropdown = $("#htstream_stats_dropdown_" + unique_id);
  dropdown.text(ele.innerText);
  dropdown.append('<span class="caret"></span>');

  // stolen from the plotting.js from multiqc ;)
  mqc_plots[plot_id].activeDatasetIdx = target;
  renderPlot(plot_id);
}
