////////////////////////////////////////////////
// Javascript for the FastQC MultiQC module
////////////////////////////////////////////////

///////////////
// Per Base Sequence Content
///////////////

// Global vars
fastqc_passfails = {}; // { <module>: { <section>: { <sample>: { data } } }
fastqc_seq_content = {}; // { <module>: { <sample>: data } }

function load_fastqc_passfails() {
  $(".fastqc_passfails").each(function (i, elem) {
    var key_value = JSON.parse(elem.innerHTML);
    fastqc_passfails[key_value[0]] = key_value[1];
  });
}

function load_fastqc_seq_content() {
  $(".fastqc_seq_content").each(function (i, elem) {
    var key_value = JSON.parse(elem.innerHTML);
    fastqc_seq_content[key_value[0]] = key_value[1];
  });
}

// Set up listeners etc on page load
$(function () {
  load_fastqc_seq_content();

  // Go through each FastQC module in case there are multiple
  // #mqc-module-section-fastqc, #mqc-module-section-fastqc-1, ...
  // or #mqc-module-section-configured-anchor, #mqc-module-section-configured-anchor-1, ...
  var fastqc_modules = $(".fastqc_passfails").closest(".mqc-module-section");
  fastqc_modules.each(function () {
    var module_element = $(this);
    var module_key = module_element.attr("id").replace(/-/g, "_").replace("mqc_module_section_", "");
    fastqc_module(module_element, module_key);
  });
});

function fastqc_module(module_element, module_key) {
  // Per-module shared vars
  var s_height = 10;
  var num_samples = 0;
  var sample_names = [];
  var sample_statuses = [];
  var labels = [];
  var c_width = 0;
  var c_height = 0;
  var ypos = 0;
  var max_bp = 0;
  var current_single_plot = undefined;

  // Make a lookup hash of sample names, in case we rename stuff later
  module_element;
  orig_s_names = {};
  for (var s_name in fastqc_seq_content[module_key]) {
    if (Object.prototype.hasOwnProperty.call(fastqc_seq_content[module_key], s_name)) {
      orig_s_names[s_name] = s_name;
    }
  }

  // Function to plot heatmap
  function fastqc_seq_content_heatmap() {
    // Get sample names, rename and skip hidden samples
    sample_names = [];
    sample_statuses = [];
    var p_data = {};
    var hidden_samples = 0;
    $.each(fastqc_seq_content[module_key], function (s_name, data) {
      // rename sample names
      var orig_s_name = s_name;
      var t_status = fastqc_passfails[module_key]["per_base_sequence_content"][s_name];
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        if (window.mqc_rename_regex_mode) {
          var re = new RegExp(f_text, "g");
          s_name = s_name.replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          s_name = s_name.replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
      orig_s_names[s_name] = orig_s_name;
      sample_statuses[s_name] = t_status;
      p_data[s_name] = JSON.parse(JSON.stringify(data)); // clone data

      var hide_sample = false;
      for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
        var f_text = window.mqc_hide_f_texts[i];
        if (window.mqc_hide_regex_mode) {
          if (s_name.match(f_text)) {
            hide_sample = true;
          }
        } else {
          if (s_name.indexOf(f_text) > -1) {
            hide_sample = true;
          }
        }
      }
      if (window.mqc_hide_mode == "show") {
        hide_sample = !hide_sample;
      }
      if (!hide_sample) {
        sample_names.push(s_name);
      } else {
        hidden_samples += 1;
      }
    });
    num_samples = sample_names.length;
    module_element
      .find("#fastqc_seq_heatmap_div .samples-hidden-warning, #fastqc_seq_heatmap_div .fastqc-heatmap-no-samples")
      .remove();
    module_element.find("#fastqc_seq_heatmap_div .hc-plot-wrapper").show();
    if (num_samples == 0) {
      module_element.find("#fastqc_seq_heatmap_div .hc-plot-wrapper").hide();
      module_element
        .find("#fastqc_seq_heatmap_div")
        .prepend('<p class="fastqc-heatmap-no-samples text-muted">No samples found.</p>');
    }
    if (hidden_samples > 0) {
      module_element.find("#fastqc_seq_heatmap_div").prepend(
        '<div class="samples-hidden-warning alert alert-warning"> \
                <span class="glyphicon glyphicon-info-sign"></span> \
                <strong>Warning:</strong> ' +
          hidden_samples +
          ' samples hidden in toolbox. \
                <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a>\
            </div>',
      );
    }
    if (num_samples == 0) {
      return;
    }

    // Convert the CSS percentage size into pixels
    c_width = module_element.find("#fastqc_seq_heatmap").parent().width() - 5; // -5 for status bar
    c_height = module_element.find("#fastqc_seq_heatmap").parent().height() - 2; // -2 for bottom line padding
    s_height = c_height / num_samples;
    // Minimum row height
    if (s_height < 2) {
      s_height = 2;
      c_height = num_samples * 2;
      module_element
        .find("#fastqc_seq_heatmap")
        .parent()
        .parent()
        .height(c_height + 10);
    }
    // Resize the canvas properties
    module_element.find("#fastqc_seq_heatmap").prop({
      width: c_width,
      height: c_height + 1,
    });
    var canvas = module_element.find("#fastqc_seq_heatmap")[0];
    if (canvas && canvas.getContext) {
      var ctx = canvas.getContext("2d");
      ctx.strokeStyle = "#666666";
      // First, do labels and get max base pairs
      max_bp = 0;
      labels = [];
      $.each(sample_names, function (idx, s_name) {
        var s = p_data[s_name];
        labels.push(s_name);
        $.each(s, function (bp, v) {
          bp = parseInt(bp);
          if (bp > max_bp) {
            max_bp = bp;
          }
        });
      });
      ypos = 0;
      $.each(sample_names, function (idx, s_name) {
        // Add a 5px wide bar indicating either status or Highlight
        var status = sample_statuses[s_name];
        var s_col = "#999999";
        if (status == "pass") {
          s_col = "#5cb85c";
        }
        if (status == "warn") {
          s_col = "#f0ad4e";
        }
        if (status == "fail") {
          s_col = "#d9534f";
        }
        // Override status colour with highlights
        $.each(window.mqc_highlight_f_texts, function (idx, f_text) {
          if (
            (window.mqc_highlight_regex_mode && s_name.match(f_text)) ||
            (!window.mqc_highlight_regex_mode && s_name.indexOf(f_text) > -1)
          ) {
            s_col = window.mqc_highlight_f_cols[idx];
          }
        });
        ctx.fillStyle = s_col;
        ctx.fillRect(0, ypos + 1, 5, s_height - 2);

        // Plot the squares for the heatmap
        var s = p_data[s_name];
        var xpos = 6;
        var last_bp = 0;
        $.each(s, function (bp, v) {
          bp = parseInt(bp);
          var this_width = (bp - last_bp) * (c_width / max_bp);
          last_bp = bp;
          // Very old versions of FastQC give counts instead of percentages
          if (v["t"] > 100) {
            var t = v["t"] + v["a"] + v["c"] + v["g"];
            v["t"] = (v["t"] / t) * 100;
            v["a"] = (v["a"] / t) * 100;
            v["c"] = (v["c"] / t) * 100;
            v["g"] = (v["g"] / t) * 100;
          }
          var r = Math.round((v["t"] / 100) * 255);
          var g = Math.round((v["a"] / 100) * 255);
          var b = Math.round((v["c"] / 100) * 255);
          ctx.fillStyle = "rgb(" + r + "," + g + "," + b + ")";
          // width+1 to avoid vertical white line gaps.
          ctx.fillRect(xpos, ypos, this_width + 1, s_height);
          xpos += this_width;
        });
        // Draw a line under this row if we don't have too many samples
        if (num_samples <= 20) {
          ctx.beginPath();
          ctx.moveTo(6, ypos);
          ctx.lineTo(c_width, ypos);
          ctx.stroke();
        }
        ypos += s_height;
      });
      // Final line under row
      ctx.beginPath();
      ctx.moveTo(6, ypos);
      ctx.lineTo(c_width, ypos);
      ctx.stroke();
    }
  }

  // Draw sequence content heatmap
  fastqc_seq_content_heatmap();

  // Add the pass / warning / fails counts to each of the FastQC submodule headings
  $.each(fastqc_passfails[module_key], function (k, vals) {
    var total = 0;
    var v = { pass: 0, warn: 0, fail: 0 };
    $.each(vals, function (s_name, status) {
      total += 1;
      v[status] += 1;
    });
    var p_bar =
      '<div class="progress fastqc_passfail_progress"> \
            <div class="progress-bar progress-bar-success" style="width: ' +
      (v["pass"] / total) * 100 +
      '%" title="' +
      v["pass"] +
      "&nbsp;/&nbsp;" +
      total +
      ' samples passed">' +
      v["pass"] +
      '</div> \
            <div class="progress-bar progress-bar-warning" style="width: ' +
      (v["warn"] / total) * 100 +
      '%" title="' +
      v["warn"] +
      "&nbsp;/&nbsp;" +
      total +
      ' samples with warnings">' +
      v["warn"] +
      '</div> \
            <div class="progress-bar progress-bar-danger" style="width: ' +
      (v["fail"] / total) * 100 +
      '%" title="' +
      v["fail"] +
      "&nbsp;/&nbsp;" +
      total +
      ' samples failed">' +
      v["fail"] +
      "</div> \
        </div>";
    module_element
      .find("[id^=fastqc_" + k + "]")
      .first()
      .append(p_bar);
  });

  // Create popovers on click
  module_element.find(".fastqc_passfail_progress .progress-bar").mouseover(function () {
    // Does this element already have a popover?
    if ($(this).attr("data-original-title")) {
      return false;
    }
    // Create it
    var pid = $(this).closest("h3").attr("id");
    var k = pid.substr(7);
    // Remove suffix when there are multiple fastqc sections
    var n = k.indexOf("-");
    k = k.substring(0, n != -1 ? n : k.length);
    var vals = fastqc_passfails[module_key][k];
    var passes = $(this).hasClass("progress-bar-success") ? true : false;
    var warns = $(this).hasClass("progress-bar-warning") ? true : false;
    var fails = $(this).hasClass("progress-bar-danger") ? true : false;
    var pclass = "";
    if (passes) {
      pclass = "success";
    }
    if (warns) {
      pclass = "warning";
    }
    if (fails) {
      pclass = "danger";
    }
    var samples = Array();
    $.each(vals, function (s_name, status) {
      if (status == "pass" && passes) {
        samples.push(s_name);
      } else if (status == "warn" && warns) {
        samples.push(s_name);
      } else if (status == "fail" && fails) {
        samples.push(s_name);
      }
    });
    $(this)
      .popover({
        title: $(this).attr("title"),
        content: samples.sort().join("<br>"),
        html: true,
        trigger: "hover click focus",
        placement: "bottom auto",
        template:
          '<div class="popover popover-fastqc-status popover-' +
          pclass +
          '" role="tooltip"> \
                <div class="arrow"></div>\
                <h3 class="popover-title"></h3>\
                <div class="fastqc-popover-intro">\
                    Click bar to fix in place <br>\
                    <a href="#" class="fastqc-status-highlight"><span class="glyphicon glyphicon-pushpin"></span> Highlight these samples</a><br>\
                    <a href="#" class="fastqc-status-hideothers"><span class="glyphicon glyphicon-eye-close"></span> Show only these samples</a>\
                </div>\
                <div class="popover-content"></div>\
            </div>',
      })
      .popover("show");
  });

  // Listener for Status highlight click
  module_element.find(".fastqc_passfail_progress").on("click", ".fastqc-status-highlight", function (e) {
    e.preventDefault();
    // Get sample names and highlight colour
    var samples = $(this).parent().parent().find(".popover-content").html().split("<br>");
    var f_col = $("#mqc_colour_filter_color").val();
    // Add sample names to the toolbox
    for (i = 0; i < samples.length; i++) {
      var f_text = samples[i];
      $("#mqc_col_filters").append(
        '<li style="color:' +
          f_col +
          ';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="' +
          f_text +
          '"/><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );
    }
    // Apply highlights and open toolbox
    apply_mqc_highlights();
    mqc_toolbox_openclose("#mqc_cols", true);
    // Update next highlight colour
    mqc_colours_idx += 1;
    if (mqc_colours_idx >= mqc_colours.length) {
      mqc_colours_idx = 0;
    }
    $("#mqc_colour_filter_color").val(mqc_colours[mqc_colours_idx]);
    // Hide the popover
    $(this).closest(".popover").popover("hide");
  });

  // Listener for Status hide others click
  module_element.find(".fastqc_passfail_progress").on("click", ".fastqc-status-hideothers", function (e) {
    e.preventDefault();
    // Get sample names
    var samples = $(this).parent().parent().find(".popover-content").html().split("<br>");
    // Check if we're already hiding anything, remove after confirm if so
    if ($("#mqc_hidesamples_filters li").length > 0) {
      if (!confirm($("#mqc_hidesamples_filters li").length + " Hide filters already exist - discard?")) {
        return false;
      } else {
        $("#mqc_hidesamples_filters").empty();
      }
    }
    // Set to "show only" and disable regex
    $('.mqc_hidesamples_showhide[value="show"]').prop("checked", true);
    $("#mqc_hidesamples .mqc_regex_mode .re_mode").removeClass("on").addClass("off").text("off");
    // Add sample names to the toolbox
    for (i = 0; i < samples.length; i++) {
      var f_text = samples[i];
      $("#mqc_hidesamples_filters").append(
        '<li><input class="f_text" value="' +
          f_text +
          '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );
    }
    // Apply highlights and open toolbox
    apply_mqc_hidesamples();
    mqc_toolbox_openclose("#mqc_hidesamples", true);
    // Hide the popover
    $(this).closest(".popover").popover("hide");
  });

  /////////
  /// SEQ CONTENT HEATMAP LISTENERS
  /////////

  // Seq Content heatmap export button
  module_element.find("#fastqc_per_base_sequence_content_export_btn").click(function (e) {
    e.preventDefault();
    // In case of repeated modules: #fastqc_per_base_sequence_content_plot, #fastqc_per_base_sequence_content_plot-1, ..
    var plot_id = module_element.find(".fastqc_per_base_sequence_content_plot").attr("id");
    // Tick only this plot in the toolbox and slide out
    $("#mqc_export_selectplots input").prop("checked", false);
    $('#mqc_export_selectplots input[value="' + plot_id + '"]').prop("checked", true);
    mqc_toolbox_openclose("#mqc_exportplots", true);
  });

  // Export plot
  module_element.find(".fastqc_per_base_sequence_content_plot").on("mqc_plotexport_image", function (e, cfg) {
    alert(
      "Apologies, it's not yet possible to export the FastQC per-base sequence content plot.\nPlease take a screengrab or export the JSON data.",
    );
  });
  module_element.find(".fastqc_per_base_sequence_content_plot").on("mqc_plotexport_data", function (e, cfg) {
    if (cfg["ft"] == "json") {
      json_str = JSON.stringify(fastqc_seq_content[module_key], null, 2);
      var blob = new Blob([json_str], { type: "text/plain;charset=utf-8" });
      saveAs(blob, cfg["fname"]);
    } else {
      alert("Apologies, the FastQC per-base sequence content plot can only be exported as JSON currently.");
    }
  });

  // Seq Content heatmap mouse rollover
  module_element.find("#fastqc_seq_heatmap").mousemove(function (e) {
    // Replace the heading above the heatmap
    var pos = findPos(this);
    var x = e.pageX - pos.x + 3;
    var y = e.pageY - pos.y;

    // Get label from y position
    var idx = Math.floor(y / s_height);
    var s_name = sample_names[idx];
    var orig_s_name = orig_s_names[sample_names[idx]];
    if (s_name === undefined) {
      return false;
    }

    // Show the pass/warn/fail status heading for this sample
    var s_status = sample_statuses[s_name];
    var s_status_class = "label-default";
    if (s_status == "pass") {
      s_status_class = "label-success";
    }
    if (s_status == "warn") {
      s_status_class = "label-warning";
    }
    if (s_status == "fail") {
      s_status_class = "label-danger";
    }
    module_element
      .find("#fastqc_per_base_sequence_content_plot_div .s_name")
      .html(
        '<span class="glyphicon glyphicon-info-sign"></span> ' +
          s_name +
          ' <span class="label s_status ' +
          s_status_class +
          '">' +
          s_status +
          "</span>",
      );

    // Update the key with the raw data for this position
    var hover_bp = Math.max(1, Math.floor((x / c_width) * max_bp));
    var thispoint = fastqc_seq_content[module_key][orig_s_name][hover_bp];
    if (!thispoint) {
      var nearestkey = 0;
      var guessdata = null;
      $.each(fastqc_seq_content[module_key][orig_s_name], function (bp, v) {
        bp = parseInt(bp);
        if (bp < hover_bp && bp > nearestkey) {
          nearestkey = bp;
          guessdata = v;
        }
      });
      if (guessdata === null) {
        console.error("Couldn't guess key for " + hover_bp);
        return false;
      } else {
        thispoint = guessdata;
      }
    }
    module_element.find("#fastqc_seq_heatmap_key_t span").text(thispoint["t"].toFixed(0) + "%");
    module_element.find("#fastqc_seq_heatmap_key_c span").text(thispoint["c"].toFixed(0) + "%");
    module_element.find("#fastqc_seq_heatmap_key_a span").text(thispoint["a"].toFixed(0) + "%");
    module_element.find("#fastqc_seq_heatmap_key_g span").text(thispoint["g"].toFixed(0) + "%");
    module_element.find("#fastqc_seq_heatmap_key_pos").text(thispoint["base"] + " bp");
  });

  // Remove sample name again when mouse leaves
  module_element.find("#fastqc_seq_heatmap").mouseout(function (e) {
    module_element
      .find("#fastqc_per_base_sequence_content_plot_div .s_name")
      .html('<span class="glyphicon glyphicon-info-sign"></span> Rollover for sample name');
    module_element.find("#fastqc_seq_heatmap_key_pos").text("-");
    module_element.find("#fastqc_seq_heatmap_key_t span").text("-");
    module_element.find("#fastqc_seq_heatmap_key_c span").text("-");
    module_element.find("#fastqc_seq_heatmap_key_a span").text("-");
    module_element.find("#fastqc_seq_heatmap_key_g span").text("-");
  });

  // Click sample
  module_element.find("#fastqc_seq_heatmap").click(function (e) {
    e.preventDefault();
    // Get label from y position
    var pos = findPos(this);
    var x = e.pageX - pos.x;
    var y = e.pageY - pos.y;
    var idx = Math.floor(y / s_height);
    var s_name = sample_names[idx];
    var orig_s_name = orig_s_names[sample_names[idx]];
    if (orig_s_name !== undefined) {
      plot_single_seqcontent(s_name);
    }
  });
  module_element.on("click", ".fastqc_seqcontent_single_prevnext", function (e) {
    e.preventDefault();
    // Find next / prev sample name
    var idx = sample_names.indexOf(current_single_plot);
    if ($(this).data("action") === "next") {
      idx++;
    } else {
      idx--;
    }
    if (idx < 0) {
      idx = sample_names.length - 1;
    }
    if (idx >= sample_names.length) {
      idx = 0;
    }
    var s_name = sample_names[idx];
    var orig_s_name = orig_s_names[sample_names[idx]];
    current_single_plot = s_name;
    // Prep the new plot data
    var plot_data = [[], [], [], []];
    var bases = Object.keys(fastqc_seq_content[module_key][orig_s_name]).sort(function (a, b) {
      return a - b;
    });
    for (i = 0; i < bases.length; i++) {
      var base = fastqc_seq_content[module_key][orig_s_name][bases[i]]["base"].toString().split("-");
      base = parseFloat(base[0]);
      plot_data[0].push([base, fastqc_seq_content[module_key][orig_s_name][bases[i]]["t"]]);
      plot_data[1].push([base, fastqc_seq_content[module_key][orig_s_name][bases[i]]["c"]]);
      plot_data[2].push([base, fastqc_seq_content[module_key][orig_s_name][bases[i]]["a"]]);
      plot_data[3].push([base, fastqc_seq_content[module_key][orig_s_name][bases[i]]["g"]]);
    }
    // Update the chart
    plot_single_seqcontent(s_name);
  });
  module_element.on("click", "#fastqc_sequence_content_single_back", function (e) {
    e.preventDefault();
    module_element.find("#fastqc_per_base_sequence_content_plot_div").slideDown();
    module_element.find("#fastqc_sequence_content_single_wrapper").slideUp(function () {
      $(this).remove();
    });
  });

  // Highlight the custom heatmap
  $(document).on("mqc_highlights mqc_hidesamples mqc_renamesamples mqc_plotresize", function (e) {
    fastqc_seq_content_heatmap();
  });
  // Seq content - window resized
  $(window).resize(function () {
    fastqc_seq_content_heatmap();
  });

  function plot_single_seqcontent(s_name) {
    current_single_plot = s_name;
    var orig_s_name = orig_s_names[s_name];
    var data = fastqc_seq_content[module_key][orig_s_name];
    var plot_data = [
      { name: "% T", data: [] },
      { name: "% C", data: [] },
      { name: "% A", data: [] },
      { name: "% G", data: [] },
    ];
    var bases = Object.keys(data).sort(function (a, b) {
      return a - b;
    });
    for (i = 0; i < bases.length; i++) {
      var d = bases[i];
      var base = data[d]["base"].toString().split("-");
      base = parseFloat(base[0]);
      plot_data[0]["data"].push({ x: base, y: data[d]["t"], name: data[d]["base"] });
      plot_data[1]["data"].push({ x: base, y: data[d]["c"], name: data[d]["base"] });
      plot_data[2]["data"].push({ x: base, y: data[d]["a"], name: data[d]["base"] });
      plot_data[3]["data"].push({ x: base, y: data[d]["g"], name: data[d]["base"] });
    }

    // Create plot div if it doesn't exist, and hide overview
    if (module_element.find("#fastqc_sequence_content_single_wrapper").length === 0) {
      var plot_div = module_element.find("#fastqc_per_base_sequence_content_plot_div");
      plot_div.slideUp();
      var newplot =
        '<div id="fastqc_sequence_content_single_wrapper"> \
            <div id="fastqc_sequence_content_single_controls">\
                <button class="btn btn-primary btn-sm" id="fastqc_sequence_content_single_back">Back to overview heatmap</button> \
                <div class="btn-group btn-group-sm"> \
                    <button class="btn btn-default fastqc_seqcontent_single_prevnext" data-action="prev">&laquo; Prev</button> \
                    <button class="btn btn-default fastqc_seqcontent_single_prevnext" data-action="next">Next &raquo;</button> \
                </div>\
            </div>\
            <div class="hc-plot-wrapper"><div id="fastqc_sequence_content_single" class="hc-plot hc-line-plot"><small>loading..</small></div></div></div>';
      $(newplot).insertAfter(plot_div).hide().slideDown();
    }

    let target = "fastqc_sequence_content_single";
    let traces = plot_data.map((d) => {
      return {
        type: "line",
        x: d["data"].map((val) => val.x),
        y: d["data"].map((val) => val.y),
        mode: "lines",
        name: d["name"],
        hovertemplate: "%{y:.1f}%",
      };
    });
    let layout = {
      title: s_name,
      colorway: ["#dc0000", "#0000dc", "#00dc00", "#404040"],
      xaxis: {
        title: "Position",
        ticksuffix: " bp",
      },
      yaxis: {
        title: "% Reads",
        range: [0, 100],
        ticksuffix: "%",
      },
      hovermode: "x unified",
    };
    let config = {
      responsive: true,
      displaylogo: false,
      displayModeBar: true,
      toImageButtonOptions: { filename: target },
      modeBarButtonsToRemove: [
        "lasso2d",
        "autoScale2d",
        "pan2d",
        "select2d",
        "zoom2d",
        "zoomIn2d",
        "zoomOut2d",
        "resetScale2d",
        "toImage",
      ],
    };
    Plotly.newPlot(target, traces, layout, config);
  }
}

// Find the position of the mouse cursor over the canvas
// http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
function findPos(obj) {
  var curleft = 0,
    curtop = 0;
  if (obj.offsetParent) {
    do {
      curleft += obj.offsetLeft;
      curtop += obj.offsetTop;
    } while ((obj = obj.offsetParent));
    return { x: curleft, y: curtop };
  }
  return undefined;
}
