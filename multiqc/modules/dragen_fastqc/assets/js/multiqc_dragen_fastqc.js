////////////////////////////////////////////////
// Javascript for the DRAGEN FastQC MultiQC module
////////////////////////////////////////////////

///////////////
// Per Base Sequence Content
///////////////

// Global vars
fastqc_seq_content = {}; // { <module>: { <sample>: data } }

function load_fastqc_seq_content() {
  $(".fastqc_seq_content").each(function (i, elem) {
    var key_value = JSON.parse(elem.innerHTML);
    fastqc_seq_content[key_value[0]] = key_value[1];
  });
}

// Set up listeners etc on page load
callAfterDecompressed.push(function (mqc_plotdata) {
  load_fastqc_seq_content();

  // Go through each DRAGEN-FastQC module in case there are multiple
  var fastqc_modules = $(".fastqc_seq_content").closest(".mqc-module-section");
  fastqc_modules.each(function () {
    var module_element = $(this);
    var module_key = module_element.data("moduleAnchor");
    fastqc_module(module_element, module_key);
  });
});

function fastqc_module(module_element, module_key) {
  // Per-module shared vars
  var s_height = 10;
  var num_samples = 0;
  var sample_names = [];
  var labels = [];
  var c_width = 0;
  var c_height = 0;
  var ypos = 0;
  var max_bp = 0;
  var current_single_plot = undefined;

  // Make a lookup hash of sample names, in case we rename stuff later
  var orig_s_names = {};
  for (var s_name in fastqc_seq_content[module_key]) {
    if (Object.prototype.hasOwnProperty.call(fastqc_seq_content[module_key], s_name)) {
      orig_s_names[s_name] = s_name;
    }
  }

  // Function to plot heatmap
  function fastqc_seq_content_heatmap() {
    // Get sample names, rename and skip hidden samples
    var p_data = {};
    var hidden_samples = 0;
    $.each(fastqc_seq_content[module_key], function (s_name, data) {
      // rename sample names
      var orig_s_name = s_name;
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        if (window.mqc_rename_regex_mode) {
          var re = new RegExp(f_text, "g");
          s_name = s_name.replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          s_name = s_name.replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
      orig_s_names[s_name] = orig_s_name;
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
      if (window.mqc_hide_mode === "show") {
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
    if (num_samples === 0) {
      module_element.find("#fastqc_seq_heatmap_div .hc-plot-wrapper").hide();
      module_element
        .find("#fastqc_seq_heatmap_div")
        .prepend('<p class="fastqc-heatmap-no-samples text-muted">No samples found.</p>');
    }
    if (hidden_samples > 0) {
      module_element.find("#fastqc_seq_heatmap_div").prepend(
        '<div class="samples-hidden-warning alert alert-warning"> \
                ⚠ <strong>Warning:</strong> ' +
          hidden_samples +
          ' samples hidden in toolbox. \
                <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a>\
            </div>',
      );
    }
    if (num_samples === 0) {
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
        // Add a 5px wide bar for highlights
        var s_col = "#999999";
        // Check for highlight colours
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

  /////////
  /// SEQ CONTENT HEATMAP LISTENERS
  /////////

  // Seq Content heatmap export button
  module_element.find("#dragen_fastqc_per_base_sequence_content_export_btn").click(function (e) {
    e.preventDefault();
    // In case of repeated modules: #dragen_fastqc_per_base_sequence_content_plot, #dragen_fastqc_per_base_sequence_content_plot-1, ..
    var plot_id = module_element.find(".dragen_fastqc_per_base_sequence_content_plot").attr("id");
    // Tick only this plot in the toolbox and slide out
    $("#mqc_export_selectplots input").prop("checked", false);
    $('#mqc_export_selectplots input[value="' + plot_id + '"]').prop("checked", true);
    mqc_toolbox_openclose("#mqc_exportplots", true);
  });

  // Export plot
  module_element.find(".dragen_fastqc_per_base_sequence_content_plot").on("mqc_plotexport_image", function (e, cfg) {
    alert(
      "Apologies, it's not yet possible to export the DRAGEN-FastQC per-base sequence content plot.\nPlease take a screengrab or export the JSON data.",
    );
  });
  module_element.find(".dragen_fastqc_per_base_sequence_content_plot").on("mqc_plotexport_data", function (e, cfg) {
    if (cfg["ft"] == "json") {
      json_str = JSON.stringify(fastqc_seq_content[module_key], null, 2);
      var blob = new Blob([json_str], { type: "text/plain;charset=utf-8" });
      saveAs(blob, cfg["fname"]);
    } else {
      alert("Apologies, the DRAGEN-FastQC per-base sequence content plot can only be exported as JSON currently.");
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

    // Show the sample name
    module_element.find("#dragen_fastqc_per_base_sequence_content_plot_div .s_name").html(s_name);

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
    module_element.find("#dragen_fastqc_per_base_sequence_content_plot_div .s_name").html("Rollover for sample name");
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
  module_element.on("click", "#fastqc_sequence_content_single_back", function (e) {
    e.preventDefault();
    module_element.find("#dragen_fastqc_per_base_sequence_content_plot_div").slideDown();
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
    if (module_element.find("#fastqc_sequence_content_single_wrapper").length == 0) {
      var plot_div = module_element.find("#dragen_fastqc_per_base_sequence_content_plot_div");
      plot_div.slideUp();
      var newplot =
        '<div id="fastqc_sequence_content_single_wrapper"> \
            <div id="fastqc_sequence_content_single_controls">\
                <button class="btn btn-primary btn-sm" id="fastqc_sequence_content_single_back">Back to overview heatmap</button> \
                <div class="btn-group btn-group-sm"> \
                    <button class="btn btn-outline-secondary fastqc_seqcontent_single_prevnext" data-action="prev">&laquo; Prev</button> \
                    <button class="btn btn-outline-secondary fastqc_seqcontent_single_prevnext" data-action="next">Next &raquo;</button> \
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
