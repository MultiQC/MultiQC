/* Javascript for the atropos MultiQC module */

// Set up listeners etc on page load
function add_atropos_listeners(phase, passfails) {
  // Add the pass / warning / fails counts to each of the atropos submodule headings
  $.each(passfails, function (k, vals) {
    let pid = "#atropos_" + k;
    let total = 0;
    let v = { pass: 0, warn: 0, fail: 0 };
    $.each(vals, function (s_name, status) {
      total += 1;
      v[status] += 1;
    });
    let p_bar =
      '<div class="progress atropos_passfail_progress"> \
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
    $(pid).append(p_bar);
  });

  // Create popovers on click
  $(".mqc-section-atropos-" + phase + " .atropos_passfail_progress .progress-bar").mouseover(function () {
    // Does this element already have a popover?
    if ($(this).attr("data-original-title")) {
      return false;
    }
    // Create it
    let pid = $(this).closest("h3").attr("id");
    let k = pid.substr(8);
    //alert(k)
    let vals = passfails[k];
    //alert(vals)
    let passes = $(this).hasClass("progress-bar-success") ? true : false;
    let warns = $(this).hasClass("progress-bar-warning") ? true : false;
    let fails = $(this).hasClass("progress-bar-danger") ? true : false;
    let pclass = "";
    if (passes) {
      pclass = "success";
    }
    if (warns) {
      pclass = "warning";
    }
    if (fails) {
      pclass = "danger";
    }
    let samples = Array();
    $.each(vals, function (s_name, status) {
      if (status === "pass" && passes) {
        samples.push(s_name);
      } else if (status === "warn" && warns) {
        samples.push(s_name);
      } else if (status === "fail" && fails) {
        samples.push(s_name);
      }
    });
    $($(this))
      .popover({
        title: $(this).attr("title"),
        content: samples.sort().join("<br>"),
        html: true,
        trigger: "hover click focus",
        placement: "bottom auto",
        template:
          '<div class="popover popover-' +
          pclass +
          '" role="tooltip"> \
                <div class="arrow"></div>\
                <h3 class="popover-title"></h3>\
                <div class="atropos-popover-intro">\
                    Click bar to fix in place <br>\
                    <a href="#" class="atropos-status-highlight"><span class="glyphicon glyphicon-pushpin"></span> Highlight these samples</a><br>\
                    <a href="#" class="atropos-status-hideothers"><span class="glyphicon glyphicon-eye-close"></span> Show only these samples</a>\
                </div>\
                <div class="popover-content"></div>\
            </div>',
      })
      .popover("show");
  });

  // Listener for Status higlight click
  $(".mqc-section-atropos .atropos_passfail_progress").on("click", ".atropos-status-highlight", function (e) {
    e.preventDefault();
    // Get sample names and highlight colour
    let samples = $(this).parent().parent().find(".popover-content").html().split("<br>");
    let f_col = mqc_colours[mqc_colours_idx];
    // Add sample names to the toolbox
    for (i = 0; i < samples.length; i++) {
      let f_text = samples[i];
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
  $(".mqc-section-atropos .atropos_passfail_progress").on("click", ".atropos-status-hideothers", function (e) {
    e.preventDefault();
    // Get sample names
    let samples = $(this).parent().parent().find(".popover-content").html().split("<br>");
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
    for (let i = 0; i < samples.length; i++) {
      let f_text = samples[i];
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
}

function SeqContentHeatmap(phase, read, data, passfails) {
  this.phase = phase;
  this.read = read;
  this.data = data;
  this.passfails = passfails;
  this.s_height = 10;
  this.num_samples = 0;
  this.sample_names = [];
  this.sample_statuses = [];
  this.labels = [];
  this.c_width = 0;
  this.c_height = 0;
  this.ypos = 0;
  this.max_bp = 0;
  this.current_single_plot = undefined;

  let _this = this;

  // Seq Content heatmap mouse rollover
  $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read).mousemove(function (e) {
    // Replace the heading above the heatmap
    let pos = findPos(this);
    let x = e.pageX - pos.x;
    let y = e.pageY - pos.y;
    // Get label from y position
    let idx = Math.floor(y / _this.s_height);
    let s_name = _this.sample_names[idx];
    if (s_name === undefined) {
      return false;
    }
    let s_status = _this.sample_statuses[s_name];
    let s_status_class = "label-default";
    if (s_status === "pass") {
      s_status_class = "label-success";
    }
    if (s_status === "warn") {
      s_status_class = "label-warning";
    }
    if (s_status === "fail") {
      s_status_class = "label-danger";
    }
    $("#atropos_bases_plot_" + _this.phase + "_" + _this.read + " .s_name").html(
      s_name + ' <span class="label s_status ' + s_status_class + '">' + s_status + "</span>",
    );

    // Show the sequence base percentages on the bar plots below
    // http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
    let ctx = this.getContext("2d");
    let p = ctx.getImageData(x, y, 1, 1).data;
    let seq_t = (p[0] / 255) * 100;
    let seq_a = (p[1] / 255) * 100;
    let seq_c = (p[2] / 255) * 100;
    let seq_g = 100 - (seq_t + seq_a + seq_c);
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_t span").text(seq_t.toFixed(0) + "%");
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_c span").text(seq_c.toFixed(0) + "%");
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_a span").text(seq_a.toFixed(0) + "%");
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_g span").text(seq_g.toFixed(0) + "%");

    // Get base pair position from x pos
    let this_bp = Math.floor((x / _this.c_width) * _this.max_bp);
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_pos").text(this_bp + " bp");
  });

  // Remove sample name again when mouse leaves
  $("#atropos_seq_heatmap").mouseout(function (e) {
    $("#atropos_bases_" + _this.phase + "_" + _this.read + "_plot .s_name").html(
      '<em class="text-muted">rollover for sample name</em>',
    );
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_pos").text("-");
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_t span").text("-");
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_c span").text("-");
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_a span").text("-");
    $("#atropos_seq_heatmap_" + _this.phase + "_" + _this.read + "_key_g span").text("-");
  });

  // Highlight the custom heatmap
  $(document).on("mqc_highlights mqc_hidesamples mqc_renamesamples mqc_plotresize", function (e) {
    _this.draw();
  });
  // Seq content - window resized
  $(window).resize(function () {
    _this.draw();
  });

  this.draw = function () {
    // Get sample names, rename and skip hidden samples
    _this.sample_names = [];
    _this.sample_statuses = [];
    let p_data = {};
    let hidden_samples = 0;
    $.each(_this.data, function (s_name, data) {
      // rename sample names
      let t_status = _this.passfails["bases_" + _this.phase][s_name];
      $.each(window.mqc_rename_f_texts, function (idx, f_text) {
        if (window.mqc_rename_regex_mode) {
          let re = new RegExp(f_text, "g");
          s_name = s_name.replace(re, window.mqc_rename_t_texts[idx]);
        } else {
          s_name = s_name.replace(f_text, window.mqc_rename_t_texts[idx]);
        }
      });
      _this.sample_statuses[s_name] = t_status;
      p_data[s_name] = JSON.parse(JSON.stringify(data)); // clone data

      let hide_sample = false;
      for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
        let f_text = window.mqc_hide_f_texts[i];
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
        _this.sample_names.push(s_name);
      } else {
        hidden_samples += 1;
      }
    });
    _this.num_samples = _this.sample_names.length;
    $(
      "#atropos_seq_heatmap_" +
        phase +
        "_" +
        read +
        "_div .samples-hidden-warning, #atropos_seq_heatmap_" +
        phase +
        "_" +
        read +
        "_div .atropos-heatmap-no-samples",
    ).remove();
    $("#atropos_seq_heatmap_" + phase + "_" + read + "_div .hc-plot-wrapper").show();
    if (_this.num_samples === 0) {
      $("#atropos_seq_heatmap_" + phase + "_" + read + "_div .hc-plot-wrapper").hide();
      $("#atropos_seq_heatmap_" + phase + "_" + read + "_div").prepend(
        '<p class="atropos-heatmap-no-samples text-muted">No samples found.</p>',
      );
    }
    if (hidden_samples > 0) {
      $("#atropos_seq_heatmap_" + phase + "_" + read + "_div").prepend(
        '<div class="samples-hidden-warning alert alert-warning"> \
                <span class="glyphicon glyphicon-info-sign"></span> \
                <strong>Warning:</strong> ' +
          hidden_samples +
          ' samples hidden in toolbox. \
                <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a>\
            </div>',
      );
    }
    if (_this.num_samples === 0) {
      return;
    }

    // Convert the CSS percentage size into pixels
    _this.c_width =
      $("#atropos_seq_heatmap_" + phase + "_" + read)
        .parent()
        .width() - 5; // -5 for status bar
    _this.c_height =
      $("#atropos_seq_heatmap_" + phase + "_" + read)
        .parent()
        .height() - 2; // -2 for bottom line padding
    _this.s_height = _this.c_height / _this.num_samples;
    // Minimum row height
    if (_this.s_height < 2) {
      _this.s_height = 2;
      _this.c_height = _this.num_samples * 2;
      $("#atropos_seq_heatmap_" + phase + "_" + read)
        .parent()
        .parent()
        .height(_this.c_height + 10);
    }
    // Resize the canvas properties
    $("#atropos_seq_heatmap_" + phase + "_" + read).prop({
      width: _this.c_width,
      height: _this.c_height + 1,
    });
    let canvas = document.getElementById("atropos_seq_heatmap_" + phase + "_" + read);
    if (canvas && canvas.getContext) {
      let ctx = canvas.getContext("2d");
      ctx.strokeStyle = "#666666";
      // First, do labels and get max base pairs
      _this.max_bp = 0;
      _this.labels = [];
      $.each(_this.sample_names, function (idx, s_name) {
        let s = p_data[s_name];
        _this.labels.push(s_name);
        $.each(s, function (bp, v) {
          bp = parseInt(bp);
          if (bp > _this.max_bp) {
            _this.max_bp = bp;
          }
        });
      });
      _this.ypos = 0;
      $.each(_this.sample_names, function (idx, s_name) {
        // Add a 5px wide bar indicating either status or Highlight
        let status = _this.sample_statuses[s_name];
        let s_col = "#999999";
        if (status === "pass") {
          s_col = "#5cb85c";
        }
        if (status === "warn") {
          s_col = "#f0ad4e";
        }
        if (status === "fail") {
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
        ctx.fillRect(0, _this.ypos + 1, 5, _this.s_height - 2);

        // Plot the squares for the heatmap
        let s = p_data[s_name];
        let xpos = 6;
        let last_bp = 0;
        $.each(s, function (bp, v) {
          bp = parseInt(bp);
          let this_width = (bp - last_bp) * (_this.c_width / _this.max_bp);
          last_bp = bp;
          let t = v["t"] + v["a"] + v["c"] + v["g"];
          v["t"] = (v["t"] / t) * 100;
          v["a"] = (v["a"] / t) * 100;
          v["c"] = (v["c"] / t) * 100;
          v["g"] = (v["g"] / t) * 100;
          let r = Math.round((v["t"] / 100) * 255);
          let g = Math.round((v["a"] / 100) * 255);
          let b = Math.round((v["c"] / 100) * 255);
          ctx.fillStyle = "rgb(" + r + "," + g + "," + b + ")";
          // width+1 to avoid vertical white line gaps.
          ctx.fillRect(xpos, _this.ypos, this_width + 1, _this.s_height);
          xpos += this_width;
        });
        // Draw a line under this row if we don't have too many samples
        if (_this.num_samples <= 20) {
          ctx.beginPath();
          ctx.moveTo(6, _this.ypos);
          ctx.lineTo(_this.c_width, _this.ypos);
          ctx.stroke();
        }
        _this.ypos += _this.s_height;
      });
      // Final line under row
      ctx.beginPath();
      ctx.moveTo(6, _this.ypos);
      ctx.lineTo(_this.c_width, _this.ypos);
      ctx.stroke();
    }
  };
}

// Find the position of the mouse cursor over the canvas
// http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
function findPos(obj) {
  let curleft = 0,
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
