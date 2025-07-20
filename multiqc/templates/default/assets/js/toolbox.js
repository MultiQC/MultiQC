////////////////////////////////////////////////
// MultiQC Report Toolbox Code
////////////////////////////////////////////////

let mqc_colours_idx = 0;
const mqc_colours = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a9a904", "#a65628", "#f781bf", "#999999"];
const zip_threshold = 8;

// Add these constants at the top of the file
const AI_PROVIDERS = {
  seqera: {
    name: "Seqera AI",
    apiKeysUrl: "https://cloud.seqera.io/tokens",
  },
  anthropic: {
    name: "Anthropic",
    defaultModel: "claude-sonnet-4-0",
    suggestedModels: ["claude-sonnet-4-0"],
    apiKeysUrl: "https://console.anthropic.com/settings/keys",
    modelsUrl: "https://docs.anthropic.com/en/docs/intro-to-claude#model-options",
  },
  openai: {
    name: "OpenAI",
    defaultModel: "gpt-4o",
    suggestedModels: ["gpt-4o", "gpt-4.1"],
    apiKeysUrl: "https://platform.openai.com/api-keys",
    modelsUrl: "https://platform.openai.com/docs/models",
  },
  aws_bedrock: {
    name: "AWS Bedrock",
    modelsUrl: "https://docs.anthropic.com/en/docs/intro-to-claude#model-options",
  },
  custom: {
    name: "Custom",
    defaultModel: "",
  },
  clipboard: {
    name: "Copy prompts",
  },
  none: {
    name: "Remove AI buttons",
  },
};
const AI_PROVIDER_GROUPS = {
  "In-report summaries": ["seqera", "anthropic", "openai", "custom"],
  Alternatives: ["clipboard", "none"],
};

//////////////////////////////////////////////////////
// TOOLBOX LISTENERS
//////////////////////////////////////////////////////
$(function () {
  // Close toolbox when Escape key is pressed
  $(document).keyup(function (e) {
    if (e.key === "Escape" && $(".mqc-toolbox").hasClass("active")) {
      mqc_toolbox_openclose(undefined, false);
    }
  });

  // Batch sample renaming buttons
  $(".mqc_sname_switches").click(function (e) {
    e.preventDefault();
    if ($(this).hasClass("active")) {
      return false;
    }
    $("#mqc_sname_switches button").removeClass("active");
    $(this).addClass("active");
    // Clear previous bulk renaming entries
    $(".mqc_sname_switches_li").remove();
    // Build new renaming entries and apply
    var j = $(this).data("index");
    if (j == 0) {
      apply_mqc_renamesamples();
    } else {
      for (i = 0; i < mqc_config["sample_names_rename"].length; i++) {
        var ft = mqc_config["sample_names_rename"][i][0];
        var tt = mqc_config["sample_names_rename"][i][j];
        $("#mqc_renamesamples_filters").append(
          '<li class="mqc_sname_switches_li"> \
          <input class="f_text from_text" value="' +
            ft +
            '" /><small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="' +
            tt +
            '" /> \
          <button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button> \
        </li>',
        );
      }
      apply_mqc_renamesamples();
    }
  });

  // Show/hide buttons
  $(".mqc_hide_switches").click(function (e) {
    e.preventDefault();
    if ($(this).hasClass("active")) {
      return false;
    }
    $("#mqc_hide_switches button").removeClass("active");
    $(this).addClass("active");

    // Clear previous show/hide group
    $("#mqc_hidesamples_filters").empty();

    // Get requested pattern and whether to show or hide the pattern
    var j = $(this).data("index");
    var pattern = mqc_config["show_hide_patterns"][j];
    var show_hide_mode = mqc_config["show_hide_mode"][j];
    var regex = mqc_config["show_hide_regex"][j];
    if (!Array.isArray(pattern)) {
      pattern = [pattern];
    }
    if (show_hide_mode === undefined) {
      show_hide_mode = "show";
    }

    // click the regex button if we want it turned on/off
    var button = document.getElementsByClassName("mqc_switch re_mode")[2];
    if (button.className.includes(" on") && !regex) {
      button.click();
    }
    if (button.className.includes(" off") && regex) {
      button.click();
    }

    // Apply the changes
    $(".mqc_hidesamples_showhide[value=" + show_hide_mode + "]").prop("checked", true);
    $(pattern).each(function (idx, val) {
      $("#mqc_hidesamples_filters").append(
        '<li><input class="f_text" value="' +
          val +
          '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );
    });
    apply_mqc_hidesamples(show_hide_mode);
  });

  // Hide toolbox when clicking outside
  $(document).mouseup(function (e) {
    if (!$(".mqc-toolbox").is(e.target) && $(".mqc-toolbox").has(e.target).length === 0) {
      if ($(".mqc-toolbox").hasClass("active")) {
        mqc_toolbox_openclose(undefined, false);
      }
    }
  });

  // Hide toolbox when a modal is shown
  $(".modal").on("show.bs.modal", function (e) {
    if ($(".mqc-toolbox").hasClass("active")) {
      mqc_toolbox_openclose(undefined, false);
    }
  });

  // Listener to re-plot graphs if config loaded
  $(document).on("mqc_config_loaded", function (e) {
    $(".hc-plot:not(.not_rendered)").each(function () {
      let target = $(this).attr("id");
      renderPlot(target);
    });
  });

  // Toolbox buttons
  $(".mqc-toolbox-buttons a").click(function (e) {
    e.preventDefault();
    var target = $(this).attr("href");
    mqc_toolbox_openclose(target);
  });

  // Download DOIs
  $(".download-citations-btn").click(function (e) {
    e.preventDefault();
    var format = $(this).data("format");
    var doi_list = { "10.1093/bioinformatics/btw354": "MultiQC" };
    $(".module-doi").each(function () {
      var module_id = $(this).closest(".mqc-module-section-first").find("h2").attr("id");
      doi_list[$(this).data("doi")] = module_id;
    });
    // Get BibTeX
    if (format == "bibtex") {
      var bibtex_string = "";
      // Kick off crossref api calls
      var ajax_promises = [];
      for (var doi in doi_list) {
        ajax_promises.push(
          $.get("https://api.crossref.org/works/" + doi + "/transform/application/x-bibtex", function (data) {
            bibtex_string += data + "\n";
          }),
        );
      }
      // Wait until all API calls are done
      $.when.apply(null, ajax_promises).then(function () {
        var blob = new Blob([bibtex_string], { type: "text/plain;charset=utf-8" });
        saveAs(blob, "multiqc_references.bib");
      });
    }
    // Download list of DOIs
    else {
      var doi_string = "";
      for (var doi in doi_list) {
        doi_string += doi + new Array(50 - doi.length).join(" ") + " # " + doi_list[doi] + "\n";
      }
      var blob = new Blob([doi_string], { type: "text/plain;charset=utf-8" });
      saveAs(blob, "multiqc_dois.txt");
    }
  });

  // Highlight colour filters
  $("#mqc_color_form").submit(function (e) {
    e.preventDefault();
    let mqc_colour_filter = $("#mqc_colour_filter");
    let mqc_colour_filter_color = $("#mqc_colour_filter_color");
    let f_text = mqc_colour_filter.val().trim();
    let f_col = mqc_colour_filter_color.val().trim();
    $("#mqc_col_filters").append(
      '<li style="color:' +
        f_col +
        ';" id="' +
        hashCode(f_text + f_col) +
        '"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="' +
        f_text +
        '" tabindex="' +
        mqc_colours_idx +
        '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
    );
    $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    mqc_colour_filter.val("");
    mqc_colours_idx += 1;
    if (mqc_colours_idx >= mqc_colours.length) {
      mqc_colours_idx = 0;
    }
    mqc_colour_filter_color.val(mqc_colours[mqc_colours_idx]);
  });
  $("#mqc_cols_apply").click(function (e) {
    if (apply_mqc_highlights()) {
      $(this).attr("disabled", true).removeClass("btn-primary").addClass("btn-default");
      mqc_auto_save_config();
    }
  });

  // Rename samples
  let mqc_renamesamples_idx = 300;
  $("#mqc_renamesamples_form").submit(function (event) {
    event.preventDefault();

    let mqc_renamesamples_from = $("#mqc_renamesamples_from");
    let mqc_renamesamples_to = $("#mqc_renamesamples_to");
    let from_text = mqc_renamesamples_from.val().trim();
    let to_text = mqc_renamesamples_to.val().trim();

    if (from_text.length === 0) {
      alert('Error - "From" text must not be blank.');
      return false;
    }

    let li =
      '<li><input class="f_text from_text" value="' + from_text + '" tabindex="' + mqc_renamesamples_idx + '" />';
    li +=
      '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="' +
      to_text +
      '" tabindex="' +
      (mqc_renamesamples_idx + 1) +
      '" />';
    li +=
      '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>';
    $("#mqc_renamesamples_filters").append(li);
    $("#mqc_rename_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");

    // Reset form
    mqc_renamesamples_from.val("");
    mqc_renamesamples_to.val("");
    mqc_renamesamples_idx += 2;
    $("#mqc_renamesamples_form input:first").focus();
  });

  $("#mqc_rename_apply").click(function (e) {
    if (apply_mqc_renamesamples()) {
      $(this).attr("disabled", true).removeClass("btn-primary").addClass("btn-default");
      mqc_auto_save_config();
    }
  });

  // Bulk rename samples
  $("#mqc_renamesamples_bulk_collapse").on("shown.bs.collapse", function () {
    $("#mqc_renamesamples_bulk_form textarea").focus();
  });
  $("#mqc_renamesamples_bulk_form").submit(function (e) {
    e.preventDefault();
    var raw = $(this).find("textarea").val();
    var lines = raw.match(/^.*([\n\r]+|$)/gm);
    $.each(lines, function (i, l) {
      var sections = l.split("\t", 2);
      if (sections.length < 2) {
        return true;
      }
      var from_text = sections[0].trim();
      var to_text = sections[1].trim();
      if (from_text.length == 0) {
        return true;
      }
      var li =
        '<li><input class="f_text from_text" value="' + from_text + '" tabindex="' + mqc_renamesamples_idx + '" />';
      li +=
        '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="' +
        to_text +
        '" tabindex="' +
        (mqc_renamesamples_idx + 1) +
        '" />';
      li +=
        '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>';
      $("#mqc_renamesamples_filters").append(li);
    });
    $("#mqc_rename_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    $(this).find("textarea").val("");
    $("#mqc_renamesamples_bulk_collapse").collapse("hide");
  });

  // Hide sample filters
  var mqc_hidesamples_idx = 200;
  $("#mqc_hidesamples_form").submit(function (e) {
    e.preventDefault();
    var f_text = $("#mqc_hidesamples_filter").val().trim();
    if (f_text.length == 0) {
      alert("Error - filter text must not be blank.");
      return false;
    }
    $("#mqc_hidesamples_filters").append(
      '<li><input class="f_text" value="' +
        f_text +
        '" tabindex="' +
        mqc_hidesamples_idx +
        '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
    );
    $("#mqc_hide_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    $("#mqc_hidesamples_filter").val("");
    mqc_hidesamples_idx += 1;
  });
  $(".mqc_hidesamples_showhide").change(function (e) {
    $("#mqc_hide_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
  });
  $("#mqc_hide_apply").click(function (e) {
    if (apply_mqc_hidesamples()) {
      $(this).attr("disabled", true).removeClass("btn-primary").addClass("btn-default");
      mqc_auto_save_config();
    }
  });

  // EXPORTING PLOTS
  // Change text on download button
  $('#mqc_exportplots a[data-toggle="tab"]').on("shown.bs.tab", function (e) {
    if ($(e.target).attr("href") === "#mqc_data_download") {
      $("#mqc-dl-plot-txt").text("Data");
    } else {
      $("#mqc-dl-plot-txt").text("Images");
    }
  });
  // Load the plot exporter
  if ($(".hc-plot").length > 0) {
    $(".hc-plot").each(function () {
      var fname = $(this).attr("id");
      $("#mqc_export_selectplots").append(
        '<div class="checkbox"><label><input type="checkbox" value="' +
          fname +
          '" checked> ' +
          fname +
          "</label></div>",
      );
    });
    // Select all / none for checkboxes
    $("#mqc_export_sall").click(function (e) {
      e.preventDefault();
      $("#mqc_export_selectplots input").prop("checked", true);
    });
    $("#mqc_export_snone").click(function (e) {
      e.preventDefault();
      $("#mqc_export_selectplots input").prop("checked", false);
    });
    // Aspect ratio fixed
    var mqc_exp_aspect_ratio = $("#mqc_exp_width").val() / $("#mqc_exp_height").val();
    $("#mqc_export_aspratio").change(function () {
      if ($(this).is(":checked")) {
        mqc_exp_aspect_ratio = $("#mqc_exp_width").val() / $("#mqc_exp_height").val();
      }
    });
    $("#mqc_exp_width").keyup(function () {
      if ($("#mqc_export_aspratio").is(":checked")) {
        $("#mqc_exp_height").val($(this).val() / mqc_exp_aspect_ratio);
      }
    });
    $("#mqc_exp_height").keyup(function () {
      if ($("#mqc_export_aspratio").is(":checked")) {
        $("#mqc_exp_width").val($(this).val() * mqc_exp_aspect_ratio);
      }
    });

    function dataUrlToBlob(dataUrl, mime) {
      // Split the data URL at the comma
      const byte_str = atob(dataUrl.split(",")[1]);
      const byte_numbers = new Array(byte_str.length);
      for (let i = 0; i < byte_str.length; i++) {
        byte_numbers[i] = byte_str.charCodeAt(i);
      }
      const byte_array = new Uint8Array(byte_numbers);
      return new Blob([byte_array], { type: mime });
    }

    // Export the plots
    $("#mqc_exportplots").submit(function (e) {
      e.preventDefault();
      let checked_plots = $("#mqc_export_selectplots input:checked");
      let zip = new JSZip();
      let promises = [];
      //////
      ////// EXPORT PLOT IMAGES
      //////
      if ($("#mqc_image_download").is(":visible")) {
        let mime = $("#mqc_export_ft").val();
        let format = mime.replace("image/", "").split("+")[0];
        let f_width = parseInt($("#mqc_exp_width").val());
        let f_height = parseInt($("#mqc_exp_height").val());
        const font_scale = parseFloat($("#mqc_export_scaling").val());
        checked_plots.each(function () {
          const target = $(this).val();

          promises.push(
            Plotly.toImage(target, {
              format: format,
              width: f_width / font_scale,
              height: f_height / font_scale,
              scale: font_scale,
            }).then(function (img) {
              if (format === "svg") {
                Plotly.Snapshot.downloadImage(target, {
                  format: format,
                  width: f_width / font_scale,
                  height: f_height / font_scale,
                  scale: font_scale,
                  filename: target,
                });
                // if (checked_plots.length <= zip_threshold) {
                //   // Not many plots to export, just trigger a download for each
                //   const data = img.replace(/^data:image\/svg\+xml;base64,/, "");
                //   const blob = new Blob([data], { type: mime });
                //   saveAs(blob, target + "." + format);
                // } else {
                //   // Lots of plots - add to a zip file for download
                //   const fname = target + "." + format;
                //   zip.file(fname, img, { base64: false });
                // }
              } else {
                // Can add logo to a PNG image
                addLogo(img, function (imageWithLogo) {
                  if (checked_plots.length <= zip_threshold) {
                    // Not many plots to export, just trigger a download for each:"
                    const blob = dataUrlToBlob(imageWithLogo, mime);
                    saveAs(blob, target + "." + format);
                  } else {
                    // Lots of plots - add to a zip file for download:
                    const fname = target + "." + format;
                    const data = imageWithLogo.replace(/^data:image\/png;base64,/, "");
                    zip.file(fname, data, { base64: true });
                  }
                });
              }
            }),
          );
        });
        if (checked_plots.length > zip_threshold) {
          // Wait for all promises to resolve
          Promise.all(promises).then(() => {
            zip.generateAsync({ type: "blob" }).then(function (content) {
              saveAs(content, "multiqc_plots.zip");
            });
          });
        }
      }
      //////
      ////// EXPORT PLOT DATA
      //////
      else if ($("#mqc_data_download").is(":visible")) {
        const format = $("#mqc_export_data_ft").val();
        console.log("Exporting data in " + format + " format");
        let skipped_plots = 0;
        checked_plots.each(function () {
          try {
            const target = $(this).val();
            const fname = target + "." + format;
            // If JSON then just dump everything
            if (format === "json") {
              const json_str = JSON.stringify(mqc_plots[target], null, 2);
              const blob = new Blob([json_str], { type: "text/plain;charset=utf-8" });
              if (checked_plots.length <= zip_threshold) {
                // Not many plots to export, just trigger a download for each
                saveAs(blob, fname);
              } else {
                // Lots of plots - add to a zip file for download
                zip.file(fname, blob);
              }
            } else if (format === "tsv" || format === "csv") {
              let plot = mqc_plots[target];
              if (plot !== undefined) {
                let text = plot.exportData(format);
                const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
                if (checked_plots.length <= zip_threshold) {
                  // Not many plots to export, just trigger a download for each
                  saveAs(blob, fname);
                } else {
                  // Lots of plots - generate a zip file for download.
                  // Add to a zip archive
                  zip.file(fname, blob);
                }
              } else {
                skipped_plots += 1;
              }
            } else {
              skipped_plots += 1;
            }
          } catch (e) {
            console.error(e);
            skipped_plots += 1;
          }
        });
        if (skipped_plots > 0) {
          alert("Warning: Could not export data from " + skipped_plots + " plots.");
        }
        // Save the zip and trigger a download
        if (checked_plots.length > zip_threshold) {
          zip.generateAsync({ type: "blob" }).then(function (content) {
            saveAs(content, "multiqc_data.zip");
          });
        }
      } else {
        alert("Error - don't know what to export!");
      }
    });
  } else {
    $("#mqc_exportplots").hide();
    $(".mqc-toolbox-buttons a[href=#mqc_exportplots]").parent().hide();
  }

  // Export plot buttons
  $(".export-plot").click(function (e) {
    e.preventDefault();
    // Get the id of the span element that was clicked
    let id = e.target.dataset.plotAnchor;
    let isTable = e.target.dataset.type === "table";
    // Tick only this plot in the toolbox and slide out
    $("#mqc_export_selectplots input").prop("checked", false);
    $('#mqc_export_selectplots input[value="' + id + '"]').prop("checked", true);
    // Special case - Table scatter plots are in a modal, need to close this first
    if (id === "table_scatter_plot") {
      $("#table_scatter_modal").modal("hide");
    }
    mqc_toolbox_openclose(
      "#mqc_exportplots",
      true,
      isTable, // no image export for table, go directly to data download
    );
  });

  /// SAVING STUFF
  // Load the saved setting names
  populate_mqc_saveselect();
  // Save config
  $("#mqc_saveconfig_form").submit(function (e) {
    e.preventDefault();
    var name = $(this).find("input").val().trim();
    if (name == "") {
      alert("Error - you must name the saved settings.");
    } else {
      mqc_save_config(name);
    }
  });
  // Load config
  $("#mqc_loadconfig_form").submit(function (e) {
    e.preventDefault();
    var name = $(this).find("select").val().trim();
    if (name == "") {
      alert("Error - No saved setting selected.");
    } else {
      load_mqc_config(name);
    }
  });
  // Delete config
  $(".mqc_config_clear").click(function (e) {
    e.preventDefault();
    var name = $("#mqc_loadconfig_form select").val().trim();
    if (name == "") {
      alert("Error - no saved settings selected.");
    } else {
      if (confirm("Delete saved settings '" + name + "'?")) {
        mqc_save_config(name, true);
      }
    }
  });
  // Set current config as default
  $(".mqc_config_set_default").click(function (e) {
    e.preventDefault();
    var name = $("#mqc_loadconfig_form select").val().trim();
    if (name == "") {
      alert("Error - no saved settings selected.");
    } else {
      load_mqc_config(name);
      mqc_save_config(name, false, true);
    }
  });
  // Clear current config default
  $(".mqc_config_clear_default").click(function (e) {
    e.preventDefault();
    mqc_clear_default_config();
  });

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

  // Use jQuery UI to make the colour filters sortable
  $("#mqc_col_filters").sortable();
  $("#mqc_col_filters").on("sortstop", function (event, ui) {
    $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
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

  /////////////////////////
  // REGEX HELP MODAL
  /////////////////////////
  $(".regex_example_buttons button").click(function (e) {
    e.preventDefault();
    $(".regex_example_demo input").val($(this).data("example"));
    regex_example_test();
  });
  $(".regex_example_demo input").keyup(function (e) {
    regex_example_test();
  });
  function regex_example_test() {
    var re = $(".regex_example_demo input").val();
    console.log("Testing " + re);
    $(".regex_example_demo pre span").each(function () {
      $(this).removeClass();
      if ($(this).text().match(re)) {
        console.log("Matches " + $(this).text());
        $(this).addClass("mark text-success");
      } else {
        console.log("Matches " + $(this).text());
        $(this).addClass("text-muted");
      }
    });
  }

  // Load auto-saved config on page load
  try {
    const autoSaveKey = AUTO_SAVE_PREFIX + window.reportUuid;
    const local_config = localStorage.getItem("mqc_config");
    if (local_config) {
      const configs = JSON.parse(local_config);
      if (configs[autoSaveKey]) {
        load_mqc_config(autoSaveKey);
      }
    }
  } catch (e) {
    console.log("Could not load auto-saved config:", e);
  }

  // Check if we have highlight patterns from local storage
  let has_highlight_filters = $("#mqc_col_filters").children().length > 0;
  let has_hide_filters = $("#mqc_hidesamples_filters").children().length > 0;
  let has_rename_filters = $("#mqc_renamesamples_filters").children().length > 0;

  // Apply pre-configured highlight patterns from config only if no local storage values
  if (!has_highlight_filters && mqc_config.highlight_patterns && mqc_config.highlight_patterns.length > 0) {
    // Set regex mode if specified
    if (mqc_config.highlight_regex) {
      $("#mqc_cols .mqc_regex_mode .re_mode").removeClass("off").addClass("on").text("on");
    }

    // Add each pattern with its color
    for (let i = 0; i < mqc_config.highlight_patterns.length; i++) {
      const pattern = mqc_config.highlight_patterns[i];
      let color =
        mqc_config.highlight_colors && mqc_config.highlight_colors[i] ? mqc_config.highlight_colors[i] : "#e41a1c";

      // Add to the filters list
      $("#mqc_col_filters").append(
        '<li style="color:' +
          color +
          '"><span class="hc_handle"><span></span><span></span><span></span></span><input class="f_text" value="' +
          pattern +
          '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );
    }

    // Apply the highlights
    apply_mqc_highlights();
  }

  // Apply pre-configured hide samples from config only if no local storage values
  if (!has_hide_filters && mqc_config.show_hide_patterns && mqc_config.show_hide_patterns.length > 0) {
    // Add each pattern
    for (let i = 1; i < mqc_config.show_hide_patterns.length; i++) {
      // Skip first (Show all)
      const pattern = mqc_config.show_hide_patterns[i];
      $("#mqc_hidesamples_filters").append(
        '<li><input class="f_text" value="' +
          pattern +
          '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );
    }

    // Set regex mode if specified for the first non-empty pattern
    for (let i = 1; i < mqc_config.show_hide_regex.length; i++) {
      if (mqc_config.show_hide_regex[i]) {
        $("#mqc_hidesamples .mqc_regex_mode .re_mode").removeClass("off").addClass("on").text("on");
        break;
      }
    }

    // Set show/hide mode based on the first non-empty pattern
    for (let i = 1; i < mqc_config.show_hide_mode.length; i++) {
      if (mqc_config.show_hide_mode[i] === "show") {
        $(".mqc_hidesamples_showhide[value=show]").prop("checked", true);
        break;
      }
    }

    // Apply the hide/show
    apply_mqc_hidesamples();
  }

  // Apply pre-configured sample renaming patterns from config only if no local storage values
  if (!has_rename_filters && mqc_config.sample_names_rename && mqc_config.sample_names_rename.length > 0) {
    let mqc_renamesamples_idx = 300;

    // Add each pattern
    for (let i = 0; i < mqc_config.sample_names_rename.length; i++) {
      const pattern = mqc_config.sample_names_rename[i];
      if (!Array.isArray(pattern) || pattern.length < 2) continue;

      const from_text = pattern[0];
      const to_text = pattern[1];

      // Add to the filters list
      $("#mqc_renamesamples_filters").append(
        '<li><input class="f_text from_text" value="' +
          from_text +
          '" tabindex="' +
          mqc_renamesamples_idx +
          '" />' +
          '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="' +
          to_text +
          '" tabindex="' +
          (mqc_renamesamples_idx + 1) +
          '" />' +
          '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );

      mqc_renamesamples_idx += 2;
    }

    // Apply the renames
    apply_mqc_renamesamples();
  }
});

//////////////////////////////////////////////////////
// UTILITY FUNCTIONS
//////////////////////////////////////////////////////
function hashCode(str) {
  var hash = 0;
  if (str.length == 0) return hash;
  for (i = 0; i < str.length; i++) {
    char = str.charCodeAt(i);
    hash = (hash << 5) - hash + char;
    hash = hash & hash; // Convert to 32bit integer
  }
  return hash;
}

// Convert RGB color to hex format
function rgbToHex(rgb) {
  // Extract numbers from rgb(r, g, b) format
  const matches = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
  if (!matches) return rgb; // Return original if not RGB format

  // Convert each component to hex
  function componentToHex(c) {
    const hex = parseInt(c).toString(16);
    return hex.length === 1 ? "0" + hex : hex;
  }

  // Combine components with # prefix
  return "#" + componentToHex(matches[1]) + componentToHex(matches[2]) + componentToHex(matches[3]);
}

//////////////////////////////////////////////////////
// GENERAL TOOLBOX FUNCTIONS
//////////////////////////////////////////////////////
function mqc_toolbox_openclose(target, open, dataTab) {
  // Hide any open tooltip so it's not left dangling
  $(".mqc-toolbox-buttons li a").tooltip("hide");
  // Find if what we clicked is already open
  let btn = $('.mqc-toolbox-buttons li a[href="' + target + '"]');
  if (open === undefined) {
    open = !btn.hasClass("active");
  }
  let already_open = $(".mqc-toolbox").hasClass("active");
  if (open) {
    if (already_open) {
      mqc_toolbox_confirmapply();
    }
    $(".mqc-toolbox, .mqc-toolbox-buttons li a, .mqc_filter_section").removeClass("active");
    btn.addClass("active");
    $(".mqc-toolbox, " + target).addClass("active");
    $(document).trigger("mqc_toolbox_open");
    let timeout = already_open ? 0 : 510;
    setTimeout(function () {
      if (target === "#mqc_cols") {
        $("#mqc_colour_filter").focus();
      }
      if (target === "#mqc_renamesamples") {
        $("#mqc_renamesamples_from").focus();
      }
      if (target === "#mqc_hidesamples") {
        $("#mqc_hidesamples_filter").focus();
      }
    }, timeout);
    if (dataTab) {
      $('#mqc_exportplots a[href="#mqc_data_download"]').tab("show");
    } else {
      $('#mqc_exportplots a[href="#mqc_image_download"]').tab("show");
    }
  } else {
    mqc_toolbox_confirmapply();
    btn.removeClass("active");
    $(".mqc-toolbox, .mqc-toolbox-buttons li a").removeClass("active");
    $(document).trigger("mqc_toolbox_close");
  }
}
function mqc_toolbox_confirmapply() {
  // Check if there's anything waiting to be applied
  if ($("#mqc_cols_apply").is(":enabled") && $("#mqc_cols").is(":visible")) {
    $.toast({
      heading: "Highlights Not Applied",
      text: "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      icon: "warning",
      position: "bottom-right",
      hideAfter: 5000,
    });
  }
  if ($("#mqc_rename_apply").is(":enabled") && $("#mqc_renamesamples").is(":visible")) {
    $.toast({
      heading: "Rename Patterns Not Applied",
      text: "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      icon: "warning",
      position: "bottom-right",
      hideAfter: 5000,
    });
  }
  if ($("#mqc_hide_apply").is(":enabled") && $("#mqc_hidesamples").is(":visible")) {
    $.toast({
      heading: "Hide Samples Not Applied",
      text: "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      icon: "warning",
      position: "bottom-right",
      hideAfter: 5000,
    });
  }
}

function validate_regexp(pattern) {
  try {
    new RegExp(pattern, "g");
    return true;
  } catch (error) {
    $.toast({
      heading: "Invalid Regular Expression!",
      text:
        "Apologies, your regular expression pattern is invalid: <code>" +
        pattern +
        "</code><br><br>" +
        'For more help and testing, try it out at <a href="https://regex101.com/" target="_blank">regex101.com</a>.',
      icon: "error",
      position: "bottom-right",
      hideAfter: 5000,
    });
    return false;
  }
}

//////////////////////////////////////////////////////
// HIGHLIGHT SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_highlights() {
  // Collect the filters into an array
  var f_texts = [];
  var f_cols = [];
  var regex_mode = $("#mqc_cols .mqc_regex_mode .re_mode").hasClass("on");
  var num_errors = 0;

  $("#mqc_col_filters li").each(function () {
    var inputElement = $(this).find(".f_text");
    var pattern = inputElement.val();
    // Validate RegExp
    $(this).removeClass("bg-danger");
    if (regex_mode && !validate_regexp(pattern)) {
      $(this).addClass("bg-danger");
      num_errors++;
    }

    // Only add pattern if it hasn't already been added
    if (pattern.length > 0 && f_texts.indexOf(pattern) < 0) {
      f_texts.push(pattern);
      f_cols.push(inputElement.css("color"));
    } else {
      f_cols[f_texts.indexOf(pattern)] = inputElement.css("color");
    }
  });
  if (num_errors > 0) {
    return false;
  }

  // Apply a 'background' highlight to remove default colouring first
  // Also highlight toolbox drawer icon
  if (f_texts.length > 0 && f_texts.indexOf("") < 0) {
    f_texts.unshift("");
    f_cols.unshift(null);
    $('.mqc-toolbox-buttons a[href="#mqc_cols"]').addClass("in_use");
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_cols"]').removeClass("in_use");
  }

  window.mqc_highlight_f_texts = f_texts;
  window.mqc_highlight_f_cols = f_cols;
  window.mqc_highlight_regex_mode = regex_mode;

  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger("mqc_highlights", [f_texts, f_cols, regex_mode]);

  return true;
}

//////////////////////////////////////////////////////
// RENAME SAMPLES
//////////////////////////////////////////////////////

function apply_mqc_renamesamples() {
  let valid_from_texts = [];
  let valid_to_texts = [];
  let regex_mode = $("#mqc_renamesamples .mqc_regex_mode .re_mode").hasClass("on");
  let num_errors = 0;
  // Collect filters
  $("#mqc_renamesamples_filters > li").each(function () {
    let from_text = $(this).find(".from_text").val();
    // Validate RegExp
    $(this).removeClass("bg-danger");
    if (regex_mode) {
      if (!validate_regexp(from_text)) {
        $(this).addClass("bg-danger");
        num_errors++;
        return;
      }
      from_text = new RegExp(from_text, "g");
    }
    valid_from_texts.push(from_text);
    let to_text = $(this).find(".to_text").val();
    valid_to_texts.push(to_text);
  });
  if (num_errors > 0) {
    return false;
  }

  // If something was renamed, highlight the toolbox icon
  if (valid_from_texts.length > 0) {
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').addClass("in_use");
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_renamesamples"]').removeClass("in_use");
  }

  window.mqc_rename_f_texts = valid_from_texts;
  window.mqc_rename_t_texts = valid_to_texts;

  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger("mqc_renamesamples", [window.mqc_rename_f_texts, window.mqc_rename_t_texts]);

  return true;
}

//////////////////////////////////////////////////////
// HIDE SAMPLES
//////////////////////////////////////////////////////
function apply_mqc_hidesamples(mode) {
  // Collect the filters into an array
  if (mode === undefined) {
    mode = $(".mqc_hidesamples_showhide:checked").val() === "show" ? "show" : "hide";
  }
  let regex_mode = $("#mqc_hidesamples .mqc_regex_mode .re_mode").hasClass("on");
  let f_texts = [];
  let num_errors = 0;
  $("#mqc_hidesamples_filters li").each(function () {
    let pattern = $(this).find(".f_text").val();
    // Validate RegExp
    $(this).removeClass("bg-danger");
    if (regex_mode && !validate_regexp(pattern)) {
      $(this).addClass("bg-danger");
      num_errors++;
    }
    f_texts.push(pattern);
  });
  if (num_errors > 0) {
    return false;
  }

  // If something was hidden, highlight the toolbox icon
  if (f_texts.length > 0) {
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').addClass("in_use");
  } else {
    $('.mqc-toolbox-buttons a[href="#mqc_hidesamples"]').removeClass("in_use");
  }

  window.mqc_hide_mode = mode;
  window.mqc_hide_f_texts = f_texts;
  window.mqc_hide_regex_mode = regex_mode;

  // Fire off a custom jQuery event for other javascript chunks to tie into
  $(document).trigger("mqc_hidesamples", [f_texts, regex_mode]);

  return true;
}

//////////////////////////////////////////////////////
// SAVE TOOLBOX SETTINGS
//////////////////////////////////////////////////////

// Add these helper functions
function getConfigObject() {
  return {
    highlights_f_texts: window.mqc_highlight_f_texts,
    highlights_f_cols: window.mqc_highlight_f_cols,
    highlight_regex: window.mqc_highlight_regex_mode,
    rename_from_texts: window.mqc_rename_f_texts,
    rename_to_texts: window.mqc_rename_t_texts,
    rename_regex: window.mqc_rename_regex_mode,
    hidesamples_mode: window.mqc_hide_mode,
    hidesamples_f_texts: window.mqc_hide_f_texts,
    hidesamples_regex: window.mqc_hide_regex_mode,
  };
}

// Save the current configuration setup
function mqc_save_config(name, clear, as_default) {
  if (name === undefined) {
    return false;
  }
  const config = getConfigObject();
  var prev_config = {};
  // Load existing configs (inc. from other reports)
  try {
    try {
      prev_config = localStorage.getItem("mqc_config");
      if (prev_config !== null && prev_config !== undefined) {
        prev_config = JSON.parse(prev_config);
      } else {
        prev_config = {};
      }

      // Update config obj with current config
      if (clear == true) {
        delete prev_config[name];
      } else {
        prev_config[name] = config;
        prev_config[name]["last_updated"] = Date();
        if (as_default) {
          for (var c in prev_config) {
            if (prev_config.hasOwnProperty(c)) {
              prev_config[c]["default"] = false;
            }
          }
        }
        prev_config[name]["default"] = as_default;
        if (as_default) console.log("Set new default config!");
      }
      localStorage.setItem("mqc_config", JSON.stringify(prev_config));
    } catch (e) {
      console.log("Could not access localStorage");
    }

    if (clear === true) {
      // Remove from load select box
      $("#mqc_loadconfig_form select option:contains('" + name + "')").remove();
      // Successfully deleted message
      $('<p class="text-danger" id="mqc-cleared-success">Settings deleted.</p>')
        .hide()
        .insertBefore($("#mqc_loadconfig_form .actions"))
        .slideDown(function () {
          setTimeout(function () {
            $("#mqc-cleared-success").slideUp(function () {
              $(this).remove();
            });
          }, 5000);
        });
    } else {
      // Remove from load select box
      $("#mqc_loadconfig_form select option:contains('" + name + "')").remove();
      // Add new name to load select box and select it
      $("#mqc_loadconfig_form select")
        .prepend("<option>" + name + (as_default ? " [default]" : "") + "</option>")
        .val(name + (as_default ? " [default]" : ""));
      // Success message
      $('<p class="text-success" id="mqc-save-success">Settings saved.</p>')
        .hide()
        .insertBefore($("#mqc_saveconfig_form"))
        .slideDown(function () {
          setTimeout(function () {
            $("#mqc-save-success").slideUp(function () {
              $(this).remove();
            });
          }, 5000);
        });
    }
  } catch (e) {
    console.log("Error updating localstorage: " + e);
  }
}

// Clear current default configuration
function mqc_clear_default_config() {
  try {
    var config = localStorage.getItem("mqc_config");
    if (!config) {
      return;
    } else {
      config = JSON.parse(config);
    }
    for (var c in config) {
      if (config.hasOwnProperty(c)) {
        config[c]["default"] = false;
      }
    }
    localStorage.setItem("mqc_config", JSON.stringify(config));
    $('<p class="text-danger" id="mqc-cleared-success">Unset default.</p>')
      .hide()
      .insertBefore($("#mqc_loadconfig_form .actions"))
      .slideDown(function () {
        setTimeout(function () {
          $("#mqc-cleared-success").slideUp(function () {
            $(this).remove();
          });
        }, 5000);
        var name = $('#mqc_loadconfig_form select option:contains("default")').text();
        $('#mqc_loadconfig_form select option:contains("default")').remove();
        name = name.replace(" [default]", "");
        $("#mqc_loadconfig_form select")
          .append("<option>" + name + "</option>")
          .val(name);
      });
  } catch (e) {
    console.log("Could not access localStorage");
  }
}

//////////////////////////////////////////////////////
// LOAD TOOLBOX SAVE NAMES
//////////////////////////////////////////////////////
function populate_mqc_saveselect() {
  var default_config = "";
  try {
    var local_config = localStorage.getItem("mqc_config");
    if (local_config !== null && local_config !== undefined) {
      local_config = JSON.parse(local_config);
      default_name = false;
      for (var name in local_config) {
        if (local_config[name]["default"]) {
          console.log("Loaded default config!");
          load_mqc_config(name);
          default_config = name;
          name = name + " [default]";
          default_name = name;
        }
        $("#mqc_loadconfig_form select")
          .append("<option>" + name + "</option>")
          .val(name);
      }
      // Set the selected select option
      if (default_name !== false) {
        $('#mqc_loadconfig_form select option:contains("' + default_name + '")').prop("selected", true);
      } else {
        $("#mqc_loadconfig_form select option:first").prop("selected", true);
      }
    }
  } catch (e) {
    console.log("Could not load local config: " + e);
    $("#mqc_saveconfig").html(
      "<h4>Error accessing localStorage</h4>" +
        '<p>This feature uses a web browser feature called "localStorage". ' +
        "We're not able to access this at the moment, which probably means that " +
        'you have the <em>"Block third-party cookies and site data"</em> setting ticked (Chrome) ' +
        "or equivalent in other browsers.</p><p>Please " +
        '<a href="https://www.google.se/search?q=Block+third-party+cookies+and+site+data" target="_blank">change this browser setting</a>' +
        " to save MultiQC report configs.</p>",
    );
  }
}

//////////////////////////////////////////////////////
// LOAD TOOLBOX SETTINGS
//////////////////////////////////////////////////////
function load_mqc_config(name, config_obj) {
  if (name === undefined) {
    return false;
  }
  var config = {};
  if (config_obj) {
    // Use provided config object
    config = config_obj;
  } else {
    // Load from localStorage
    try {
      var local_config = localStorage.getItem("mqc_config");
      if (local_config !== null && local_config !== undefined) {
        local_config = JSON.parse(local_config);
        config = local_config[name];
      }
    } catch (e) {
      console.log("Could not access localStorage");
    }
  }

  // Apply config - rename samples
  if (notEmptyObj(config["rename_regex"])) {
    if (config["rename_regex"] === true) {
      $("#mqc_renamesamples .mqc_regex_mode .re_mode").removeClass("off").addClass("on").text("on");
      window.mqc_rename_regex_mode = true;
    }
  }
  if (notEmptyObj(config["rename_from_texts"]) && notEmptyObj(config["rename_to_texts"])) {
    window.mqc_rename_f_texts = [];
    window.mqc_rename_t_texts = [];
    $.each(config["rename_from_texts"], function (idx, from_text) {
      var to_text = config["rename_to_texts"][idx];
      if (from_text.length === 0) {
        return true;
      }
      var li = '<li><input class="f_text from_text" value="' + from_text + '" />';
      li +=
        '<small class="glyphicon glyphicon-chevron-right"></small><input class="f_text to_text" value="' +
        to_text +
        '" />';
      li +=
        '<button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>';
      window.mqc_rename_f_texts.push(from_text);
      window.mqc_rename_t_texts.push(to_text);
      $("#mqc_renamesamples_filters").append(li);
    });
    $(document).trigger("mqc_renamesamples", [
      window.mqc_rename_f_texts,
      window.mqc_rename_t_texts,
      config["rename_regex"],
    ]);
  }

  // Apply config - highlights
  if (notEmptyObj(config["highlight_regex"])) {
    if (config["highlight_regex"] === true) {
      $("#mqc_cols .mqc_regex_mode .re_mode").removeClass("off").addClass("on").text("on");
      window.mqc_highlight_regex_mode = true;
    }
  }
  if (notEmptyObj(config["highlights_f_texts"]) && notEmptyObj(config["highlights_f_cols"])) {
    window.mqc_highlight_f_texts = [];
    window.mqc_highlight_f_cols = [];
    $.each(config["highlights_f_texts"], function (idx, f_text) {
      var f_col = config["highlights_f_cols"][idx];
      $("#" + hashCode(f_text + f_col)).remove();
      $("#mqc_col_filters").append(
        '<li style="color:' +
          f_col +
          ';" id="' +
          hashCode(f_text + f_col) +
          '"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="' +
          f_text +
          '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );
      window.mqc_highlight_f_texts.push(f_text);
      window.mqc_highlight_f_cols.push(f_col);
      mqc_colours_idx += 1;
    });
    $("#mqc_colour_filter_color").val(mqc_colours[mqc_colours_idx]);
    $(document).trigger("mqc_highlights", [
      window.mqc_highlight_f_texts,
      window.mqc_highlight_f_cols,
      config["highlight_regex"],
    ]);
  }

  // Apply config - hide samples
  if (notEmptyObj(config["hidesamples_regex"])) {
    if (config["hidesamples_regex"] == true) {
      $("#mqc_hidesamples .mqc_regex_mode .re_mode").removeClass("off").addClass("on").text("on");
      window.mqc_hide_regex_mode = true;
    }
  }
  if (notEmptyObj(config["hidesamples_mode"])) {
    if (config["hidesamples_mode"] == "show") {
      $(".mqc_hidesamples_showhide").prop("checked", false);
      $(".mqc_hidesamples_showhide[val=show]").prop("checked", true);
      window.mqc_hide_mode = "show";
    }
  }
  if (notEmptyObj(config["hidesamples_f_texts"])) {
    window.mqc_hide_f_texts = [];
    $.each(config["hidesamples_f_texts"], function (idx, f_text) {
      if (f_text.length == 0) {
        return true;
      }
      $("#mqc_hidesamples_filters").append(
        '<li><input class="f_text" value="' +
          f_text +
          '" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>',
      );
      window.mqc_hide_f_texts.push(f_text);
    });
    $(document).trigger("mqc_hidesamples", [window.mqc_hide_f_texts, config["hidesamples_regex"]]);
  }

  // Trigger loaded event to initialise plots
  $(document).trigger("mqc_config_loaded");
}

//////////////////////////////////////////////////////
// SAVE SETTINGS TO FILE
//////////////////////////////////////////////////////
function downloadConfigFile(name) {
  const config = getConfigObject();
  const blob = new Blob([JSON.stringify(config, null, 2)], { type: "application/json" });
  const filename = `multiqc_config_${name || "settings"}.json`;
  saveAs(blob, filename);
}

function loadConfigFromFile(file) {
  const reader = new FileReader();
  reader.onload = function (e) {
    try {
      const config = JSON.parse(e.target.result);

      // Validate that this config matches the current report
      if (config.report_id && config.report_id !== window.reportUuid) {
        if (!confirm("This configuration is from a different report. Load it anyway?")) {
          return;
        }
      }

      // Load the config
      load_mqc_config("custom_file_load", config);

      // Auto-save the newly loaded config
      mqc_auto_save_config();

      // Show success message
      $.toast({
        heading: "Configuration Loaded",
        text: "Settings have been loaded from file successfully",
        icon: "success",
        position: "bottom-right",
      });
    } catch (err) {
      $.toast({
        heading: "Error Loading Configuration",
        text: "Could not parse the configuration file: " + err.message,
        icon: "error",
        position: "bottom-right",
      });
    }
  };
  reader.readAsText(file);
}

$(function () {
  // // Add file input handler
  // $("#mqc_load_config_file").on("change", function (e) {
  //   if (this.files && this.files[0]) {
  //     loadConfigFromFile(this.files[0]);
  //   }
  // });
  // // Add download handler
  // $("#mqc_download_config").click(function (e) {
  //   e.preventDefault();
  //   const name = $("#mqc_saveconfig_form input").val().trim() || "settings";
  //   downloadConfigFile(name);
  // });

  // Lazy load file input handler
  let fileInputInitialized = false;
  $("#mqc_saveconfig").on("mouseenter", function () {
    if (!fileInputInitialized) {
      // Initialize file input handler only when user interacts with the save section
      $("#mqc_load_config_file_wrapper").append('<input type="file" id="mqc_load_config_file" accept=".json">');
      $("#mqc_load_config_file").on("change", function (e) {
        if (this.files && this.files[0]) {
          loadConfigFromFile(this.files[0]);
        }
      });
      fileInputInitialized = true;
    }
  });
  // Add download handler
  $("#mqc_download_config").click(function (e) {
    e.preventDefault();
    const name = $("#mqc_saveconfig_form input").val().trim() || "settings";
    downloadConfigFile(name);
  });
});

//////////////////////////////////////////////////////
// SAVE SETTINGS AUTOMATICALLY
//////////////////////////////////////////////////////
const AUTO_SAVE_PREFIX = "autosave_";

function mqc_auto_save_config() {
  // Get the current config
  var config = getConfigObject();
  config.last_updated = Date();

  try {
    // Save to localStorage with report UUID
    var prev_config = {};
    try {
      prev_config = localStorage.getItem("mqc_config");
      prev_config = prev_config ? JSON.parse(prev_config) : {};
    } catch (e) {
      console.log("Could not access localStorage");
    }

    // Add/update the auto-saved config
    const autoSaveKey = AUTO_SAVE_PREFIX + window.reportUuid;
    prev_config[autoSaveKey] = config;
    localStorage.setItem("mqc_config", JSON.stringify(prev_config));
  } catch (e) {
    console.log("Error auto-saving config:", e);
  }
}

// Storage helper functions
function saveToLocalStorage(key, value) {
  try {
    localStorage.setItem(key, value);
  } catch (e) {
    console.warn("Failed to save to localStorage:", e);
  }
}

function getFromLocalStorage(key) {
  try {
    return localStorage.getItem(key);
  } catch (e) {
    console.warn("Failed to read from localStorage:", e);
    return null;
  }
}

function updatePanel(providerId) {
  const provider = AI_PROVIDERS[providerId];

  // Add label dynamically
  let aiModelInfo = "";
  let aiApiKeyInfo = "";

  // Update model field with stored or default value for new provider
  if (providerId === "seqera") {
    $(".ai-generate-button-wrapper").show();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_provider_info").html("");
    $("#ai_model_group").hide();
    $("#ai_endpoint_group").hide();
    $("#ai_query_options_group").hide();
    $("#ai_api_key_group").show();
    $("#ai_api_key_info").html(
      `<span id="ai_api_key_info_required">Field is required.</span> You can find your API key in the <a style='text-decoration: underline;' href='${provider.apiKeysUrl}' target='_blank'>${provider.name} dashboard</a>`,
    );
    $("#ai_provider_logo").html(
      '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 386.21 80" height="22" width="112" class="_logo_1pric_80"><defs><clipPath id="Logo_svg__a"><path d="M0 0h386.21v80H0z" style="fill:none;stroke-width:0"></path></clipPath><clipPath id="Logo_svg__b"><path d="M0 0h79.68v80H0z" style="fill:none;stroke-width:0"></path></clipPath></defs><g style="clip-path:url(#Logo_svg__a)"><g style="clip-path:url(#Logo_svg__b)"><path d="m67.65 41.42.11.04.11.04.11.05.11.06.1.06.1.06.1.07.09.07.09.08.09.08.08.09.08.09.08.1.07.1.07.1.06.1.05.11.05.11.04.11.04.11.03.12.03.11.02.12.02.11v.12l.02.12v.36l-.02.12-.02.12-.02.12-.03.12-.03.11-.04.12-.04.11-.05.11-.05.11-.06.1-.07.1-.07.1-.08.09-.08.09-.08.08-.09.09-.09.08-.09.07-.1.07-.1.07-.1.06-.11.06-.11.05-.11.05-.11.04-.11.04-.11.03-.11.03-.12.02h-.11l-.12.03h-.02.12l.22-.01h.22l.22-.03.21-.03.22-.03.21-.04.21-.04.21-.05.21-.05.21-.06.21-.07.2-.07.2-.08.2-.08.2-.09.2-.09.19-.1.19-.1.19-.11.18-.12.18-.12.18-.12.18-.13.17-.13.17-.13.17-.14.16-.15.16-.15.16-.16.15-.16.14-.16.14-.17.14-.17.13-.18.13-.18.12-.18.12-.19.11-.19.11-.19.11-.19.1-.19.09-.2.09-.2.08-.2.08-.2.07-.21.06-.21.06-.21.05-.21.05-.21.04-.22.04-.21.03-.21.02-.22.02-.22.02-.22v-.88l-.02-.22-.02-.22-.02-.21-.03-.22-.04-.22-.04-.21-.05-.21-.04-.16H48.51l19.02 5.25.13.04Z" style="fill:#160f26;stroke-width:0"></path><path d="m79.67 43.32-.02-.28-.02-.27-.03-.28-.03-.27-.04-.27-.04-.27-.05-.27-.05-.27-.06-.27-.07-.27-.07-.27-.08-.26-.08-.26-.09-.26-.1-.26-.1-.25-.11-.26-.11-.25-.11-.25-.12-.25-.13-.24-.13-.24-.14-.24-.14-.24-.14-.23-.15-.23-.15-.22-.16-.22-.05-.07h-2.38l.04.16.05.21.04.21.04.22.03.22.02.21.02.22.02.22v.88l-.02.22-.02.22-.02.22-.03.21-.04.21-.04.22-.05.21-.05.21-.06.21-.06.21-.07.21-.08.2-.08.2-.09.2-.09.2-.1.19-.11.19-.11.19-.11.19-.12.19-.12.18-.13.18-.13.18-.14.17-.14.17-.14.16-.15.16-.16.16-.16.15-.16.15-.17.14-.17.13-.17.13-.18.13-.18.12-.18.12-.18.12-.19.11-.19.1-.19.1-.2.09-.2.09-.2.08-.2.08-.2.07-.21.07-.21.06-.21.05-.21.05-.21.04-.21.04-.22.03-.21.03-.22.02h-.22l-.22.02H0v10.32h67.38l.27-.02.27-.02.27-.03.27-.03.27-.04.27-.04.27-.05.27-.06.27-.06.26-.06.26-.07.26-.08.26-.08.26-.09.26-.09.26-.1.25-.1.25-.11.25-.12.25-.12.24-.12.24-.13.24-.13.24-.14.23-.14.23-.15.23-.16.22-.16.22-.17.21-.17.21-.18.21-.18.2-.18.2-.19.19-.19.19-.2.18-.2.18-.2.18-.21.17-.21.17-.21.17-.22.16-.22.15-.22.15-.23.14-.23.14-.24.14-.24.13-.24.13-.24.12-.25.11-.25.11-.25.11-.25.1-.25.1-.26.09-.26.08-.26.08-.27.07-.26.07-.27.06-.27.05-.27.05-.27.04-.27.04-.27.03-.27.03-.27.02-.28.02-.27v-.27l.01-.28v-.28l-.01-.27Z" style="fill:#fa6863;stroke-width:0"></path><path d="m12.04 38.57-.11-.04-.11-.04-.11-.05-.11-.06-.1-.06-.1-.06-.1-.07-.09-.07-.09-.08-.09-.08-.08-.09-.08-.09-.08-.1-.07-.1-.07-.1-.06-.1-.05-.11-.05-.11-.04-.11-.04-.11-.03-.12-.03-.11-.03-.12-.02-.11v-.12l-.02-.12v-.36l.02-.12.02-.12.03-.12.03-.12.03-.11.04-.12.04-.11.05-.11.05-.1.06-.11.07-.1.07-.1.08-.09.08-.09.08-.09.09-.08.09-.08.09-.07.1-.07.1-.07.1-.06.11-.05.11-.05.11-.05.11-.04.11-.04.11-.03.11-.03.12-.02h.11l.12-.03h.02-.34l-.22.03-.22.02-.21.03-.22.03-.22.04-.21.04-.21.05-.21.05-.21.06-.21.07-.2.07-.21.07-.2.08-.2.09-.2.09-.19.1-.19.1-.18.11-.18.12-.18.12-.18.12-.18.13-.18.13-.17.14-.17.14-.16.15-.16.15-.16.15-.15.16-.14.16-.14.17-.13.17-.14.17-.12.18-.12.18-.12.18-.11.19-.11.19-.11.2-.1.19-.09.2-.09.2-.08.2-.08.2-.07.21-.06.21-.06.21-.05.21-.05.21-.04.21-.04.22-.03.21-.02.22-.02.21-.02.22v.88l.02.22.02.22.02.22.03.22.04.21.04.21.05.21.04.17h26.23l-19.02-5.25-.13-.04Z" style="fill:#160f26;stroke-width:0"></path><path d="M.01 36.67v.28l.04.27.03.27.03.27.04.27.04.27.05.27.05.27.06.27.07.27.07.27.08.26.08.26.09.26.1.26.1.25.11.26.11.25.11.25.12.24.13.24.13.24.14.24.14.23.14.23.15.23.16.23.16.22.05.07h2.38l-.04-.17-.05-.21-.04-.21-.04-.21-.03-.22-.02-.22-.02-.22-.02-.22v-.88l.02-.22.02-.21.02-.22.03-.21.04-.22.04-.21.05-.21.05-.21.06-.21.06-.21.07-.21.08-.2.08-.2.09-.2.09-.2.1-.19.11-.2.11-.19.11-.19.12-.18.12-.18.12-.18.14-.17.13-.17.14-.17.14-.16.15-.16.16-.15.16-.15.16-.15.17-.14.17-.14.18-.13.18-.13.18-.12.18-.12.18-.12.18-.11.19-.1.19-.1.2-.09.2-.09.2-.08.21-.07.2-.07.21-.07.21-.06.21-.05.21-.05.21-.04.22-.04.22-.03.21-.03.22-.02.22-.02h.22l.22-.01H79.7V23.23H12.31l-.27.02-.27.02-.27.03-.27.03-.27.04-.27.04-.27.05-.27.05-.26.06-.26.07-.26.07-.26.08-.26.08-.26.09-.26.09-.26.1-.25.1-.25.11-.25.12-.25.12-.24.12-.24.13-.24.13-.24.14-.23.15-.23.15-.22.16-.22.16-.22.16-.21.17-.21.18-.21.18-.2.18-.2.19-.19.19-.19.2-.19.2-.18.2-.18.21-.17.21-.17.21-.17.22-.16.22-.16.22-.15.23-.14.23-.14.24-.14.24-.13.24-.13.24-.12.25-.11.25-.11.25-.11.25-.1.26-.1.26-.09.26-.08.26-.08.26-.07.27-.07.27-.06.27-.05.27-.05.27-.04.27-.04.27-.03.27-.03.27-.02.27v.27l-.03.28v.82Z" style="fill:#f18046;stroke-width:0"></path><path d="m12.04 15.35-.11-.04-.11-.05-.11-.05-.1-.06-.11-.06-.1-.06-.1-.07-.1-.07-.09-.08-.09-.08-.08-.09-.08-.09-.08-.1-.07-.1-.06-.1-.06-.1-.05-.11-.05-.11-.04-.11-.04-.12-.04-.11-.03-.11-.03-.12-.02-.12-.02-.12v-.48l.02-.12.02-.12.03-.12.03-.12.04-.12.04-.11.04-.11.05-.11.05-.11.06-.1.06-.1.07-.1.08-.09.08-.09.08-.09.09-.08.09-.08.1-.07.1-.07.1-.07.11-.06.1-.06.11-.05.11-.04.11-.04.11-.04.11-.03.12-.02.11-.02h.12l.12-.02h.04-.36l-.22.03-.22.02-.21.03-.22.03-.22.04-.21.04-.21.05-.21.05-.21.06-.21.07-.2.07-.21.08-.2.08-.2.09-.2.09-.19.1-.19.1-.18.11-.19.12-.18.12-.18.12-.18.13-.18.13-.17.14-.16.14-.16.14-.16.15-.15.15-.15.16-.14.16-.14.17-.13.17-.14.18-.12.18-.12.18-.12.18-.11.19-.11.19-.11.2-.1.19-.09.2-.09.2-.08.2-.08.2-.07.21-.06.21-.06.21-.05.21-.05.21-.04.21-.04.22-.03.21-.03.22-.02.21-.02.22v.88l.02.22.02.22.03.21.03.22.04.22.04.21.05.21.05.21.06.21.06.21.07.2.08.2.08.2.09.2.09.2.1.2.06.11.27-.14.27-.13.27-.13.27-.12.27-.12.28-.11.28-.1.28-.11.28-.09.29-.09.29-.08.29-.08.29-.07.29-.06.29-.06.29-.06.3-.05.3-.04.3-.04.3-.03.3-.03.3-.02h.29l.31-.02h18.63l-19.02-5.25-.13-.04ZM66.83 69.68h.44l.22-.03.22-.02.22-.03.22-.03.21-.04.21-.04.21-.05.21-.05.21-.06.21-.07.21-.07.2-.08.2-.08.2-.09.2-.09.19-.1.19-.1.19-.11.18-.12.18-.12.18-.12.18-.13.17-.13.17-.13.17-.14.16-.15.16-.15.15-.16.15-.16.15-.16.14-.17.14-.17.13-.18.12-.18.12-.18.12-.19.11-.19.11-.19.11-.19.1-.2.09-.2.08-.2.09-.2.07-.2.07-.21.07-.21.06-.21.06-.21.05-.21.04-.22.04-.21.03-.21.03-.22.02-.21v-.22l.02-.22v-.66l-.02-.22-.02-.22-.03-.21-.03-.22-.04-.22-.04-.21-.05-.21-.06-.21-.06-.21-.07-.21-.07-.21-.07-.2-.09-.2-.08-.2-.09-.2-.1-.2-.06-.12-.27.14-.27.13-.27.13-.27.12-.28.12-.28.11-.28.1-.28.1-.29.1-.28.09-.28.08-.29.08-.29.07-.29.07-.29.06-.29.05-.29.05-.3.04-.3.04-.3.03-.3.03-.3.02-.3.02h-.3l-.3.01H48.51l19.02 5.25.13.04.11.04.11.04.11.05.11.06.1.06.1.06.1.07.09.07.09.08.09.08.08.09.08.09.08.1.07.1.07.1.06.1.05.11.05.11.04.11.04.11.03.12.03.11.02.12.02.12v.12l.02.12v.36l-.02.12-.02.12-.02.12-.03.12-.03.11-.04.12-.04.11-.05.11-.05.1-.06.11-.07.1-.07.1-.08.09-.08.09-.08.09-.09.08-.09.08-.09.07-.1.07-.1.07-.1.06-.11.05-.11.05-.11.05-.11.04-.11.04-.11.03-.11.03-.12.02h-.11l-.12.03h-.12" style="fill:#160f26;stroke-width:0"></path><path d="M79.67 66.55v-.27l-.03-.27-.03-.28-.03-.27-.04-.27-.04-.27-.05-.27-.05-.27-.06-.27-.07-.26-.07-.27-.08-.27-.08-.26-.09-.26-.1-.26-.1-.25-.11-.26-.11-.25-.11-.25-.12-.25-.13-.24-.13-.24-.14-.24-.14-.24-.14-.23-.15-.23-.16-.22-.16-.22-.17-.22-.17-.21-.18-.21-.18-.2-.18-.2-.19-.2-.19-.2-.2-.19-.2-.19-.2-.18-.21-.18-.21-.18h-.01l-.3.17-.26.15-.27.15-.04.02.07.13.1.2.09.2.08.2.09.2.07.2.07.21.07.21.06.21.06.21.05.21.04.21.04.22.03.22.03.21.02.22v.22l.02.22v.65l-.02.22-.02.21-.03.22-.03.21-.04.21-.04.22-.05.21-.06.21-.06.21-.07.21-.07.21-.07.2-.09.2-.08.2-.09.2-.1.2-.11.19-.11.19-.11.19-.12.19-.12.18-.12.18-.13.18-.14.17-.14.17-.15.16-.15.16-.15.16-.16.15-.16.15-.17.14-.17.13-.17.13-.18.13-.18.12-.18.12-.18.12-.19.11-.19.1-.19.1-.2.09-.2.09-.2.08-.2.08-.21.07-.21.07-.21.06-.21.05-.21.05-.21.04-.21.04-.22.03-.22.03-.22.02h-.22l-.22.03H0V80h67.38l.27-.02.27-.02.27-.02.27-.03.27-.04.27-.04.27-.05.27-.06.26-.06.26-.06.26-.07.26-.08.26-.08.26-.09.26-.09.25-.1.25-.1.25-.11.25-.12.24-.12.25-.12.24-.13.24-.13.23-.14.23-.15.23-.15.22-.16.22-.16.22-.16.21-.17.21-.18.21-.18.2-.18.2-.19.2-.19.19-.19.19-.2.18-.2.18-.21.18-.21.17-.21.17-.22.16-.22.16-.22.15-.23.14-.23.14-.24.14-.24.13-.24.13-.24.12-.25.11-.25.11-.25.11-.25.1-.25.1-.26.09-.26.08-.26.08-.26.07-.26.07-.27.06-.27.05-.27.05-.27.04-.27.04-.27.03-.27.03-.27.02-.28v-.27l.03-.27v-.83Z" style="fill:#3d95fd;stroke-width:0"></path><path d="M.01 13.45v.28l.03.27.03.28.03.27.04.27.04.27.05.27.05.27.06.27.07.27.07.26.08.27.09.26.09.26.1.26.1.25.11.26.11.25.11.25.12.24.13.24.13.24.13.24.14.24.14.23.15.23.16.23.16.22.17.22.17.21.17.21.18.2.18.2.19.2.19.2.19.19.2.19.2.18.21.18.21.17.3-.17.26-.15.27-.15.04-.02-.07-.13-.1-.2-.09-.2-.09-.2-.08-.2-.08-.2-.07-.2-.06-.21-.06-.21-.05-.21-.05-.21-.04-.21-.04-.22-.03-.22-.03-.21-.02-.22-.02-.22v-.88l.02-.22.02-.21.03-.22.03-.21.04-.22.04-.21.05-.21.05-.21.06-.21.06-.21.07-.21.08-.2.08-.2.09-.2.09-.2.1-.19.11-.2.11-.19.11-.19.12-.18.12-.18.12-.18.14-.18.13-.17.14-.17.14-.16.15-.16.15-.15.16-.15.16-.14.16-.14.17-.14.18-.13.18-.13.18-.12.18-.12.19-.12.18-.11.19-.1.19-.1.2-.09.2-.09.2-.08.21-.08.2-.07.21-.07.21-.06.21-.05.21-.05.21-.04.22-.04.22-.03.21-.03.22-.02.22-.02h67.27V0H12.31l-.27.03-.27.02-.27.02-.27.03-.27.04-.27.04-.27.05-.27.05-.26.06-.26.06-.26.07-.26.08-.26.08-.26.09-.26.09-.26.1-.25.1-.25.11-.25.12-.25.12-.24.13-.24.13-.24.13-.24.14-.23.15-.23.15-.22.16-.22.16-.22.17-.21.17-.21.18-.21.18-.2.18-.2.19-.19.19-.19.2-.19.2-.18.2-.18.21-.17.21-.17.21-.17.22-.16.22-.16.22-.15.23-.14.23-.14.23-.13.24-.13.24-.13.24-.12.25-.11.25-.11.25-.11.25-.1.25-.1.26-.09.26-.09.26-.08.26-.07.27-.07.27-.06.27-.05.27-.05.27-.04.27-.04.27-.03.27-.03.27-.02.28v.27l-.03.27v.82Z" style="fill:#0dc09d;stroke-width:0"></path></g><path d="M115.57 65.34H95.13v-8.07h20.44c3.17 0 5.6-.61 7.32-1.85 1.71-1.23 2.57-2.82 2.57-4.77s-.81-3.37-2.42-4.28c-1.62-.91-4.01-1.65-7.17-2.24l-3.3-.58c-3.23-.58-6.17-1.43-8.82-2.53s-4.75-2.63-6.3-4.57c-1.55-1.95-2.33-4.44-2.33-7.49 0-4.54 1.68-8.06 5.04-10.55 3.36-2.5 7.82-3.74 13.37-3.74h19.96v8.17h-19.96c-2.71 0-4.84.5-6.39 1.51-1.55 1-2.33 2.42-2.33 4.23 0 1.95.76 3.37 2.28 4.28q2.28 1.365 6.15 2.04l3.39.58c3.42.58 6.56 1.4 9.4 2.43 2.84 1.04 5.09 2.53 6.74 4.47 1.65 1.95 2.47 4.54 2.47 7.78 0 4.8-1.78 8.53-5.33 11.19q-5.325 3.99-14.34 3.99M238.44 79.99v-22.5h-1.55c-.71 1.3-1.74 2.55-3.1 3.74-1.36 1.2-3.09 2.19-5.18 2.97q-3.15 1.17-7.8 1.17c-4.01 0-7.69-.97-11.05-2.92s-6.04-4.75-8.04-8.42c-2-3.66-3-8.09-3-13.28v-1.46c0-5.19 1.02-9.61 3.05-13.28 2.03-3.66 4.73-6.47 8.09-8.42s7.01-2.92 10.95-2.92c4.65 0 8.22.84 10.71 2.53s4.34 3.6 5.57 5.74h1.55v-8.27h9.79v65.3h-9.98Zm-14.82-23.38c4.39 0 7.98-1.4 10.75-4.18s4.17-6.78 4.17-11.96v-.88c0-5.12-1.41-9.08-4.21-11.87-2.81-2.79-6.38-4.18-10.71-4.18s-7.8 1.4-10.61 4.18c-2.81 2.79-4.21 6.75-4.21 11.87v.88c0 5.19 1.4 9.18 4.21 11.96 2.81 2.79 6.35 4.18 10.61 4.18M190.43 39.01c0-4.86-.97-9.11-2.91-12.74s-4.64-6.47-8.09-8.51c-3.46-2.04-7.48-3.07-12.06-3.07s-8.85 1.04-12.4 3.11a22.3 22.3 0 0 0-4.76 3.72c-.18.18-.35.36-.52.55-.51.56-.98 1.15-1.44 1.77-.6.82-1.16 1.7-1.67 2.62-2.03 3.7-3.05 8.04-3.05 13.04v1.17c0 5 1.02 9.34 3.05 13.04s4.86 6.57 8.48 8.61c3.61 2.04 7.85 3.07 12.69 3.07 4.26 0 7.82-.68 10.66-2.04s5.1-3.05 6.78-5.06q.705-.84 1.32-1.65c1.07-1.41 1.96-2.77 2.65-4.09l-8.33-4.28c-1.03 2.21-2.51 4.15-4.41 5.84-1.91 1.69-4.73 2.53-8.48 2.53-4.01 0-7.35-1.25-10.03-3.74-.17-.16-.33-.31-.48-.48-.16-.17-.31-.33-.46-.5-.36-.41-.69-.84-.99-1.3-1.3-1.96-2.05-4.3-2.24-7.01-.02-.26-.03-.52-.04-.78h36.72v-3.8Zm-30.04-3.89-5.1 5.12h-1.68c-.03-1.48 0-3.37.19-5.12.1-.67.22-1.31.39-1.93.16-.62.35-1.21.58-1.77.75-1.88 1.87-3.47 3.34-4.77 2.36-2.08 5.41-3.11 9.16-3.11s6.78 1.04 9.11 3.11c2.33 2.08 3.65 4.9 3.97 8.46h-19.96ZM303.6 39.01c0-4.86-.97-9.11-2.91-12.74s-4.64-6.47-8.09-8.51c-3.46-2.04-7.48-3.07-12.06-3.07s-8.85 1.04-12.4 3.11a22.3 22.3 0 0 0-4.76 3.72c-.17.18-.35.36-.52.55q-.765.84-1.44 1.77c-.6.82-1.16 1.7-1.67 2.62-2.03 3.7-3.05 8.04-3.05 13.04v1.17c0 5 1.02 9.34 3.05 13.04s4.86 6.57 8.48 8.61c3.61 2.04 7.85 3.07 12.69 3.07 4.26 0 7.81-.68 10.66-2.04 2.84-1.36 5.1-3.05 6.78-5.06q.705-.84 1.32-1.65c1.07-1.41 1.96-2.77 2.66-4.09l-8.33-4.28c-1.03 2.21-2.51 4.15-4.41 5.84-1.91 1.69-4.73 2.53-8.48 2.53-4.01 0-7.35-1.25-10.03-3.74-.17-.16-.33-.31-.48-.48-.16-.17-.31-.33-.46-.5-.36-.41-.69-.84-.99-1.3-1.3-1.96-2.05-4.3-2.24-7.01-.02-.26-.03-.52-.04-.78h36.73v-3.8Zm-30.04-3.89-5.1 5.12h-1.68c-.03-1.48 0-3.37.19-5.12.1-.67.22-1.31.39-1.93.16-.62.35-1.21.58-1.77.75-1.88 1.87-3.47 3.34-4.77 2.36-2.08 5.41-3.11 9.16-3.11s6.78 1.04 9.11 3.11c2.33 2.08 3.65 4.9 3.97 8.46h-19.96ZM312.21 65.37V14.69H322v5.64h1.55c.77-2.01 2.02-3.48 3.73-4.43s3.83-1.41 6.35-1.41h5.72v9.05h-6.1c-3.23 0-5.88.89-7.95 2.68s-3.1 4.52-3.1 8.22v30.93h-9.98ZM386.21 32.86c0-5.86-1.78-10.39-5.35-13.58s-8.47-4.79-14.7-4.79c-4.09 0-7.54.65-10.36 1.96s-5.09 3.03-6.81 5.18-2.97 4.5-3.75 7.03l9.34 3.03c.58-2.6 1.78-4.72 3.6-6.35s4.44-2.44 7.88-2.44 6.1.86 7.79 2.59c1.69 1.72 2.53 3.96 2.53 6.69v8.33h-1.62l-5.19-5.21h-7.79c-5.26 0-9.62 1.25-13.09 3.76s-5.21 6.2-5.21 11.09c0 3.26.8 6.02 2.38 8.3 1.59 2.28 3.73 4.01 6.42 5.18 2.7 1.17 5.76 1.76 9.19 1.76s5.97-.47 7.98-1.42c2.01-.94 3.54-2.05 4.57-3.32 1.04-1.27 1.81-2.43 2.34-3.47h1.56v8.21h8.27V32.85Zm-9.83 11.24c0 4.04-1.25 7.22-3.75 9.53-2.49 2.31-5.72 3.47-9.68 3.47-2.92 0-5.22-.65-6.91-1.96s-2.53-3.06-2.53-5.28.81-3.89 2.43-5.03c1.63-1.14 3.76-1.71 6.42-1.71h14.01v.98Z" style="fill:#160f26;stroke-width:0"></path></g></svg>',
    );
    $("#mqc_anonymize_samples_switch").show();
    $("#ai_context_window_group").hide();
  } else if (providerId === "clipboard") {
    $(".ai-generate-button-wrapper").hide();
    $(".ai-copy-button-wrapper").show();
    $("#ai_model_group").hide();
    $("#ai_api_key_group").hide();
    $("#ai_provider_info").html(
      `Copy AI LLM prompt, including report data, to the clipboard. For use with any 3rd party AI tools.`,
    );
    $("#ai_endpoint_group").hide();
    $("#ai_query_options_group").hide();
    $("#ai_provider_logo").html("");
    $("#mqc_anonymize_samples_switch").show();
    $("#ai_context_window_group").hide();
  } else if (providerId === "none") {
    $(".ai-generate-button-wrapper").hide();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_provider_info").html("Remove 'Summarize' buttons from report if you're not going to use AI summaries.");
    $("#ai_provider_logo").html("");
    $("#ai_model_group").hide();
    $("#ai_api_key_group").hide();
    $("#ai_endpoint_group").hide();
    $("#ai_query_options_group").hide();
    $("#mqc_anonymize_samples_switch").hide();
    $("#ai_context_window_group").hide();
  } else {
    // OpenAI, Anthropic, custom - showing logos and more helpful guides
    $(".ai-generate-button-wrapper").show();
    $(".ai-copy-button-wrapper").hide();
    $("#ai_provider_info").html("");
    $("#ai_provider_logo").html("");
    $("#ai_model_group").show();
    $("#ai_api_key_group").show();
    $("#mqc_anonymize_samples_switch").show();
    if (providerId === "custom") {
      $("#ai_endpoint_group").show();
      $("#ai_query_options_group").show();
      $("#ai_context_window_group").show();
      aiApiKeyInfo = `<span id="ai_api_key_info_required">Field is required.</span></a>`;
      $("#ai_model_info").hide();
      $("#ai_api_key_info").hide();
    } else {
      $("#ai_endpoint_group").hide();
      $("#ai_query_options_group").hide();
      $("#ai_context_window_group").hide();
      aiApiKeyInfo = `<span id="ai_api_key_info_required">Field is required.</span> You can find your API key in the <a style='text-decoration: underline;' href='${provider.apiKeysUrl}' target='_blank'>${provider.name} console</a>`;
      // Add clickable model suggestions if available
      let suggestedModels = provider.suggestedModels || [];
      aiModelInfo = `You can find the available models for ${provider.name} in the <a style='text-decoration: underline;' href='${provider.modelsUrl}' target='_blank'>${provider.name} documentation</a>.`;
      if (suggestedModels.length > 0) {
        aiModelInfo += " For example: ";
        aiModelInfo += suggestedModels
          .map((model) => `<a href="#" class="ai-model-suggestion" data-model="${model}"><code>${model}</code></a>`)
          .join(", ");
      }
      $("#ai_model_info").html(aiModelInfo).show();
      $("#ai_api_key_info").html(aiApiKeyInfo).show();
    }
    // Doing it here again because model depends on provider
    const storedModel = getStoredModelName(providerId);
    const defaultModel = provider.defaultModel;
    $("#ai-model").val(storedModel || defaultModel);

    if (providerId === "openai") {
      $("#ai_provider_logo").html(
        '<svg width="112" viewBox="0 0 1180 320" xmlns="http://www.w3.org/2000/svg"><path d="m367.44 153.84c0 52.32 33.6 88.8 80.16 88.8s80.16-36.48 80.16-88.8-33.6-88.8-80.16-88.8-80.16 36.48-80.16 88.8zm129.6 0c0 37.44-20.4 61.68-49.44 61.68s-49.44-24.24-49.44-61.68 20.4-61.68 49.44-61.68 49.44 24.24 49.44 61.68z"/><path d="m614.27 242.64c35.28 0 55.44-29.76 55.44-65.52s-20.16-65.52-55.44-65.52c-16.32 0-28.32 6.48-36.24 15.84v-13.44h-28.8v169.2h28.8v-56.4c7.92 9.36 19.92 15.84 36.24 15.84zm-36.96-69.12c0-23.76 13.44-36.72 31.2-36.72 20.88 0 32.16 16.32 32.16 40.32s-11.28 40.32-32.16 40.32c-17.76 0-31.2-13.2-31.2-36.48z"/><path d="m747.65 242.64c25.2 0 45.12-13.2 54-35.28l-24.72-9.36c-3.84 12.96-15.12 20.16-29.28 20.16-18.48 0-31.44-13.2-33.6-34.8h88.32v-9.6c0-34.56-19.44-62.16-55.92-62.16s-60 28.56-60 65.52c0 38.88 25.2 65.52 61.2 65.52zm-1.44-106.8c18.24 0 26.88 12 27.12 25.92h-57.84c4.32-17.04 15.84-25.92 30.72-25.92z"/><path d="m823.98 240h28.8v-73.92c0-18 13.2-27.6 26.16-27.6 15.84 0 22.08 11.28 22.08 26.88v74.64h28.8v-83.04c0-27.12-15.84-45.36-42.24-45.36-16.32 0-27.6 7.44-34.8 15.84v-13.44h-28.8z"/><path d="m1014.17 67.68-65.28 172.32h30.48l14.64-39.36h74.4l14.88 39.36h30.96l-65.28-172.32zm16.8 34.08 27.36 72h-54.24z"/><path d="m1163.69 68.18h-30.72v172.32h30.72z"/><path d="m297.06 130.97c7.26-21.79 4.76-45.66-6.85-65.48-17.46-30.4-52.56-46.04-86.84-38.68-15.25-17.18-37.16-26.95-60.13-26.81-35.04-.08-66.13 22.48-76.91 55.82-22.51 4.61-41.94 18.7-53.31 38.67-17.59 30.32-13.58 68.54 9.92 94.54-7.26 21.79-4.76 45.66 6.85 65.48 17.46 30.4 52.56 46.04 86.84 38.68 15.24 17.18 37.16 26.95 60.13 26.8 35.06.09 66.16-22.49 76.94-55.86 22.51-4.61 41.94-18.7 53.31-38.67 17.57-30.32 13.55-68.51-9.94-94.51zm-120.28 168.11c-14.03.02-27.62-4.89-38.39-13.88.49-.26 1.34-.73 1.89-1.07l63.72-36.8c3.26-1.85 5.26-5.32 5.24-9.07v-89.83l26.93 15.55c.29.14.48.42.52.74v74.39c-.04 33.08-26.83 59.9-59.91 59.97zm-128.84-55.03c-7.03-12.14-9.56-26.37-7.15-40.18.47.28 1.3.79 1.89 1.13l63.72 36.8c3.23 1.89 7.23 1.89 10.47 0l77.79-44.92v31.1c.02.32-.13.63-.38.83l-64.41 37.19c-28.69 16.52-65.33 6.7-81.92-21.95zm-16.77-139.09c7-12.16 18.05-21.46 31.21-26.29 0 .55-.03 1.52-.03 2.2v73.61c-.02 3.74 1.98 7.21 5.23 9.06l77.79 44.91-26.93 15.55c-.27.18-.61.21-.91.08l-64.42-37.22c-28.63-16.58-38.45-53.21-21.95-81.89zm221.26 51.49-77.79-44.92 26.93-15.54c.27-.18.61-.21.91-.08l64.42 37.19c28.68 16.57 38.51 53.26 21.94 81.94-7.01 12.14-18.05 21.44-31.2 26.28v-75.81c.03-3.74-1.96-7.2-5.2-9.06zm26.8-40.34c-.47-.29-1.3-.79-1.89-1.13l-63.72-36.8c-3.23-1.89-7.23-1.89-10.47 0l-77.79 44.92v-31.1c-.02-.32.13-.63.38-.83l64.41-37.16c28.69-16.55 65.37-6.7 81.91 22 6.99 12.12 9.52 26.31 7.15 40.1zm-168.51 55.43-26.94-15.55c-.29-.14-.48-.42-.52-.74v-74.39c.02-33.12 26.89-59.96 60.01-59.94 14.01 0 27.57 4.92 38.34 13.88-.49.26-1.33.73-1.89 1.07l-63.72 36.8c-3.26 1.85-5.26 5.31-5.24 9.06l-.04 89.79zm14.63-31.54 34.65-20.01 34.65 20v40.01l-34.65 20-34.65-20z"/></svg>',
      );
    } else if (providerId === "anthropic") {
      $("#ai_provider_logo").html(
        '<svg width="112" version="1.1" id="Layer_1" xmlns:x="ns_extend;" xmlns:i="ns_ai;" xmlns:graph="ns_graphs;" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" viewBox="0 0 578.9 65" style="enable-background:new 0 0 578.9 65;" xml:space="preserve"><style type="text/css">.st0{fill:#181818;}</style><metadata><sfw xmlns="ns_sfw;"><slices></slices><sliceSourceBounds bottomLeftOrigin="true" height="65" width="578.9" x="-307.7" y="-207.8"></sliceSourceBounds></sfw></metadata><g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,18.299999237060547,0.27000001072883606)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M99.6,44.8l-28.3-44H56v62.8h13v-44l28.3,44h15.3V0.8h-13V44.8L99.6,44.8z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,34.869998931884766,0.27000001072883606)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M106.8,12.9h21.1v50.7h13.5V12.9h21.1V0.8h-55.7V12.9L106.8,12.9z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,51.22999954223633,0.27000001072883606)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M200,25.9h-29.6v-25h-13.5v62.8h13.5V38H200v25.7h13.5V0.8H200V25.9L200,25.9z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,69.23999786376953,0.27000001072883606)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M225.5,12.9h16.6c6.6,0,10.1,2.4,10.1,7c0,4.6-3.5,7-10.1,7h-16.6V12.9L225.5,12.9z M265.7,20c0-11.9-8.7-19.1-23-19.1H212v62.8h13.5V39.1h15L254,63.7h14.9L254,37.2C261.5,34.3,265.7,28.3,265.7,20L265.7,20z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,84.98999786376953,0)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M291.2,52.4c-10.6,0-17.1-7.5-17.1-19.8c0-12.5,6.5-20,17.1-20c10.5,0,16.9,7.5,16.9,20C308.1,44.9,301.7,52.4,291.2,52.4L291.2,52.4z M291.2,0c-18.1,0-31,13.5-31,32.6c0,18.9,12.8,32.4,31,32.4c18,0,30.8-13.5,30.8-32.4C322,13.5,309.3,0,291.2,0L291.2,0z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,103.29000091552734,0.27000001072883606)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M346.4,28.7h-16.6V12.9h16.6c6.6,0,10.1,2.7,10.1,7.9S353.1,28.7,346.4,28.7L346.4,28.7z M347,0.8h-30.7v62.8h13.5V40.9H347c14.3,0,23-7.5,23-20C370,8.4,361.3,0.8,347,0.8L347,0.8z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,128.0399932861328,0)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M436.5,42.8c-2.3,6.1-7,9.6-13.4,9.6c-10.6,0-17.1-7.5-17.1-19.8c0-12.5,6.5-20,17.1-20c6.4,0,11,3.5,13.4,9.6h14.3C447.2,8.7,436.7,0,423.1,0c-18.1,0-31,13.5-31,32.6c0,18.9,12.8,32.4,31,32.4c13.7,0,24.2-8.8,27.7-22.2H436.5L436.5,42.8z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,117.83000183105469,0.27000001072883606)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M360.9,0.8l25.1,62.8h13.7L374.6,0.8H360.9L360.9,0.8z"></path></g></g></g><g transform="matrix(1,0,0,1,0,0)"><g transform="matrix(1,0,0,1,0,0.27000001072883606)"><g transform="matrix(1,0,0,1,0,0)"><path class="st0" d="M23.7,38.8l8.6-22.1l8.6,22.1H23.7L23.7,38.8z M25.1,0.8L0,63.7h14l5.1-13.2h26.2l5.1,13.2h14L39.4,0.8H25.1L25.1,0.8z"></path></g></g></g></g></svg>',
      );
    } else {
      $("#ai_provider_logo").html("");
    }
  }
}

//////////////////////////////////////////////////////
// AI settings handlers
//////////////////////////////////////////////////////
$(function () {
  // Populate provider dropdown dynamically
  const aiProviderSelect = $("#ai-provider");
  aiProviderSelect.empty();
  Object.entries(AI_PROVIDER_GROUPS).forEach(([groupName, groupProviders]) => {
    let optgroup = $("<optgroup>", { label: groupName });
    groupProviders.forEach((groupProviderID) => {
      optgroup.append(
        $("<option>", {
          value: groupProviderID,
          text: AI_PROVIDERS[groupProviderID].name,
        }),
      );
    });
    aiProviderSelect.append(optgroup);
  });

  // Set initial values from storage or values from Python
  const providerId = getStoredProvider() || aiConfigProviderId || "seqera";
  aiProviderSelect.val(providerId);
  const provider = AI_PROVIDERS[providerId];
  $("#ai-api-key").val(getStoredApiKey(providerId) || "");

  let model = getStoredModelName(providerId);
  if (model === null && aiConfigModel !== "None") model = aiConfigModel;
  if (model === null && provider.defaultModel) model = provider.defaultModel;
  $("#ai-model").val(model);

  let endpoint = getStoredEndpoint();
  if (endpoint === null) endpoint = aiConfigCustomEndpoint;
  $("#ai-endpoint").val(endpoint);

  let storedOptions = getStoredQueryOptions();
  if (storedOptions === null) storedOptions = aiConfigExtraQueryOptions;
  $("#ai-query-options").val(JSON.stringify(storedOptions));

  let storedContextWindow = getStoredContextWindow();
  if (storedContextWindow === null && aiConfigCustomContextWindow !== "None") {
    try {
      storedContextWindow = parseInt(aiConfigCustomContextWindow);
    } catch {}
  }
  storedContextWindow = storedContextWindow ?? 128000;
  $("#ai-context-window").val(storedContextWindow);

  updatePanel(providerId);

  // Save provider changes
  aiProviderSelect.change(function () {
    const selectedProviderId = $(this).val();
    storeProvider(selectedProviderId);

    // Update API key field
    const storedKey = getStoredApiKey(selectedProviderId);
    $("#ai-api-key").val(storedKey || "");

    updatePanel(selectedProviderId);
  });

  // Save model changes
  $("#ai-model").change(function () {
    const providerId = $("#ai-provider").val();
    const model = $(this).val();
    storeModelName(providerId, model);
  });

  $("#ai-endpoint").change(function () {
    const endpoint = $(this).val();
    storeEndpoint(endpoint);
  });

  // Save API key changes
  $("#ai-api-key").change(function () {
    const providerId = $("#ai-provider").val();
    const apiKey = $(this).val();
    storeApiKey(providerId, apiKey);
  });

  // Save query options
  $("#ai-query-options").change(function () {
    const queryOptions = $(this).val() || "{}";
    try {
      let str = JSON.parse(queryOptions);
      storeQueryOptions(str);
    } catch {}
  });

  // Also save on blur (when field loses focus)
  $("#ai-model, #ai-endpoint, #ai-api-key").blur(function () {
    $(this).trigger("change");
  });

  // Add click handlers for model suggestions (using event delegation for dynamic content)
  $(document).on("click", ".ai-model-suggestion", function (e) {
    e.preventDefault();
    const modelName = $(this).data("model");
    $("#ai-model").val(modelName).trigger("change");
  });

  // Initialize anonymize samples switch
  const anonymizeSamplesEnabled = getStoredSampleAnonymizationEnabled() ?? aiAnonymizeSamples === "true";
  if (anonymizeSamplesEnabled) {
    $(".mqc_switch.anonymize_samples").removeClass("off").addClass("on").text("on");
  }

  // Handle anonymize samples switch clicks
  $(".mqc_switch_wrapper.mqc_anonymize_samples").click(function (e) {
    e.preventDefault();
    const switchElem = $(this).find(".mqc_switch");
    const isEnabled = switchElem.hasClass("on");
    if (isEnabled) {
      switchElem.removeClass("on").addClass("off").text("off");
      storeSampleAnonymizationEnabled(false);
    } else {
      switchElem.removeClass("off").addClass("on").text("on");
      storeSampleAnonymizationEnabled(true);
    }
  });
});

// Read the JWT local cookie in in the seqera.io domain - that's our API key:
$(function () {
  const jwt = document.cookie.split("; ").find((row) => row.startsWith("jwt="));
  if (jwt) {
    const jwtValue = jwt.split("=")[1];
    saveToLocalStorage(`mqc_ai_key_${provider}`, jwtValue);
  }
});

// Storing user settings
function getStoredProvider() {
  let storedProviderId = getFromLocalStorage("mqc_ai_provider");
  if (storedProviderId && AI_PROVIDERS[storedProviderId]) return storedProviderId;
  return null;
}
function storeProvider(providerId) {
  saveToLocalStorage("mqc_ai_provider", providerId);
}
function getStoredApiKey(providerId) {
  return getFromLocalStorage(`mqc_ai_key_${providerId}`);
}
function storeApiKey(providerId, apiKey) {
  saveToLocalStorage(`mqc_ai_key_${providerId}`, apiKey);
}
function getStoredModelName(providerId) {
  return getFromLocalStorage(`mqc_ai_model_${providerId}`);
}
function storeModelName(providerId, modelName) {
  saveToLocalStorage(`mqc_ai_model_${providerId}`, modelName);
}
function getStoredEndpoint() {
  return getFromLocalStorage(`mqc_ai_endpoint`);
}
function storeEndpoint(endpoint) {
  saveToLocalStorage(`mqc_ai_endpoint`, endpoint);
}
function getStoredQueryOptions() {
  let data = getFromLocalStorage(`mqc_ai_query_options`);
  if (data) {
    try {
      return JSON.parse(data);
    } catch (e) {
      console.error("Error parsing extra query options", e);
      return null;
    }
  }
  return null;
}
function storeQueryOptions(options) {
  saveToLocalStorage(`mqc_ai_query_options`, JSON.stringify(options));
}
function getStoredSampleAnonymizationEnabled() {
  const stored = getFromLocalStorage("mqc_ai_anonymize_samples");
  return stored === null ? null : stored === "true";
}
function storeSampleAnonymizationEnabled(value) {
  saveToLocalStorage("mqc_ai_anonymize_samples", value.toString());
}
function getStoredContextWindow() {
  return getFromLocalStorage(`mqc_ai_context_window`);
}
function storeContextWindow(value) {
  saveToLocalStorage(`mqc_ai_context_window`, value.toString());
}

function addLogo(imageDataUrl, callback) {
  // Append "watermark" to the image
  let plotlyImage = new Image();
  plotlyImage.onload = function () {
    let canvas = document.createElement("canvas");
    let ctx = canvas.getContext("2d");

    // Set canvas size to double for retina display
    canvas.width = plotlyImage.width;
    canvas.height = plotlyImage.height; // additional height for the text

    ctx.drawImage(plotlyImage, 0, 0, plotlyImage.width, plotlyImage.height);

    // Text properties
    ctx.font = "9px Lucida Grande"; // Set the desired font-size and type
    ctx.textAlign = "right";
    ctx.fillStyle = "#9f9f9f"; // Semi-transparent black text

    ctx.fillText("Created with MultiQC", plotlyImage.width - 15, plotlyImage.height - 6);
    // Callback with the combined image
    callback(canvas.toDataURL());
  };
  plotlyImage.src = imageDataUrl;
}
