////////////////////////////////////////////////
// MultiQC Report Toolbox - Export Functionality (ECharts backend)
//
// Copy of multiqc/templates/default/src/js/toolbox/export.js with only the
// Plotly image-export block replaced by an ECharts equivalent. See
// "In-report PNG export" / "export finding" in
// multiqc-echarts-exploration/BUILD_PLAN.md and 03-interactivity-testing.md:
// a single live ECharts instance cannot emit both PNG and SVG (getDataURL's
// output format is locked to the renderer, not the `type` arg), so PNG uses a
// throwaway canvas-renderer instance and SVG uses a throwaway SSR svg-renderer
// instance via renderToSVGString().
////////////////////////////////////////////////

// Render `plot`'s current option (same assembly renderPlot() uses, via the shared
// buildCurrentOption() in echarts-plotting.js, so exports can never drift from what's
// on screen) into a data URL / SVG string.
function echartsImage(plot, format, width, height, scale) {
  let option = window.buildCurrentOption(plot);
  if (!option) {
    return Promise.reject(new Error("Cannot export " + plot.anchor + ": all series are hidden."));
  }

  if (format === "svg") {
    let chart = echarts.init(null, null, { renderer: "svg", ssr: true, width: width, height: height });
    chart.setOption(option);
    let svg = chart.renderToSVGString();
    chart.dispose();
    return Promise.resolve(new Blob([svg], { type: "image/svg+xml;charset=utf-8" }));
  }

  // PNG: needs a real DOM element to init a canvas-renderer instance against, so it can't
  // use ssr:true like the SVG path. Render into a hidden, offscreen container instead.
  let el = document.createElement("div");
  el.style.position = "fixed";
  el.style.top = "-10000px";
  el.style.left = "-10000px";
  el.style.width = width + "px";
  el.style.height = height + "px";
  document.body.appendChild(el);

  let chart = echarts.init(el, null, { renderer: "canvas" });
  chart.setOption(option);
  let dataUrl = chart.getDataURL({ type: "png", pixelRatio: scale });
  chart.dispose();
  document.body.removeChild(el);
  return Promise.resolve(dataUrl);
}

// Make functions available globally
window.initExport = function () {
  // Change text on download button
  $('#mqc_exportplots a[data-bs-toggle="tab"]').on("shown.bs.tab", function (e) {
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
        `<div class="form-check">
          <input type="checkbox" class="form-check-input" value="${fname}" id="mqc_export_plot_${fname}" checked>
          <label class="form-check-label" for="mqc_export_plot_${fname}">${fname}</label>
        </div>`,
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

    // Export the plots
    $("#mqc_exportplots").submit(function (e) {
      e.preventDefault();
      let checked_plots = $("#mqc_export_selectplots input:checked");
      let zip = new JSZip();
      let promises = [];
      //////
      // EXPORT PLOT IMAGES
      //////
      if ($("#mqc_image_download").is(":visible")) {
        let mime = $("#mqc_export_ft").val();
        let format = mime.replace("image/", "").split("+")[0];
        let f_width = parseInt($("#mqc_exp_width").val());
        let f_height = parseInt($("#mqc_exp_height").val());
        const font_scale = parseFloat($("#mqc_export_scaling").val());
        checked_plots.each(function () {
          const target = $(this).val();
          let plot = mqc_plots[target];
          if (plot === undefined) return;

          promises.push(
            echartsImage(plot, format, f_width / font_scale, f_height / font_scale, font_scale).then(function (img) {
              if (format === "svg") {
                // img is a Blob (image/svg+xml) here, unlike the Plotly path's data URL.
                if (checked_plots.length <= zip_threshold) {
                  // Not many plots to export, just trigger a download for each
                  saveAs(img, target + "." + format);
                } else {
                  // Lots of plots - add to a zip file for download
                  zip.file(target + "." + format, img);
                }
              } else {
                // Can add logo to a PNG image
                return new Promise(function (resolve) {
                  addLogo(img, function (imageWithLogo) {
                    if (checked_plots.length <= zip_threshold) {
                      // Not many plots to export, just trigger a download for each:
                      const blob = dataUrlToBlob(imageWithLogo, mime);
                      saveAs(blob, target + "." + format);
                    } else {
                      // Lots of plots - add to a zip file for download:
                      const fname = target + "." + format;
                      const data = imageWithLogo.replace(/^data:image\/png;base64,/, "");
                      zip.file(fname, data, { base64: true });
                    }
                    resolve();
                  });
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
      // EXPORT PLOT DATA
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
};
