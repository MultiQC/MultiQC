////////////////////////////////////////////////
// MultiQC Table code
////////////////////////////////////////////////

// Execute when page load has finished loading
$(function () {
  if ($(".mqc_table").length > 0) {
    // Enable tablesorter on MultiQC tables
    let get_sort_val = function (node) {
      // if val is defined, use it
      if (node.getAttribute("val") !== null) {
        let val = node.getAttribute("val");
        let floatVal = parseFloat(val);
        if (!isNaN(floatVal)) {
          return val; // expected to return a string
        }
        return val;
      }

      const text = node.innerText;

      // If first char is a digit, strip non-numeric
      // This is to handle cases of e.g. 300X and 9.0X.
      if (text.length > 0 && text[0].match(/\d/)) {
        return text.replace(/[^\d.]/g, "");
      }

      return text;
    };
    $(".mqc_table").tablesorter({
      sortInitialOrder: "desc",
      textExtraction: get_sort_val,
      cancelSelection: false,
      headers: null, // can revert when https://github.com/Mottie/tablesorter/pull/1851 is merged
    });

    // Update tablesorter if samples renamed
    $(document).on("mqc_renamesamples", function (e, f_texts, t_texts, regex_mode) {
      $(".mqc_table").trigger("update");
    });

    $(".mqc-table-to-violin").click(function (e) {
      e.preventDefault();
      let tableId = $(this).data("table-id");
      let violinId = $(this).data("violin-id");
      $("#mqc_violintable_wrapper_" + tableId).hide();
      $("#mqc_violintable_wrapper_" + violinId).show();
      renderPlot(violinId);
    });

    $(".mqc-violin-to-table").click(function (e) {
      e.preventDefault();
      let tableId = $(this).data("table-id");
      let violinId = $(this).data("violin-id");
      $("#mqc_violintable_wrapper_" + tableId).show();
      $("#mqc_violintable_wrapper_" + violinId).hide();
    });

    $(".mqc_table_copy_btn").click(function () {
      let btn = $(this);
      let table = $(btn.data("clipboard-target"))[0];

      const range = document.createRange();
      range.selectNode(table);
      window.getSelection().removeAllRanges();
      window.getSelection().addRange(range);

      try {
        document.execCommand("copy");
        window.getSelection().removeAllRanges();

        btn.addClass("active").html('<span class="glyphicon glyphicon-copy"></span> Copied!');
        setTimeout(() => {
          btn.removeClass("active").html('<span class="glyphicon glyphicon-copy"></span> Copied!');
        }, 2000);
      } catch (err) {
        console.error("Failed to copy table: ", err);
      }
    });

    // Make table headers fixed when table body scrolls (use CSS transforms)
    // http://stackoverflow.com/a/25902860/713980
    $(".mqc-table-responsive").scroll(function () {
      $(this)
        .find("thead")
        .css("transform", "translate(0," + $(this).scrollTop() + "px)");
    });

    // Table header-specific bootstrap tooltips
    $(".mqc_table_tooltip").tooltip({ container: "body" });

    // Expand tables to full height
    $(".mqc-table-expand").click(function () {
      if ($(this).find("span").hasClass("glyphicon-chevron-down")) {
        $(this).parent().find(".mqc-table-responsive").css("max-height", "none");
        $(this).find("span").removeClass("glyphicon-chevron-down").addClass("glyphicon-chevron-up");
      } else {
        $(this).parent().find(".mqc-table-responsive").css("max-height", "400px");
        $(this).find("span").removeClass("glyphicon-chevron-down").addClass("glyphicon-chevron-down");
      }
    });

    /////// COLUMN CONFIG
    // show + hide columns
    $(".mqc_table_col_visible").change(function () {
      let tableId = $(this).data("table-id");
      let violinId = $(this).data("violin-id");
      mqc_table_col_updateVisible(tableId, violinId);
    });
    // Bulk set visible / hidden
    $(".mqc_configModal_bulkVisible").click(function (e) {
      e.preventDefault();
      let tableId = $(this).data("table-id");
      let violinId = $(this).data("violin-id");
      let visible = $(this).data("action") === "showAll";
      $("#" + tableId + "_configModal_table tbody .mqc_table_col_visible").prop("checked", visible);
      mqc_table_col_updateVisible(tableId, violinId);
    });
    function mqc_table_col_updateVisible(tableId, violinId) {
      let target = "#" + tableId;

      let metricsHidden = {};
      $(target + "_configModal_table .mqc_table_col_visible").each(function () {
        let metric = $(this).val();
        metricsHidden[metric] = !$(this).is(":checked");
      });

      Object.entries(metricsHidden).map(([metric, hidden]) => {
        if (hidden) {
          $(target + " ." + metric).addClass("hidden");
          $(target + "_configModal_table ." + metric).addClass("text-muted");
        } else {
          $(target + " ." + metric).removeClass("hidden");
          $(target + "_configModal_table ." + metric).removeClass("text-muted");
        }
      });
      // Hide empty rows
      $(target + " tbody tr").show();
      $(target + " tbody tr").each(function () {
        let hasVal = false;
        $(this)
          .find("td")
          .each(function () {
            if (!$(this).hasClass("sorthandle") && $(this).text() !== "") {
              hasVal = true;
            }
          });
        if (!hasVal) {
          $(this).hide();
        }
      });
      // Update counts
      $(target + "_numrows").text($(target + " tbody tr:visible").length);
      $(target + "_numcols").text($(target + " thead th:visible").length - 1);

      // Also update the violin plot
      if (violinId !== undefined) {
        let plot = mqc_plots[violinId];
        plot.datasets.map((dataset) => {
          dataset["metrics"].map((metric) => {
            dataset["header_by_metric"][metric]["hidden"] = metricsHidden[metric];
          });
        });
        renderPlot(violinId);
      }
    }

    // Make rows in MultiQC "Configure Columns" tables sortable
    $(".mqc_configModal").on("show.bs.modal", function (e) {
      $(e.target)
        .find(".mqc_table.mqc_sortable tbody")
        .sortable({
          handle: ".sorthandle",
          helper: function fixWidthHelper(e, ui) {
            ui.children().each(function () {
              $(this).width($(this).width());
            });
            return ui;
          },
        });
    });

    // Change order of columns
    $(".mqc_configModal_table").on("sortstop", function (e, ui) {
      change_mqc_table_col_order($(this));
    });
    $(".mqc_configModal_table").bind("sortEnd", function () {
      change_mqc_table_col_order($(this));
    });

    // TOOLBOX LISTENERS

    // highlight samples
    $(document).on("mqc_highlights", function (e, f_texts, f_cols, regex_mode) {
      $(".mqc_table_sortHighlight").hide();
      $(".mqc_table tbody th").removeClass("highlighted").removeData("highlight");
      $(".mqc_table tbody th").each(function (i) {
        let th = $(this);
        let thtext = $(this).text();
        let thiscol = "#333";
        $.each(f_texts, function (idx, f_text) {
          if ((regex_mode && thtext.match(f_text)) || (!regex_mode && thtext.indexOf(f_text) > -1)) {
            thiscol = f_cols[idx] ?? "#cccccc";
            th.addClass("highlighted").data("highlight", idx);
            $(".mqc_table_sortHighlight").show();
          }
        });
        $(this).css("color", thiscol);
      });
    });

    // Sort MultiQC tables by highlight
    $(".mqc_table_sortHighlight").click(function (e) {
      e.preventDefault();
      let target = $(this).data("target");
      // collect highlighted rows
      let hrows = $(target + " tbody th.highlighted")
        .parent()
        .detach();
      hrows = hrows.sort(function (a, b) {
        return $(a).find("th").data("highlight") - $(b).find("th").data("highlight");
      });
      if ($(this).data("direction") === "desc") {
        hrows = hrows.get().reverse();
        $(target + " tbody").prepend(hrows);
        $(this).data("direction", "asc");
      } else {
        $(target + " tbody").append(hrows);
        $(this).data("direction", "desc");
      }
    });

    // Rename samples
    $(document).on("mqc_renamesamples", function (e, f_texts, t_texts, regex_mode) {
      $(".mqc_table tbody th").each(function () {
        let s_name = String($(this).data("original-sn"));
        $.each(f_texts, function (idx, f_text) {
          if (regex_mode) {
            let re = new RegExp(f_text, "g");
            s_name = s_name.replace(re, t_texts[idx]);
          } else {
            s_name = s_name.replace(f_text, t_texts[idx]);
          }
        });
        $(this).text(s_name);
      });
    });

    // Hide samples
    $(document).on("mqc_hidesamples", function (e, f_texts, regex_mode) {
      // Hide rows in MultiQC tables
      $(".mqc_table tbody th").each(function () {
        let match = false;
        let hfilter = $(this).text();
        $.each(f_texts, function (idx, f_text) {
          if ((regex_mode && hfilter.match(f_text)) || (!regex_mode && hfilter.indexOf(f_text) > -1)) {
            match = true;
          }
        });
        if (window.mqc_hide_mode === "show") {
          match = !match;
        }
        if (match) {
          $(this).parent().hide().addClass("hidden");
        } else {
          $(this).parent().show().removeClass("hidden");
        }
      });
      $(".mqc_table_numrows").each(function () {
        let tid = $(this).attr("id").replace("_numrows", "");
        $(this).text($("#" + tid + " tbody tr:visible").length);
      });

      // Hide empty columns
      $(".mqc_table").each(function () {
        let table = $(this);
        let gsthidx = 0;
        table.find("thead th, tbody tr td").show();
        table.find("thead th").each(function () {
          if (gsthidx === 0) {
            gsthidx += 1;
            return true;
          }
          let count = 0;
          let empties = 0;
          table
            .find("tbody tr td:nth-child(" + (gsthidx + 2) + ")")
            .filter(":visible")
            .each(function () {
              count += 1;
              if ($(this).text() === "") {
                empties += 1;
              }
            });
          if (count > 0 && count === empties) {
            $(this).hide();
            table.find("tbody tr td:nth-child(" + (gsthidx + 2) + ")").hide();
          }
          gsthidx += 1;
        });
      });
      $(".mqc_table_numcols").each(function () {
        let tid = $(this).attr("id").replace("_numcols", "");
        $(this).text($("#" + tid + " thead th:visible").length - 1);
      });
    });
  } // End of check for table

  // Table Scatter Modal
  $("#tableScatterForm").submit(function (e) {
    e.preventDefault();
  });
  $(".mqc_table_makeScatter").click(function (e) {
    // Reset dropdowns
    if ($("#tableScatter_tid").val() !== $(this).data("table")) {
      $("#tableScatter_col1, #tableScatter_col2").html('<option value="">Select Column</option>');
      // Add columns to dropdowns
      $($(this).data("table") + " thead tr th").each(function (e) {
        let c_id = $(this).attr("id");
        if (c_id !== undefined) {
          let c_name = $(this).attr("data-namespace") + ": " + $(this).text();
          $("#tableScatter_col1, #tableScatter_col2").append('<option value="' + c_id + '">' + c_name + "</select>");
        }
      });
      $("#tableScatter_tid").val($(this).data("table"));
      $("#tableScatterPlot").html("<small>Please select two table columns.</small>").addClass("not_rendered");
    }
  });
  $("#tableScatterForm select").change(function (e) {
    let tid = $("#tableScatter_tid").val();
    let col1 = $("#tableScatter_col1").val().replace("header_", "");
    let col2 = $("#tableScatter_col2").val().replace("header_", "");
    let col1_name = $("#tableScatter_col1 option:selected").text();
    let col2_name = $("#tableScatter_col2 option:selected").text();
    let col1_max = parseFloat($(tid + " thead th#header_" + col1).data("dmax"));
    let col1_min = parseFloat($(tid + " thead th#header_" + col1).data("dmin"));
    let col2_max = parseFloat($(tid + " thead th#header_" + col2).data("dmax"));
    let col2_min = parseFloat($(tid + " thead th#header_" + col2).data("dmin"));
    if (isNaN(col1_max)) {
      col1_max = undefined;
    }
    if (isNaN(col1_min)) {
      col1_min = undefined;
    }
    if (isNaN(col2_max)) {
      col2_max = undefined;
    }
    if (isNaN(col2_min)) {
      col2_min = undefined;
    }
    if (col1 !== "" && col2 !== "") {
      $("#tableScatterPlot").html("<small>loading..</small>");
      if ($(tid).attr("data-title")) {
        plot_title = $(tid).attr("data-title");
      } else {
        plot_title = tid.replace(/^#/, "").replace(/_/g, " ");
      }
      // Get the data values
      mqc_plots["tableScatterPlot"] = {
        plot_type: "scatter",
        config: {
          id: "tableScatter_" + tid,
          title: plot_title,
          xlab: col1_name,
          ylab: col2_name,
          xmin: col1_min,
          xmax: col1_max,
          ymin: col2_min,
          ymax: col2_max,
        },
        datasets: [[]],
      };
      $(tid + " tbody tr").each(function (e) {
        let s_name = $(this).children("th.rowheader").text();
        let val_1 = $(this)
          .children("td." + col1)
          .text()
          .replace(/[^\d\.]/g, "");
        let val_2 = $(this)
          .children("td." + col2)
          .text()
          .replace(/[^\d\.]/g, "");
        if (!isNaN(parseFloat(val_1)) && isFinite(val_1) && !isNaN(parseFloat(val_2)) && isFinite(val_2)) {
          mqc_plots["tableScatterPlot"]["datasets"][0].push({
            name: s_name,
            x: parseFloat(val_1),
            y: parseFloat(val_2),
          });
        }
      });
      if (Object.keys(mqc_plots["tableScatterPlot"]["datasets"][0]).length > 0) {
        let target = "tableScatterPlot";
        let traces = mqc_plots["tableScatterPlot"]["datasets"][0].map(function (point) {
          return {
            type: "scatter",
            x: [point.x],
            y: [point.y],
            name: point.name,
            text: [point.name],
          };
        });
        let layout = {
          title: plot_title,
          xaxis: {
            title: col1_name,
            range: [col1_min, col1_max],
          },
          yaxis: {
            title: col2_name,
            range: [col2_min, col2_max],
          },
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
        let plot = Plotly.newPlot(target, traces, layout, config);
        if (!plot) {
          $("#tableScatterPlot").html("<small>Error: Something went wrong when plotting the scatter plot.</small>");
          $("#tableScatterPlot").addClass("not_rendered");
        } else {
          $("#tableScatterPlot").removeClass("not_rendered");
        }
      } else {
        $("#tableScatterPlot").html("<small>Error: No data pairs found for these columns.</small>");
        $("#tableScatterPlot").addClass("not_rendered");
      }
    } else {
      $("#tableScatterPlot").html("<small>Please select two table columns.</small>");
      $("#tableScatterPlot").addClass("not_rendered");
    }
  });
});

// Reorder columns in MultiQC tables.
// Note: Don't have to worry about floating headers, as 'Configure Columns'
// button is only visible when this is hidden. Ace!
function change_mqc_table_col_order(table) {
  // Find the targets of this sorting
  let elemId = table.attr("id");
  let tableId = elemId.replace("_configModal_table", "");

  // Collect the desired order of columns
  let classes = [];
  $("#" + elemId + " tbody tr").each(function () {
    classes.push($(this).attr("class"));
  });
  // Go through each row
  $("#" + tableId + " tr").each(function () {
    let cols = {};
    let row = $(this);
    // Detach any cell that matches a known class from above
    row.find("td, th").each(function () {
      let cell = $(this);
      $.each(classes, function (idx, c) {
        if (cell.hasClass(c)) {
          cols[c] = cell.detach();
        }
      });
    });
    // Insert detached cells back in the order given in the sorted table
    for (let idx in classes) {
      let c = classes[idx];
      if (cols[c] !== undefined) {
        row.append(cols[c]);
      }
    }
  });
}
