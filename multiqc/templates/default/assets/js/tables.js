////////////////////////////////////////////////
// MultiQC Table code
////////////////////////////////////////////////

// Execute when page load has finished loading
$(function () {
  if ($(".mqc_per_sample_table").length > 0) {
    // Enable tablesorter on MultiQC tables
    let getSortVal = function (node) {
      // if val is defined, use it
      let val = $(node).data("sorting-val");
      if (val !== null && val !== undefined) {
        if (val === "") {
          return val;
        }
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
    $(".mqc_per_sample_table").tablesorter({
      sortInitialOrder: "desc",
      textExtraction: getSortVal,
      cancelSelection: false,
      headers: null, // can revert when https://github.com/Mottie/tablesorter/pull/1851 is merged
    });

    // Update tablesorter if samples renamed
    $(document).on("mqc_renamesamples", function (e, f_texts, t_texts, regex_mode) {
      $(".mqc_per_sample_table").trigger("update");
    });

    $(".mqc-table-to-violin").click(function (e) {
      e.preventDefault();
      let tableAnchor = $(this).data("table-anchor");
      let violinAnchor = $(this).data("violin-anchor");
      $("#mqc_violintable_wrapper_" + tableAnchor).hide();
      $("#mqc_violintable_wrapper_" + violinAnchor).show();
      renderPlot(violinAnchor);
    });

    $(".mqc-violin-to-table").click(function (e) {
      e.preventDefault();
      let tableAnchor = $(this).data("table-anchor");
      let violinAnchor = $(this).data("violin-anchor");
      $("#mqc_violintable_wrapper_" + tableAnchor).show();
      $("#mqc_violintable_wrapper_" + violinAnchor).hide();
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
      let tableAnchor = $(this).data("table-anchor");
      let violinAnchor = $(this).data("violin-anchor");
      mqc_table_col_updateVisible(tableAnchor, violinAnchor);
    });
    // Bulk set visible / hidden
    $(".mqc_config_modal_bulk_visible").click(function (e) {
      e.preventDefault();
      let tableAnchor = $(this).data("table-anchor");
      let violinAnchor = $(this).data("violin-anchor");
      let visible = $(this).data("action") === "showAll";
      $("#" + tableAnchor + "_config_modal_table tbody .mqc_table_col_visible").prop("checked", visible);
      mqc_table_col_updateVisible(tableAnchor, violinAnchor);
    });
    function mqc_table_col_updateVisible(tableAnchor, violinAnchor) {
      let target = "#" + tableAnchor;

      let metricsHidden = {};
      $(target + "_config_modal_table .mqc_table_col_visible").each(function () {
        let metric = $(this).val();
        metricsHidden[metric] = !$(this).is(":checked");
      });

      Object.entries(metricsHidden).map(([metric, hidden]) => {
        if (hidden) {
          $(target + " ." + metric).addClass("column-hidden");
          $(target + "_config_modal_table ." + metric).addClass("text-muted");
        } else {
          $(target + " ." + metric).removeClass("column-hidden");
          $(target + "_config_modal_table ." + metric).removeClass("text-muted");
        }
      });
      // Hide empty rows
      $(target + " tbody tr").each(function () {
        let trIsEmpty = true;
        let tr = $(this);
        tr.find("td").each(function () {
          let td = $(this);
          if (!td.hasClass("column-hidden") && !td.hasClass("sorthandle") && td.text() !== "") {
            trIsEmpty = false;
          }
        });
        if (trIsEmpty) {
          tr.addClass("row-empty");
        } else {
          tr.removeClass("row-empty");
        }
      });
      // Update counts
      $(target + "_numrows").text($(target + " tbody tr:visible").length);
      $(target + "_numcols").text($(target + " thead th:visible").length - 1);

      // Also update the violin plot
      if (violinAnchor !== undefined) {
        let plot = mqc_plots[violinAnchor];
        plot.datasets.map((dataset) => {
          dataset["metrics"].map((metric) => {
            dataset["header_by_metric"][metric]["hidden"] = metricsHidden[metric];
          });
        });
        renderPlot(violinAnchor);
      }
    }

    // Make rows in MultiQC "Configure Columns" tables sortable
    $(".mqc_config_modal").on("show.bs.modal", function (e) {
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
    $(".mqc_config_modal_table").on("sortstop", function (e, ui) {
      change_mqc_table_col_order($(this));
    });
    $(".mqc_config_modal_table").bind("sortEnd", function () {
      change_mqc_table_col_order($(this));
    });

    // TOOLBOX LISTENERS

    // highlight samples
    $(document).on("mqc_highlights", function (e, f_texts, f_cols, regex_mode) {
      $(".mqc_table_sortHighlight").hide();
      $(".mqc_per_sample_table tbody th").removeClass("highlighted").removeData("highlight");
      $(".mqc_per_sample_table tbody th").each(function (i) {
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
      let tableAnchor = $(this).data("table-anchor");
      // collect highlighted rows
      let hrows = $("#" + tableAnchor + " tbody th.highlighted")
        .parent()
        .detach();
      hrows = hrows.sort(function (a, b) {
        return $(a).find("th").data("highlight") - $(b).find("th").data("highlight");
      });
      if ($(this).data("direction") === "desc") {
        hrows = hrows.get().reverse();
        $("#" + tableAnchor + " tbody").prepend(hrows);
        $(this).data("direction", "asc");
      } else {
        $("#" + tableAnchor + " tbody").append(hrows);
        $(this).data("direction", "desc");
      }
    });

    // Rename samples
    $(document).on("mqc_renamesamples", function (e, f_texts, t_texts, regex_mode) {
      $(".mqc_per_sample_table tbody th span.th-sample-name").each(function () {
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
      $(".mqc_per_sample_table tbody tr").each(function () {
        let tr = $(this);
        let th = tr.find("th");
        let match = false;
        let s_name = th.find(".th-sample-name").text();
        let g_name = tr.data("sample-group");
        $.each(f_texts, function (idx, f_text) {
          if (regex_mode) {
            if (s_name.match(f_text) || g_name.match(f_text)) {
              match = true;
            }
          } else {
            if (s_name.indexOf(f_text) > -1 || g_name.indexOf(f_text) > -1) {
              match = true;
            }
          }
        });
        if (window.mqc_hide_mode === "show") {
          match = !match;
        }
        if (match) {
          tr.addClass("sample-hidden");
        } else {
          tr.removeClass("sample-hidden");
        }
      });
      $(".mqc_table_numrows").each(function () {
        let tid = $(this).attr("id").replace("_numrows", "");
        $(this).text($("#" + tid + " tbody tr:visible").length);
      });

      // Hide empty columns
      $(".mqc_per_sample_table").each(function () {
        let table = $(this);
        let gsthidx = 0;
        table.find("thead th").each(function () {
          let th = $(this);
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
              let td = $(this);
              count += 1;
              if (td.text() === "") {
                empties += 1;
              }
            });
          if (count > 0 && count === empties) {
            th.addClass("column-hidden");
            table.find("tbody tr td:nth-child(" + (gsthidx + 2) + ")").addClass("column-hidden");
          }
          gsthidx += 1;
        });
      });
      $(".mqc_table_numcols").each(function () {
        let tid = $(this).attr("id").replace("_numcols", "");
        $(this).text($("#" + tid + " thead th:visible").length - 1);
      });
    });

    // Support expanding grouped samples in table
    $(".expandable-row-primary").click(function (e) {
      e.preventDefault();

      // if the user is selecting text, do not expand the row
      if (window.getSelection().toString().length > 0) {
        return;
      }
      // let th = $(this);
      // final most parent table object
      let tr = $(this);
      let table = tr.closest("table");
      let tableId = tr.data("table-id");
      // find all rows with the same data-group-id
      let groupId = tr.data("sample-group");
      let otherRows = table.find(
        "tbody tr.expandable-row-secondary[data-sample-group='" + groupId + "'][data-table-id='" + tableId + "']",
      );
      // toggle the visibility of the rows and type of arrow
      otherRows.toggleClass("expandable-row-secondary-hidden");
      tr.toggleClass("expanded");
    });

    // We want to allow user select sample name text without expanding the row
    $(".th-sample-name").click(function (e) {
      // stop propagation to prevent row expansion
      // e.stopPropagation();
    });
  } // End of check for table

  // Table Scatter Modal
  $("#table_scatter_form").submit(function (e) {
    e.preventDefault();
  });
  $(".mqc_table_make_scatter").click(function (e) {
    // Reset dropdowns
    let tableAnchor = $(this).data("table-anchor");
    let table_scatter_table_anchor_el = $("#table_scatter_table_anchor");
    if (table_scatter_table_anchor_el.val() !== tableAnchor) {
      $("#table_scatter_col1, #table_scatter_col2").html('<option value="">Select Column</option>');
      // Add columns to dropdowns
      $("#" + tableAnchor + " thead tr th").each(function (e) {
        let thId = $(this).attr("id");
        if (thId !== undefined) {
          let name = $(this).text();
          let namespace = $(this).data("namespace");
          if (namespace) {
            name = namespace + ": " + name;
          }
          $("#table_scatter_col1, #table_scatter_col2").append('<option value="' + thId + '">' + name + "</select>");
        }
      });
      table_scatter_table_anchor_el.val(tableAnchor);
      $("#table_scatter_plot").html("<small>Please select two table columns.</small>").addClass("not_rendered");
    }
  });
  $("#table_scatter_form select").change(function (e) {
    let tableAnchor = $("#table_scatter_table_anchor").val();
    let col1 = $("#table_scatter_col1").val().replace("header_", "");
    let col2 = $("#table_scatter_col2").val().replace("header_", "");
    let col1_name = $("#table_scatter_col1 option:selected").text();
    let col2_name = $("#table_scatter_col2 option:selected").text();
    let col1_min = parseFloat($("#" + tableAnchor + " thead th#header_" + col1).data("dmin"));
    let col1_max = parseFloat($("#" + tableAnchor + " thead th#header_" + col1).data("dmax"));
    let col2_max = parseFloat($("#" + tableAnchor + " thead th#header_" + col2).data("dmax"));
    let col2_min = parseFloat($("#" + tableAnchor + " thead th#header_" + col2).data("dmin"));
    if (!isNaN(col1_min) && !isNaN(col1_max) && col1_max != 0) {
      col1_max += (col1_max - col1_min) * 0.05;
    }
    if (!isNaN(col1_min) && !isNaN(col1_max) && col1_min != 0) {
      col1_min -= (col1_max - col1_min) * 0.05;
    }
    if (!isNaN(col2_min) && !isNaN(col2_max) && col2_max != 0) {
      col2_max += (col2_max - col2_min) * 0.05;
    }
    if (!isNaN(col2_min) && !isNaN(col2_max) && col2_min != 0) {
      col2_min -= (col2_max - col2_min) * 0.05;
    }
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
    let plotDiv = $("#table_scatter_plot");
    if (col1 !== "" && col2 !== "") {
      plotDiv.html("<small>loading..</small>");
      let plotTitle;
      if ($("#" + tableAnchor).data("title")) {
        plotTitle = $("#" + tableAnchor).data("title");
      } else {
        plotTitle = tableAnchor.replace(/_/g, " ");
      }
      let plotDataset = [];
      $("#" + tableAnchor + " tbody tr").each(function (e) {
        let sName = $(this).children("th.rowheader").text();
        let val_1 = $(this)
          .children("td." + col1)
          .data("sorting-val");
        let val_2 = $(this)
          .children("td." + col2)
          .data("sorting-val");
        if (!isNaN(parseFloat(val_1)) && isFinite(val_1) && !isNaN(parseFloat(val_2)) && isFinite(val_2)) {
          plotDataset.push({
            name: sName,
            x: parseFloat(val_1),
            y: parseFloat(val_2),
          });
        }
      });
      if (Object.keys(plotDataset).length > 0) {
        let target = "table_scatter_plot";
        let traces = plotDataset.map(function (point) {
          return {
            type: "scatter",
            x: [point.x],
            y: [point.y],
            name: point.name,
            text: [point.name],
          };
        });
        let layout = {
          title: plotTitle,
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
          plotDiv.html("<small>Error: Something went wrong when plotting the scatter plot.</small>");
          plotDiv.addClass("not_rendered");
        } else {
          plotDiv.removeClass("not_rendered");
        }
      } else {
        plotDiv.html("<small>Error: No data pairs found for these columns.</small>");
        plotDiv.addClass("not_rendered");
      }
    } else {
      plotDiv.html("<small>Please select two table columns.</small>");
      plotDiv.addClass("not_rendered");
    }
  });
});

// Reorder columns in MultiQC tables.
// Note: Don't have to worry about floating headers, as 'Configure Columns'
// button is only visible when this is hidden. Ace!
function change_mqc_table_col_order(table) {
  // Find the targets of this sorting
  let elemId = table.attr("id");
  let tableAnchor = elemId.replace("_config_modal_table", "");

  // Collect the desired order of columns
  let classes = [];
  $("#" + elemId + " tbody tr").each(function () {
    classes.push($(this).attr("class"));
  });
  // Go through each row
  $("#" + tableAnchor + " tr").each(function () {
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
