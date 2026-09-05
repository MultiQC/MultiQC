////////////////////////////////////////////////
// MultiQC Report Toolbox - Highlight Samples
////////////////////////////////////////////////

// Make functions available globally
window.apply_mqc_highlights = function () {
  // Collect the filters into an array
  var f_texts = [];
  var f_cols = [];
  var regex_mode = $("#mqc_cols .mqc_regex_mode input").prop("checked");
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
};

// Make functions available globally
window.initHighlights = function () {
  // Initialize highlight functionality

  // Highlight colour filters
  $("#mqc_color_form").submit(function (e) {
    e.preventDefault();
    let f_text = $("#mqc_colour_filter").val().trim();
    let f_col = $("#mqc_colour_filter_color").val().trim();
    $("#mqc_col_filters").append(window.make_colorsamples_filter(f_text, f_col));
    $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    $("#mqc_colour_filter").val("");
    $("#mqc_colour_filter_color").val(mqc_colours[window.mqc_colours_idx % mqc_colours.length]);
  });

  $("#mqc_cols_apply").click(function (e) {
    if (apply_mqc_highlights()) {
      $(this).attr("disabled", true).removeClass("btn-primary").addClass("btn-default");
      mqc_auto_save_config();
    }
  });

  // Use jQuery UI to make the colour filters sortable
  $("#mqc_col_filters").sortable();
  $("#mqc_col_filters").on("sortstop", function (event, ui) {
    $("#mqc_cols_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
  });

  // Apply pre-configured highlight patterns from config only if no local storage values
  let has_highlight_filters = $("#mqc_col_filters").children().length > 0;
  if (!has_highlight_filters && mqc_config.highlight_patterns && mqc_config.highlight_patterns.length > 0) {
    // Set regex mode if specified
    if (mqc_config.highlight_regex) {
      $("#mqc_cols .mqc_regex_mode input").prop("checked", true);
    }

    // Add each pattern with its color
    for (let i = 0; i < mqc_config.highlight_patterns.length; i++) {
      const pattern = mqc_config.highlight_patterns[i];
      let color =
        mqc_config.highlight_colors && mqc_config.highlight_colors[i] ? mqc_config.highlight_colors[i] : "#e41a1c";

      // Add to the filters list
      $("#mqc_col_filters").append(window.make_colorsamples_filter(pattern, color));
    }

    // Apply the highlights
    apply_mqc_highlights();
  }
};

window.mqc_colours_idx = 0;
window.make_colorsamples_filter = function (f_text, f_col) {
  let row = `
  <li style="color: ${f_col};" id="${hashCode(f_text + f_col)}" class="d-flex justify-content-between">
    <span class="hc_handle">&boxV;</span>
    <input class="f_text flex-grow-1" value="${f_text}" tabindex="${mqc_colours_idx}" />
    <button type="button" class="btn-close py-2 mt-1" aria-label="Remove"></button>
  </li>`;
  window.mqc_colours_idx += 1;
  return row;
};
