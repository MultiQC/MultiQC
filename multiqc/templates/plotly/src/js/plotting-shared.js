////////////////////////////////////////////////
// Engine-neutral plotting code (shared between Plotly and other rendering engines)
////////////////////////////////////////////////

// Global plot data variable. Accessed in many other JavaScript files.
window.mqc_plots = {};

// Initialise the toolbox filters
window.mqc_highlight_f_texts = [];
window.mqc_highlight_f_cols = [];
window.mqc_highlight_regex_mode = false;
window.mqc_rename_f_texts = [];
window.mqc_rename_t_texts = [];
window.mqc_hide_mode = "hide";
window.mqc_hide_f_texts = [];
window.mqc_hide_regex_mode = false;

function getPseudonym(sampleName) {
  // If anonymization is enabled for AI prompts, use the aiPseudonymMap built in Python runtime,
  // which defines pseudonyms for original sample names (before renaming)

  // See if toolbox switch is on
  const anonymizationEnabled = getStoredSampleAnonymizationEnabled();
  if (!anonymizationEnabled) return undefined;

  // aiPseudonymMap is built in Python runtime and defined in head.html
  if (!aiPseudonymMap) return undefined;

  // Check the map. Exact match?
  if (aiPseudonymMap[sampleName]) return aiPseudonymMap[sampleName];

  // Try replacing partial matches for cases like sample="SAMPLE1-SAMPLE2"
  // Start with the longest original name to avoid situations when one sample is a prefix of another
  let sortedOriginals = Object.keys(aiPseudonymMap).sort((a, b) => b.length - a.length);
  let result = sampleName;
  for (let original of sortedOriginals) {
    if (result.includes(original)) {
      result = result.replace(original, aiPseudonymMap[original]);
    }
  }
  return result;
}

class Sample {
  constructor(name) {
    this.originalName = name;
    this.name = name;
    this.highlight = null;
    this.hidden = false;
    this.pseudonym = getPseudonym(name);
  }
}

// Highlighting, hiding and renaming samples. Takes a list of samples, returns
// a list of objects: {"name": "new_name", "highlight": "#cccccc", "hidden": false}
function applyToolboxSettings(samples, plotAnchor) {
  // init object with default values, apply pseudonymization
  let objects = samples.map((name) => new Sample(name));

  // Rename samples
  if (window.mqc_rename_f_texts.length > 0) {
    objects.map((obj) => {
      for (let patternIdx = 0; patternIdx < window.mqc_rename_f_texts.length; patternIdx++) {
        let pattern = window.mqc_rename_f_texts[patternIdx];
        let new_text = window.mqc_rename_t_texts[patternIdx];
        obj.name = obj.name.replace(pattern, new_text);
      }
    });
  }

  // Highlight samples
  if (window.mqc_highlight_f_texts.length > 0) {
    objects.map((obj) => {
      for (let i = 0; i < window.mqc_highlight_f_texts.length; i++) {
        const f_text = window.mqc_highlight_f_texts[i];
        const f_col = window.mqc_highlight_f_cols[i];
        let match = false;
        if (window.mqc_highlight_regex_mode) {
          if (obj.name.match(f_text)) match = true;
        } else {
          if (obj.name.indexOf(f_text) > -1) match = true;
        }
        if (match) obj.highlight = f_col;
      }
    });
  }

  // Hide samples
  if (window.mqc_hide_f_texts.length > 0) {
    let groupDiv = $("#" + plotAnchor).closest(".mqc_hcplot_plotgroup");
    groupDiv.parent().find(".samples-hidden-warning").remove();
    groupDiv.show();

    objects.map((obj) => {
      let match = false;
      for (let i = 0; i < window.mqc_hide_f_texts.length; i++) {
        const f_text = window.mqc_hide_f_texts[i];
        if (window.mqc_hide_regex_mode) {
          if (obj.name.match(f_text)) match = true;
        } else {
          if (obj.name.indexOf(f_text) > -1) match = true;
        }
      }
      if (window.mqc_hide_mode === "show") match = !match;
      if (match) obj.hidden = true;
    });

    // Some series hidden. Show a warning text string.
    let nHidden = objects.filter((obj) => obj.hidden).length;
    if (nHidden > 0) {
      const alert = `
      <div class="samples-hidden-warning alert alert-warning">
        ⚠ <strong>Warning:</strong> ${nHidden} samples hidden.
        <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose('#mqc_hidesamples', true); return false;">See toolbox.</a>
      </div>`;
      groupDiv.before(alert);
    }
    // All series hidden. Hide the graph.
    if (nHidden === objects.length) {
      groupDiv.hide();
      return objects;
    }
  }

  // Return the object indexed by sample names
  return objects;
}

// Make functions available globally
window.applyToolboxSettings = applyToolboxSettings;

function updateObject(target, source, nullOnly = false) {
  // Iterate through all keys in the source object
  for (const key in source) {
    // Check if the value is not null
    // Check if the value is an object and not an array
    if (Array.isArray(source[key])) {
      // Recursively update the array
      if (!nullOnly || target[key] === undefined || target[key] === null) {
        target[key] = [...source[key]];
      }
    } else if (typeof source[key] === "object") {
      // If the target doesn't have this key, or it's not an object, initialize it
      if (!target[key] || typeof target[key] !== "object") {
        target[key] = {};
      }
      // Recursively update the object
      updateObject(target[key], source[key], nullOnly);
    } else {
      if (!nullOnly || target[key] === undefined || target[key] === null) {
        // Directly update the value
        target[key] = source[key];
      }
    }
  }
}

// Make updateObject globally available
window.updateObject = updateObject;

let loadingWarning;

$(function () {
  // Show loading warning
  loadingWarning = $(".mqc_loading_warning").show();
});

window.callAfterDecompressed.push(function (mqc_plotdata) {
  window.mqc_plots = Object.fromEntries(
    Object.values(mqc_plotdata).map((data) => [data.anchor, window.initPlot(data)]),
  );

  let shouldLoad = $(".hc-plot.not_loaded:visible");

  // Show plots on page load: either render, or show the "Show Plot" button
  shouldLoad.each(function () {
    let anchor = $(this).attr("id");
    let plot = mqc_plots[anchor];
    setTimeout(function () {
      // Deferring each plot call prevents browser from locking up
      if (plot.deferRender) {
        $("#" + anchor)
          .removeClass("not_loaded")
          .html('<button class="btn btn-outline-secondary btn-lg render_plot">Show plot</button>');
      } else {
        window.renderPlot(anchor);
      }
      if ($(".hc-plot.not_loaded:visible").length === 0)
        // All plots loaded successfully (rendered or deferred with "Show Plot"), so hiding the warning
        $(".mqc_loading_warning").hide();
    }, 50);
  });

  // All plots loaded successfully, so hiding the warning
  if (shouldLoad.length === 0) loadingWarning.hide();

  // Render a plot when clicked (heavy plots are not automatically rendered by default)
  $("body").on("click", ".render_plot", function (e) {
    let plotAnchor = $(this).parent().attr("id");
    window.renderPlot(plotAnchor);
  });

  // Button "Render all plots" clicked, so rendering everything, and hiding the parent button object
  $("#mqc-render-all-plots").click(function () {
    $(".hc-plot").each(function () {
      window.renderPlot($(this).attr("id"));
    });
    $(this).parent().hide();
  });

  // Replot graphs when something changed in filters
  $(document).on("mqc_highlights mqc_renamesamples mqc_hidesamples", function () {
    // Replot graphs
    $(".hc-plot:not(.not_rendered)").each(function () {
      window.renderPlot($(this).attr("id"));
    });
  });

  // A "Percentages" button above a plot is clicked
  $("button.interactive-switch-group.percent-switch").click(function (e) {
    e.preventDefault();
    let plotAnchor = $(this).data("plot-anchor");

    // Toggling flags
    mqc_plots[plotAnchor].pActive = !$(this).hasClass("active");
    $(this).toggleClass("active");

    if (mqc_plots[plotAnchor].rendered) {
      window.renderPlot(plotAnchor); // re-render
    }
  });

  // A "Log" button above a plot is clicked
  $("button.interactive-switch-group.log10-switch").click(function (e) {
    e.preventDefault();
    let plotAnchor = $(this).data("plot-anchor");

    // Toggling flags
    mqc_plots[plotAnchor].lActive = !$(this).hasClass("active");
    $(this).toggleClass("active");

    if (mqc_plots[plotAnchor].rendered) {
      window.renderPlot(plotAnchor); // re-render
    }
  });

  // Switch data source
  $(".interactive-switch-group.dataset-switch-group button").click(function (e) {
    e.preventDefault();
    let el = $(this);
    if (el.hasClass("active")) return;
    el.siblings("button.active").removeClass("active");
    el.addClass("active");
    let plotAnchor = el.data("plot-anchor");
    let activeDatasetIdx = mqc_plots[plotAnchor].activeDatasetIdx;
    let newDatasetIdx = el.data("datasetIndex");
    mqc_plots[plotAnchor].activeDatasetIdx = newDatasetIdx;
    if (activeDatasetIdx === newDatasetIdx) return;

    if (mqc_plots[plotAnchor].rendered) {
      window.renderPlot(plotAnchor); // re-render
    }
  });

  // Make divs height-draggable
  // http://jsfiddle.net/Lkwb86c8/
  $(".hc-plot:not(.no-handle)").each(function () {
    let el = $(this);
    if (!el.parent().hasClass("hc-plot-wrapper")) {
      el.wrap('<div class="hc-plot-wrapper"></div>');
    }
    if (!el.siblings().hasClass("hc-plot-handle")) {
      el.after('<div class="hc-plot-handle"><span></span><span></span><span></span></div>');
    }
    el.css({ height: "auto", top: 0, bottom: "6px", position: "absolute" });
  });

  $(".hc-plot-handle").on("mousedown", function (e) {
    let wrapper = $(this).parent();
    let plotAnchor = wrapper.children(".hc-plot")[0].id;
    let startHeight = wrapper.height();
    let pY = e.pageY;

    let doc = $(document);
    doc.on("mouseup", function () {
      // Clear listeners now that we've let go
      doc.off("mousemove");
      doc.off("mouseup");
      // Fire off a custom jQuery event for other javascript chunks to tie into
      // Bind to the plot div, which should have a custom ID
      $(wrapper.parent().find(".hc-plot, .beeswarm-plot")).trigger("mqc_plotresize");
    });

    $(document).on("mousemove", function (me) {
      let newHeight = startHeight + (me.pageY - pY) + 2; // 2 px for the border or something
      wrapper.css("height", newHeight);
      if (mqc_plots[plotAnchor] !== undefined) mqc_plots[plotAnchor].resize(newHeight - 7); // 7 is the height of the handle overlapping the plot wrapper
    });
  });
});
