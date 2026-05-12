////////////////////////////////////////////////
// MultiQC Report Toolbox - Save/Load Configurations
////////////////////////////////////////////////

// Add these helper functions
window.getConfigObject = function () {
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
};

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
        .prepend(`<option>${name + as_default ? " [default]" : ""}</option>`)
        .val(name + (as_default ? " [default]" : ""));
      // Success message
      $('<p class="bg-success-subtle text-success-emphasis p-2" id="mqc-save-success">Settings saved.</p>')
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
    $('<p class="bg-danger-subtle text-danger-emphasis p-2" id="mqc-cleared-success">Unset default.</p>')
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
        $("#mqc_loadconfig_form select").append(`<option>${name}</option>`).val(name);
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
      $("#mqc_renamesamples .mqc_regex_mode input").prop("checked", true);
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
      window.mqc_rename_f_texts.push(from_text);
      window.mqc_rename_t_texts.push(to_text);
      $("#mqc_renamesamples_filters").append(window.make_renamesamples_filter(from_text, to_text));
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
      $("#mqc_cols .mqc_regex_mode input").prop("checked", true);
      window.mqc_highlight_regex_mode = true;
    }
  }
  if (notEmptyObj(config["highlights_f_texts"]) && notEmptyObj(config["highlights_f_cols"])) {
    window.mqc_highlight_f_texts = [];
    window.mqc_highlight_f_cols = [];
    $.each(config["highlights_f_texts"], function (idx, f_text) {
      if (f_text.length === 0) {
        return true;
      }
      var f_col = config["highlights_f_cols"][idx];
      $("#" + hashCode(f_text + f_col)).remove();
      $("#mqc_col_filters").append(window.make_colorsamples_filter(f_text, f_col));
      window.mqc_highlight_f_texts.push(f_text);
      window.mqc_highlight_f_cols.push(f_col);
    });
    $("#mqc_colour_filter_color").val(mqc_colours[window.mqc_colours_idx % mqc_colours.length]);
    $(document).trigger("mqc_highlights", [
      window.mqc_highlight_f_texts,
      window.mqc_highlight_f_cols,
      config["highlight_regex"],
    ]);
  }

  // Apply config - hide samples
  if (notEmptyObj(config["hidesamples_regex"])) {
    if (config["hidesamples_regex"] == true) {
      $("#mqc_hidesamples .mqc_regex_mode input").prop("checked", true);
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
      if (f_text.length === 0) {
        return true;
      }
      $("#mqc_hidesamples_filters").append(window.make_hidesamples_filter(f_text));
      window.mqc_hide_f_texts.push(f_text);
    });
    $(document).trigger("mqc_hidesamples", [window.mqc_hide_f_texts, config["hidesamples_regex"]]);
  }

  // Trigger loaded event to initialise plots
  $(document).trigger("mqc_config_loaded");
  return true;
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
      showToast("Configuration Loaded", "Settings have been loaded from file successfully", "success");
    } catch (err) {
      showToast("Error Loading Configuration", "Could not parse the configuration file: " + err.message, "error");
    }
  };
  reader.readAsText(file);
}

//////////////////////////////////////////////////////
// SAVE SETTINGS AUTOMATICALLY
//////////////////////////////////////////////////////
window.mqc_auto_save_config = function () {
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
};

// Make functions available globally
window.initSaveLoad = function () {
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
      if (load_mqc_config(name)) {
        // Show success message
        showToast("Configuration Loaded", "Settings have been loaded successfully", "success");
      } else {
        showToast("Error Loading Configuration", "Could not load the configuration", "error");
      }
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

  // Lazy load file input handler
  let fileInputInitialized = false;
  $("#mqc_saveconfig").on("mouseenter", function () {
    if (!fileInputInitialized) {
      // Initialize file input handler only when user interacts with the save section
      $("#mqc_load_config_file_wrapper").append(
        '<input type="file" class="form-control" id="mqc_load_config_file" accept=".json">',
      );
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
};
