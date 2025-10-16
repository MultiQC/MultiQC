////////////////////////////////////////////////
// MultiQC Report Toolbox - Hide Samples
////////////////////////////////////////////////

// Make functions available globally
window.apply_mqc_hidesamples = function (mode) {
  // Collect the filters into an array
  if (mode === undefined) {
    mode = $(".mqc_hidesamples_showhide:checked").val() === "show" ? "show" : "hide";
  }
  let regex_mode = $("#mqc_hidesamples .mqc_regex_mode input").prop("checked");
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
};

window.initHideSamples = function () {
  // Initialize hide samples functionality

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
    $("#mqc_hidesamples_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
    $("#mqc_hidesamples_filter").val("");
    mqc_hidesamples_idx += 1;
  });

  $(".mqc_hidesamples_showhide").change(function (e) {
    $("#mqc_hidesamples_apply").attr("disabled", false).removeClass("btn-default").addClass("btn-primary");
  });

  $("#mqc_hidesamples_apply").click(function (e) {
    if (apply_mqc_hidesamples()) {
      $(this).attr("disabled", true).removeClass("btn-primary").addClass("btn-default");
      mqc_auto_save_config();
    }
  });

  // Apply pre-configured hide samples from config only if no local storage values
  let has_hide_filters = $("#mqc_hidesamples_filters").children().length > 0;
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
        $("#mqc_hidesamples .mqc_regex_mode input").prop("checked", true);
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
};
