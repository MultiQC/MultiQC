////////////////////////////////////////////////
// MultiQC Report Toolbox - Rename Samples
////////////////////////////////////////////////

// Make functions available globally
window.apply_mqc_renamesamples = function () {
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
};

window.initRename = function () {
  // Initialize rename functionality

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

  // Apply pre-configured sample renaming patterns from config only if no local storage values
  let has_rename_filters = $("#mqc_renamesamples_filters").children().length > 0;
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
};
