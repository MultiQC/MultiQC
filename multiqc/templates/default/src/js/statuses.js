////////////////////////////////////////////////
// MultiQC Module Status Bars
// Generic status bar implementation for any module
////////////////////////////////////////////////

// Global vars for status data
let mqc_status_data = {}; // { <module_key>: { <section_key>: { <sample_name>: "pass"|"warn"|"fail" } } }

// Load status data from embedded JSON
function load_mqc_status_data() {
  $(".mqc-status-data").each(function (i, elem) {
    const data = JSON.parse(elem.innerHTML);
    const [module_key, section_key, sample_statuses] = data;

    if (!mqc_status_data[module_key]) {
      mqc_status_data[module_key] = {};
    }
    mqc_status_data[module_key][section_key] = sample_statuses;
  });
}

// Initialize status bars on page load
$(function () {
  callAfterDecompressed.push(function (mqc_plotdata) {
    load_mqc_status_data();

    // Initialize all status bars
    $(".mqc-status-progress-wrapper").each(function () {
      const $wrapper = $(this);
      const module_key = $wrapper.data("module-key");
      const section_key = $wrapper.data("section-key");

      init_status_bar($wrapper, module_key, section_key);
    });
  });
});

// Initialize a single status bar with popovers
function init_status_bar($wrapper, module_key, section_key) {
  const vals = mqc_status_data[module_key]?.[section_key];
  if (!vals) return;

  let statusHoverTimeout = null;
  let statusHideTimeout = null;

  $wrapper.find(".mqc-status-progress > .progress").each(function () {
    const element = this;
    const $progressBar = $(this).find(".progress-bar");
    const passes = $progressBar.hasClass("bg-success");
    const warns = $progressBar.hasClass("bg-warning");
    const fails = $progressBar.hasClass("bg-danger");

    let pclass = "";
    if (passes) pclass = "success";
    if (warns) pclass = "warning";
    if (fails) pclass = "danger";

    // Collect sample names for this status
    let samples = [];
    $.each(vals, function (s_name, status) {
      if ((status === "pass" && passes) || (status === "warn" && warns) || (status === "fail" && fails)) {
        samples.push(s_name);
      }
    });

    // Build popover content
    const popoverContent = `
      <div class="mqc-status-popover-intro mb-2">
        <p class="mb-1 text-center">Click bar to fix in place</p>
        <div class="row">
          <div class="col">
            <button class="mqc-status-highlight btn btn-sm btn-outline-secondary w-100">
              Highlight
            </button>
          </div>
          <div class="col">
            <button class="mqc-status-hideothers btn btn-sm btn-outline-secondary w-100">
              Filter
            </button>
          </div>
        </div>
      </div>
      <ul class="list-unstyled mb-0">${samples
        .sort()
        .map((s) => `<li>${s}</li>`)
        .join("")}</ul>`;

    // Create popover with Bootstrap 5
    const popover = new bootstrap.Popover(element, {
      title: $(this).attr("aria-label"),
      content: popoverContent,
      html: true,
      sanitize: false,
      trigger: "manual",
      placement: "bottom",
      customClass: `popover-mqc-status popover-${pclass}`,
    });

    // Track if popover is pinned
    $(element).data("popover-pinned", false);

    // Show on mouseover (if not pinned)
    $(element).on("mouseenter", function () {
      const el = this;
      if (statusHideTimeout) {
        clearTimeout(statusHideTimeout);
        statusHideTimeout = null;
      }
      if (!$(el).data("popover-pinned")) {
        popover.show();
        // Attach hover handlers to popover
        setTimeout(function () {
          const popoverEl = document.querySelector(".popover-mqc-status");
          if (popoverEl) {
            $(popoverEl)
              .off("mouseenter mouseleave")
              .on("mouseenter", function () {
                if (statusHideTimeout) {
                  clearTimeout(statusHideTimeout);
                  statusHideTimeout = null;
                }
              })
              .on("mouseleave", function () {
                if (!$(el).data("popover-pinned")) {
                  statusHideTimeout = setTimeout(function () {
                    if (!$(el).data("popover-pinned")) {
                      popover.hide();
                    }
                  }, 200);
                }
              });
          }
        }, 50);
      }
    });

    // Hide on mouseout (if not pinned)
    $(element).on("mouseleave", function () {
      const el = this;
      if (!$(el).data("popover-pinned")) {
        // Check if cursor moved to the popover
        setTimeout(function () {
          if (!$(".popover-mqc-status:hover").length) {
            statusHideTimeout = setTimeout(function () {
              if (!$(el).data("popover-pinned")) {
                popover.hide();
              }
            }, 200);
          }
        }, 50);
      }
    });

    // Pin on click
    $(element).on("click", function (e) {
      e.stopPropagation();
      if (statusHideTimeout) {
        clearTimeout(statusHideTimeout);
        statusHideTimeout = null;
      }
      $(this).data("popover-pinned", true);
      popover.show();
    });
  });

  // Close pinned popovers when clicking outside
  $(document).on("click", function (e) {
    if (!$(e.target).closest(".mqc-status-progress .progress, .popover").length) {
      $wrapper.find(".mqc-status-progress > .progress").each(function () {
        const popover = bootstrap.Popover.getInstance(this);
        if (popover) {
          popover.hide();
          $(this).data("popover-pinned", false);
        }
      });
    }
  });
}

// Global event handlers for highlight/filter buttons
$(document).on("click", ".mqc-status-highlight", function (e) {
  e.preventDefault();
  // Get sample names from list items
  let samples = $(this)
    .closest(".popover-body")
    .find("ul li")
    .map(function () {
      return $(this).text();
    })
    .get();

  // Hide the popover first
  const popoverEl = $(this).closest(".popover-mqc-status");
  if (popoverEl.length) {
    $(".mqc-status-progress > .progress").each(function () {
      const popoverInstance = bootstrap.Popover.getInstance(this);
      if (popoverInstance) {
        popoverInstance.hide();
        $(this).data("popover-pinned", false);
      }
    });
  }

  let f_col = $("#mqc_colour_filter_color").val();
  // Add sample names to the toolbox
  for (let i = 0; i < samples.length; i++) {
    let f_text = samples[i];
    $("#mqc_col_filters").append(window.make_colorsamples_filter(f_text, f_col));
  }
  // Apply highlights and open toolbox
  apply_mqc_highlights();
  mqc_toolbox_openclose("#mqc_cols", true);
  // Increment color index and update next highlight colour
  window.mqc_colours_idx += 1;
  $("#mqc_colour_filter_color").val(mqc_colours[window.mqc_colours_idx % mqc_colours.length]);
});

$(document).on("click", ".mqc-status-hideothers", function (e) {
  e.preventDefault();
  // Get sample names from list items
  let samples = $(this)
    .closest(".popover-body")
    .find("ul li")
    .map(function () {
      return $(this).text();
    })
    .get();

  // Hide the popover first
  const popoverEl = $(this).closest(".popover-mqc-status");
  if (popoverEl.length) {
    $(".mqc-status-progress > .progress").each(function () {
      const popoverInstance = bootstrap.Popover.getInstance(this);
      if (popoverInstance) {
        popoverInstance.hide();
        $(this).data("popover-pinned", false);
      }
    });
  }

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
  $("#mqc_hidesamples .mqc_regex_mode input").prop("checked", false);
  // Add sample names to the toolbox
  for (let i = 0; i < samples.length; i++) {
    let f_text = samples[i];
    $("#mqc_hidesamples_filters").append(window.make_hidesamples_filter(f_text));
  }
  // Apply highlights and open toolbox
  apply_mqc_hidesamples();
  mqc_toolbox_openclose("#mqc_hidesamples", true);
});
