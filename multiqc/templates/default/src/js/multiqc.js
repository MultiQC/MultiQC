////////////////////////////////////////////////
// Base JS for MultiQC Reports
////////////////////////////////////////////////

// Collect functions to be called after plot data is decompressed.
// Includes functions in plotting.js, and any module-specific JS like multiqc_fastqc.js
window.callAfterDecompressed = [];

// Helper config - is defined and object length > 0?
window.notEmptyObj = function (obj) {
  try {
    if (obj === undefined) {
      return false;
    }
    if (obj.length === 0) {
      return false;
    }
  } catch (e) {
    return false;
  }
  return true;
};

// Bootstrap toast utility function
window.showToast = function (heading, text, icon = "info", hideAfter = 5000) {
  const template = document.getElementById("toast-template");
  const container = document.getElementById("toast-container");

  // Clone the template
  const toastElement = template.cloneNode(true);

  // Create unique ID for this toast
  const toastId = "toast-" + Date.now() + "-" + Math.random().toString(36).substr(2, 9);
  toastElement.id = toastId;

  // Remove hidden class and show the toast
  toastElement.classList.remove("d-none");

  // Map icon types to Bootstrap classes and colors
  const iconMap = {
    success: { bg: "text-bg-success", icon: "✓" },
    warning: { bg: "text-bg-warning", icon: "⚠" },
    error: { bg: "text-bg-danger", icon: "✗" },
    info: { bg: "text-bg-info", icon: "ℹ" },
  };

  const iconConfig = iconMap[icon] || iconMap.info;

  // Update the toast content
  const header = toastElement.querySelector(".toast-header");
  const iconSpan = toastElement.querySelector(".toast-icon");
  const headingElement = toastElement.querySelector(".toast-heading");
  const bodyElement = toastElement.querySelector(".toast-body");

  header.className = `toast-header ${iconConfig.bg}`;
  iconSpan.textContent = iconConfig.icon;
  headingElement.textContent = heading;
  bodyElement.innerHTML = text;

  // Add toast to container
  container.appendChild(toastElement);

  // Initialize and show the toast
  const bsToast = new bootstrap.Toast(toastElement, {
    autohide: hideAfter !== false,
    delay: hideAfter,
  });

  // Show the toast
  bsToast.show();

  // Remove from DOM when hidden to prevent accumulation
  toastElement.addEventListener("hidden.bs.toast", () => {
    toastElement.remove();
  });
};

$(function () {
  // Enable the bootstrap tooltip hovers
  $('[data-toggle="tooltip"]').tooltip();

  // Side nav expansion
  $("#side-nav-handle").click(function (e) {
    $(".mainpage, .side-nav, .footer").toggleClass("hidden-nav");
    $("#side-nav-handle span").toggleClass("glyphicon-triangle-left glyphicon-triangle-right");
    // send resize trigger for replotting after css animation
    setTimeout(function () {
      $(document).resize();
    }, 510);
  });

  // Hide welcome alert if setting saved
  try {
    let hide_welcome = localStorage.getItem("mqc_hide_welcome");
    if (hide_welcome !== "true") {
      $("#mqc_header_hr").show();
      $("#mqc_welcome").show();
    }
    $("#mqc_hide_welcome_btn").click(function (e) {
      localStorage.setItem("mqc_hide_welcome", "true");
    });
  } catch (e) {
    console.log(
      "Could not access localStorage: " +
        e +
        "\nPlease disable 'Block third-party cookies and site data' or browser equivalent.",
    );
  }
  $("#mqc_hide_welcome_btn, #mqc_welcome .close").click(function (e) {
    $("#mqc_header_hr").show();
  });

  // Initialise module DOI popovers
  $(".module-doi").click(function (e) {
    // Don't follow the link
    e.preventDefault();
    let el = $(this);

    // Check if we already have a popover
    if (el.data("bs.popover")) {
      return;
    }

    // Check if we already failed to get a popover
    if (el.hasClass("no-doi-details")) {
      window.open(el.attr("href"), "_blank");
    }

    // Get full paper details
    var doi = $(this).data("doi");
    $.get("https://api.crossref.org/works/" + doi, function (data) {
      // Prepare fields
      var title = data.message.title[0];
      var journal = data.message["short-container-title"][0];
      var year = data.message.published["date-parts"][0][0];
      var authors = [];
      $.each(data.message.author, function (idx, author) {
        authors.push(author.given);
      });
      var content = "<p><em>" + journal + "</em> (" + year + ")</p>";
      content += "<p class='small'>" + authors.join(", ") + "</p>";
      content += '<p><a href="https://doi.org/' + doi + '" class="btn btn-primary" target="_blank">View paper</a></p>';
      // Make the popover
      el.popover({
        title: title,
        content: content,
        html: true,
        trigger: "focus",
        placement: "auto left",
      });
      el.popover("show");
    }).fail(function () {
      el.addClass("no-doi-details");
      window.open(el.attr("href"), "_blank");
    });
  });
});
