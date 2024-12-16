////////////////////////////////////////////////
// Base JS for MultiQC Reports
////////////////////////////////////////////////

// Collect functions to be called after plot data is decompressed.
// Includes functions in plotting.js, and any module-specific JS like multiqc_fastqc.js
let callAfterDecompressed = [];

// Helper config - is defined and object length > 0?
function notEmptyObj(obj) {
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
}

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
