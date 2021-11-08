////////////////////////////////////////////////
// Base JS for MultiQC Reports
////////////////////////////////////////////////

// Helper config - is defined and object length > 0?
function notEmptyObj(obj) {
  try {
    if (obj === undefined) {
      return false;
    }
    if (obj.length == 0) {
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
    var hide_welcome = localStorage.getItem("mqc_hide_welcome");
    if (hide_welcome !== "true") {
      $("#mqc_header_hr").slideUp();
      $("#mqc_welcome").slideDown();
    }
    $("#mqc_hide_welcome_btn").click(function (e) {
      localStorage.setItem("mqc_hide_welcome", "true");
    });
  } catch (e) {
    console.log(
      "Could not access localStorage: " +
        e +
        "\nPlease disable 'Block third-party cookies and site data' or browser equivalent."
    );
  }
  $("#mqc_hide_welcome_btn, #mqc_welcome .close").click(function (e) {
    $("#mqc_header_hr").show();
  });

  // Initialise module DOI popovers
  $(".module-doi").each(function () {
    var el = $(this);
    // Get full paper details
    var doi = $(this).data("doi");
    $.get("https://api.crossref.org/works/" + doi, function (data) {
      console.log(data);
      // Authors
      var authors = "";
      $.each(data.message.author, function (idx, author) {
        authors += author.given + ", ";
      });
      // Make the popover
      el.popover({
        title: data.message.title[0],
        content:
          "<p><em>" +
          data.message["short-container-title"][0] +
          "</em> (" +
          data.message.published["date-parts"][0][0] +
          ")</p>" +
          "<p class='small'>" +
          authors +
          '</p><p><a href="https://doi.org/' +
          doi +
          '" class="btn btn-primary" target="_blank">View paper</a></p>',
        html: true,
        trigger: "focus",
      });
    });
  });
  // Stop the DOI link from working, as we have popovers
  $(".module-doi").click(function (e) {
    e.preventDefault();
  });
});
