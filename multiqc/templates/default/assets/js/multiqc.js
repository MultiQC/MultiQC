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

$(document).ready(function () {
  $("#ai-continue-in-chat").click(function (e) {
    e.preventDefault();
    let el = $(this);
    let website = el.data("website");
    let encodedSystemMessage = el.data("encoded-system-message");
    let encodedChatMessages = el.data("encoded-chat-messages");
    let reportUuid = el.data("report-uuid");

    let url = website + "/ask-ai/";
    if (reportUuid) {
      url += "?multiqc-report-uuid=" + reportUuid;
    }

    const chatWindow = window.open(url, "_blank");

    if (encodedSystemMessage && encodedChatMessages) {
      function sendMessage() {
        chatWindow.postMessage(
          {
            type: "chatInitialMessages",
            content: {
              encodedSystemMessage: encodedSystemMessage,
              encodedChatMessages: encodedChatMessages,
            },
          },
          website,
        );
      }
      setTimeout(sendMessage, 2000);
    }
  });

  // Add "Show More" button to AI summary
  $(".ai-summary").each(function () {
    const $details = $(this).find("details");
    const $showMoreBtn = $(this).find(".ai-summary-expand");
    $showMoreBtn.on("click", function (e) {
      if ($details.prop("open")) {
        $details.prop("open", false);
      } else {
        $details.prop("open", true);
      }
    });

    // Do no expand when clicked on the whole area
    $details.on("click", function (e) {
      e.preventDefault();
    });
  });

  $(".ai-summary sample").hover(function () {
    $(this).css("opacity", 0.9);
  });

  $(".btn-ai-generate-summary").click(async function (e) {
    e.preventDefault();
    let el = $(this);

    // Disable button and show spinner
    el.prop("disabled", true);
    el.html("Generating...");

    try {
      let moduleAnchor = el.data("module-anchor");
      let sectionAnchor = el.data("section-anchor");
      let plotAnchor = el.data("plot-anchor");
      let url = el.data("seqera-ai-url");
      let plot = mqc_plots[plotAnchor];
      let data = plot.prepData();

      let moduleWrapper = $("#mqc-module-section-" + moduleAnchor);
      let moduleName = moduleWrapper.find(".mqc-module-title").text();
      let moduleSummary = moduleWrapper.find(".mqc-module-section-first > p").text();
      let moduleComment = moduleWrapper.find(".mqc-section-comment").text();

      let sectionName = el.data("section-name");
      let sectionWrapper = $("#mqc-section-wrapper-" + sectionAnchor);
      let description = sectionWrapper.find(".mqc-section-description").text();
      let comment = sectionWrapper.find(".mqc-section-comment").text();
      let helptext = sectionWrapper.find(".mqc-section-helptext").text();

      let prompt = `
        You are an expert in bioinformatics, sequencing technologies, and genomics data analysis.
        You are given plot data a MultiQC report that was generated by a bioinformatics workflow.
        Your task is to analyse the data, and give a very short and concise overall summary for the results.
        Don't waste words: mention only the important QC issues. If there are no issues, just say so.
        Limit it to 1-2 sentences.
      `;

      let report = `
        QC tool that generated data for the plot: ${moduleName} (${moduleSummary}). ${
          moduleComment ? `\nComment: ${moduleComment}` : ""
        }

        Section: ${sectionName} ${description ? `(${description})` : ""}. ${comment ? `\nComment: ${comment}` : ""} ${
          helptext ? `\nHelptext: ${helptext}` : ""
        }

        Plot type: ${plot.plotType}. ${plot.plotDescription ? `\nPlot description: ${plot.plotDescription}` : ""} ${
          plot.plotHelptext ? `\nPlot helptext: ${plot.plotHelptext}` : ""
        }

        Data: ${JSON.stringify(data)}
      `;

      // Request the AI summary from http://127.0.0.1:8000//interpret-multiqc-report
      const response = await fetch(`${url}/interpret-multiqc-report`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          system_message: prompt,
          report: report,
        }),
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const result = await response.json();
      console.log(result);

      let interp = result.interpretation;

      const betaIcon = `
        <span style="vertical-align: middle">
          <svg width="30" height="12" viewBox="0 0 49 19" fill="none" xmlns="http://www.w3.org/2000/svg">
              <rect x="1.4375" y="0.5" width="47" height="18" rx="9" fill="#160F26" fill-opacity="0.1"/>
              <rect x="1.4375" y="0.5" width="47" height="18" rx="9" stroke="#160F26"/>
              <path d="M13.0392 14V5.27273H16.0904C16.6983 5.27273 17.1998 5.37784 17.5946 5.58807C17.9895 5.79545 18.2836 6.07528 18.4767 6.42756C18.6699 6.77699 18.7665 7.16477 18.7665 7.59091C18.7665 7.96591 18.6998 8.27557 18.5662 8.51989C18.4355 8.7642 18.2623 8.95739 18.0463 9.09943C17.8333 9.24148 17.6017 9.34659 17.3517 9.41477V9.5C17.6188 9.51705 17.8873 9.6108 18.1571 9.78125C18.427 9.9517 18.6529 10.196 18.8347 10.5142C19.0165 10.8324 19.1074 11.2216 19.1074 11.6818C19.1074 12.1193 19.008 12.5128 18.8091 12.8622C18.6103 13.2116 18.2963 13.4886 17.8674 13.6932C17.4384 13.8977 16.8801 14 16.1926 14H13.0392ZM14.0961 13.0625H16.1926C16.883 13.0625 17.373 12.929 17.6628 12.6619C17.9554 12.392 18.1017 12.0653 18.1017 11.6818C18.1017 11.3864 18.0265 11.1136 17.8759 10.8636C17.7253 10.6108 17.5108 10.4091 17.2324 10.2585C16.954 10.1051 16.6245 10.0284 16.2438 10.0284H14.0961V13.0625ZM14.0961 9.10795H16.0563C16.3745 9.10795 16.6614 9.04545 16.9171 8.92045C17.1756 8.79545 17.3801 8.61932 17.5307 8.39205C17.6841 8.16477 17.7608 7.89773 17.7608 7.59091C17.7608 7.20739 17.6273 6.8821 17.3603 6.61506C17.0932 6.34517 16.6699 6.21023 16.0904 6.21023H14.0961V9.10795ZM23.4519 14.1364C22.8212 14.1364 22.2772 13.9972 21.8198 13.7188C21.3652 13.4375 21.0144 13.0455 20.7672 12.5426C20.5229 12.0369 20.4007 11.4489 20.4007 10.7784C20.4007 10.108 20.5229 9.51705 20.7672 9.00568C21.0144 8.49148 21.3581 8.09091 21.7985 7.80398C22.2417 7.5142 22.7587 7.36932 23.3496 7.36932C23.6905 7.36932 24.0272 7.42614 24.3596 7.53977C24.6919 7.65341 24.9945 7.83807 25.2672 8.09375C25.54 8.34659 25.7573 8.68182 25.9192 9.09943C26.0811 9.51705 26.1621 10.0312 26.1621 10.642V11.0682H21.1167V10.1989H25.1394C25.1394 9.82955 25.0655 9.5 24.9178 9.21023C24.7729 8.92045 24.5655 8.69176 24.2956 8.52415C24.0286 8.35653 23.7132 8.27273 23.3496 8.27273C22.949 8.27273 22.6025 8.37216 22.3098 8.57102C22.0201 8.76705 21.7971 9.02273 21.6408 9.33807C21.4846 9.65341 21.4064 9.99148 21.4064 10.3523V10.9318C21.4064 11.4261 21.4917 11.8452 21.6621 12.1889C21.8354 12.5298 22.0755 12.7898 22.3823 12.9688C22.6891 13.1449 23.0456 13.233 23.4519 13.233C23.7161 13.233 23.9547 13.196 24.1678 13.1222C24.3837 13.0455 24.5698 12.9318 24.726 12.7812C24.8823 12.6278 25.003 12.4375 25.0882 12.2102L26.0598 12.483C25.9576 12.8125 25.7857 13.1023 25.5442 13.3523C25.3027 13.5994 25.0044 13.7926 24.6493 13.9318C24.2942 14.0682 23.8951 14.1364 23.4519 14.1364ZM30.5385 7.45455V8.30682H27.1465V7.45455H30.5385ZM28.1351 5.88636H29.1408V12.125C29.1408 12.4091 29.182 12.6222 29.2644 12.7642C29.3496 12.9034 29.4576 12.9972 29.5882 13.0455C29.7218 13.0909 29.8624 13.1136 30.0101 13.1136C30.1209 13.1136 30.2118 13.108 30.2828 13.0966C30.3539 13.0824 30.4107 13.071 30.4533 13.0625L30.6578 13.9659C30.5897 13.9915 30.4945 14.017 30.3723 14.0426C30.2502 14.071 30.0953 14.0852 29.9078 14.0852C29.6238 14.0852 29.3453 14.0241 29.0726 13.902C28.8027 13.7798 28.5783 13.5938 28.3993 13.3438C28.2232 13.0938 28.1351 12.7784 28.1351 12.3977V5.88636ZM33.9775 14.1534C33.5627 14.1534 33.1863 14.0753 32.8482 13.919C32.5101 13.7599 32.2417 13.5312 32.0428 13.233C31.8439 12.9318 31.7445 12.5682 31.7445 12.142C31.7445 11.767 31.8184 11.4631 31.9661 11.2301C32.1138 10.9943 32.3113 10.8097 32.5584 10.6761C32.8056 10.5426 33.0783 10.4432 33.3766 10.3778C33.6777 10.3097 33.9803 10.2557 34.2843 10.2159C34.682 10.1648 35.0044 10.1264 35.2516 10.1009C35.5016 10.0724 35.6834 10.0256 35.7971 9.96023C35.9135 9.89489 35.9718 9.78125 35.9718 9.61932V9.58523C35.9718 9.16477 35.8567 8.83807 35.6266 8.60511C35.3993 8.37216 35.0542 8.25568 34.5911 8.25568C34.111 8.25568 33.7346 8.3608 33.4618 8.57102C33.1891 8.78125 32.9973 9.00568 32.8865 9.24432L31.932 8.90341C32.1025 8.50568 32.3297 8.19602 32.6138 7.97443C32.9007 7.75 33.2132 7.59375 33.5513 7.50568C33.8922 7.41477 34.2275 7.36932 34.557 7.36932C34.7672 7.36932 35.0087 7.39489 35.2814 7.44602C35.557 7.49432 35.8226 7.59517 36.0783 7.74858C36.3368 7.90199 36.5513 8.13352 36.7218 8.44318C36.8922 8.75284 36.9775 9.16761 36.9775 9.6875V14H35.9718V13.1136H35.9206C35.8525 13.2557 35.7388 13.4077 35.5797 13.5696C35.4206 13.7315 35.209 13.8693 34.9448 13.983C34.6806 14.0966 34.3581 14.1534 33.9775 14.1534ZM34.1309 13.25C34.5286 13.25 34.8638 13.1719 35.1365 13.0156C35.4121 12.8594 35.6195 12.6577 35.7587 12.4105C35.9007 12.1634 35.9718 11.9034 35.9718 11.6307V10.7102C35.9292 10.7614 35.8354 10.8082 35.6905 10.8509C35.5485 10.8906 35.3837 10.9261 35.1962 10.9574C35.0115 10.9858 34.8311 11.0114 34.655 11.0341C34.4817 11.054 34.3411 11.071 34.2331 11.0852C33.9718 11.1193 33.7275 11.1747 33.5002 11.2514C33.2757 11.3253 33.0939 11.4375 32.9547 11.5881C32.8184 11.7358 32.7502 11.9375 32.7502 12.1932C32.7502 12.5426 32.8794 12.8068 33.138 12.9858C33.3993 13.1619 33.7303 13.25 34.1309 13.25Z" fill="#160F26"/>
            </svg>
        </span>`;

      const sparkleIcon = `<span style="vertical-align: middle">
            <svg width="12" height="12" viewBox="0 0 16 16" fill="none" xmlns="http://www.w3.org/2000/svg">
                <path d="M6.4375 7L7.9375 1.5L9.4375 7L14.9375 8.5L9.4375 10.5L7.9375 15.5L6.4375 10.5L0.9375 8.5L6.4375 7Z" stroke="black" stroke-width="0.75" stroke-linejoin="round"/>
                <path d="M13.1786 2.82143L13.5 4L13.8214 2.82143L15 2.5L13.8214 2.07143L13.5 1L13.1786 2.07143L12 2.5L13.1786 2.82143Z" stroke="#160F26" stroke-width="0.5" stroke-linejoin="round"/>
            </svg>
        </span>`;

      let summary = `<div style="display: flex; justify-content: space-between; align-items: center">
            <b>Section AI Summary ${betaIcon}</b>
            <button class='btn btn-default btn-sm' id='ai-continue-in-chat'>Continue with ${sparkleIcon} <strong>Seqera AI</strong></button>
        </div>`;

      summary += interp.summary;
      $("#" + sectionAnchor + "_ai_summary summary").html(summary);
      $("#" + sectionAnchor + "_ai_summary .ai-summary-extra-content").html(`
          ${interp.detailed_summary ? `<p>${interp.detailed_summary}</p>` : ""}
          ${interp.recommendations ? `<p>${interp.recommendations}</p>` : ""}
        `);
      $("#" + sectionAnchor + "_ai_summary").show();

      // Hide the button after successful generation
      el.hide();
    } catch (error) {
      console.error("Error:", error);
      // Restore button on error
      el.prop("disabled", false);
      el.html(
        '<span style="vertical-align: baseline"><svg width="10" height="10" viewBox="0 0 16 14" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M6.4375 7L7.9375 1.5L9.4375 7L14.9375 8.5L9.4375 10.5L7.9375 15.5L6.4375 10.5L0.9375 8.5L6.4375 7Z" stroke="black" stroke-width="0.75" stroke-linejoin="round"/><path d="M13.1786 2.82143L13.5 4L13.8214 2.82143L15 2.5L13.8214 2.07143L13.5 1L13.1786 2.07143L12 2.5L13.1786 2.82143Z" stroke="#160F26" stroke-width="0.5" stroke-linejoin="round"/></svg></span> AI summary',
      );
    }
  });
});
