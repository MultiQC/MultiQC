// Initialise module DOI popovers
$(function () {
  const HOVER_DELAY = 500;
  const HIDE_DELAY = 1000;
  const POPOVER_MOUNT_DELAY = 100;
  const MOUSELEAVE_CHECK_DELAY = 50;

  let doiHoverTimeout = null;
  let doiHideTimeout = null;

  // Helper functions
  const clearTimer = (timer) => timer && clearTimeout(timer);
  const isClickedOpen = (el) => $(el).data("popover-clicked");
  const setClickedOpen = (el, value) => $(el).data("popover-clicked", value);

  function scheduleHide(el) {
    clearTimer(doiHideTimeout);
    doiHideTimeout = setTimeout(() => {
      const popover = bootstrap.Popover.getInstance(el);
      if (popover && !isClickedOpen(el)) {
        popover.hide();
      }
    }, HIDE_DELAY);
  }

  function attachPopoverHoverHandlers(el) {
    const popoverEl = document.querySelector(".popover");
    if (!popoverEl) return;

    $(popoverEl)
      .off("mouseenter mouseleave")
      .on("mouseenter", () => clearTimer(doiHideTimeout))
      .on("mouseleave", () => !isClickedOpen(el) && scheduleHide(el));
  }

  function buildPopoverContent(data) {
    const { title, "container-title": containerTitle, published, author } = data.message;
    const journal = containerTitle?.[0];
    const year = published?.["date-parts"]?.[0]?.[0];
    const authors = author?.map((a) => a.given).filter(Boolean) || [];
    const doi = $(this).data("doi");

    const parts = [];
    if (journal || year) {
      parts.push(`<p>${journal ? `<em>${journal}</em>` : ""}${year ? ` (${year})` : ""}</p>`);
    }
    if (authors.length) {
      parts.push(`<p class='small'>${authors.join(", ")}</p>`);
    }
    parts.push(`<p><a href="https://doi.org/${doi}" class="btn btn-primary" target="_blank">View paper</a></p>`);

    return { title: title?.[0], content: parts.join("") };
  }

  function createPopover(el, popoverContent, viaClick) {
    const popover = new bootstrap.Popover(el, {
      ...popoverContent,
      html: true,
      trigger: "manual",
      placement: "top",
      // Don't move it around trying to keep it on screen during scroll
      fallbackPlacements: [],
    });
    popover.show();
    if (viaClick) setClickedOpen(el, true);
    setTimeout(() => attachPopoverHoverHandlers(el), POPOVER_MOUNT_DELAY);
    return popover;
  }

  function loadDoiPopover(el, showImmediately = false, viaClick = false) {
    const existing = bootstrap.Popover.getInstance(el);
    if (existing) {
      if (showImmediately) {
        existing.show();
        if (viaClick) setClickedOpen(el, true);
      }
      return;
    }

    if ($(el).hasClass("no-doi-details")) return;

    const doi = $(el).data("doi");
    $.get(`https://api.crossref.org/works/${doi}`)
      .done((data) => {
        const popoverContent = buildPopoverContent.call(el, data);
        if (showImmediately) {
          createPopover(el, popoverContent, viaClick);
        }
      })
      .fail(() => $(el).addClass("no-doi-details"));
  }

  // Event handlers
  $(".module-doi")
    .on("mouseenter", function () {
      clearTimer(doiHideTimeout);
      doiHoverTimeout = setTimeout(() => loadDoiPopover(this, true, false), HOVER_DELAY);
    })
    .on("mouseleave", function () {
      clearTimer(doiHoverTimeout);
      if (isClickedOpen(this)) return;

      setTimeout(() => {
        if (!$(".popover:hover").length) {
          scheduleHide(this);
        }
      }, MOUSELEAVE_CHECK_DELAY);
    })
    .on("click", function (e) {
      e.preventDefault();
      clearTimer(doiHoverTimeout);
      clearTimer(doiHideTimeout);

      if ($(this).hasClass("no-doi-details")) {
        window.open($(this).attr("href"), "_blank");
        return;
      }

      loadDoiPopover(this, true, true);
    });

  // Close on outside click
  $(document).on("click", (e) => {
    if (!$(e.target).closest(".module-doi, .popover").length) {
      $(".module-doi").each(function () {
        const popover = bootstrap.Popover.getInstance(this);
        if (popover) {
          popover.hide();
          setClickedOpen(this, false);
        }
      });
    }
  });
});
