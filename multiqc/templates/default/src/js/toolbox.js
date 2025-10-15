////////////////////////////////////////////////
// MultiQC Report Toolbox Code - Main Coordinator
////////////////////////////////////////////////

//////////////////////////////////////////////////////
// TOOLBOX LISTENERS
//////////////////////////////////////////////////////
$(function () {
  const toolboxOffcanvasDiv = document.getElementById("mqc-toolbox");
  const toolboxOffcanvas = new bootstrap.Offcanvas("#mqc-toolbox");

  // Show the toolbox when a button is clicked
  $(".mqc-toolbox-buttons a").click(function (e) {
    if (!toolboxOffcanvasDiv.classList.contains("show")) {
      toolboxOffcanvas.show();
    }
    // Show the tab
    const tabTrigger = new bootstrap.Tab(this);
    tabTrigger.show();
  });

  // Handle mobile navigation toolbox links
  $(".mobile-nav-toolbox-btns a[href^='#mqc_']").click(function (e) {
    e.preventDefault();
    const targetTab = $(this).attr("href");

    // Show the toolbox
    if (!toolboxOffcanvasDiv.classList.contains("show")) {
      toolboxOffcanvas.show();
    }

    // Find the actual tab trigger button in the toolbox that corresponds to this tab
    const actualTabTrigger = $(`.mqc-toolbox-buttons a[href="${targetTab}"]`)[0];
    const tabTrigger = new bootstrap.Tab(actualTabTrigger);
    tabTrigger.show();

    // Close mobile navigation
    const mobileNav = document.getElementById("mqc-nav-collapse");
    if (mobileNav && mobileNav.classList.contains("show")) {
      const bsCollapse = new bootstrap.Collapse(mobileNav);
      bsCollapse.hide();
    }
  });

  // Listener when toolbox is hidden
  toolboxOffcanvasDiv.addEventListener("hidden.bs.offcanvas", (event) => {
    // Show toast if if unsaved changes
    mqc_toolbox_confirmapply();
    // Remove active class from all tabs
    $(".mqc-toolbox-buttons .list-group-item").removeClass("active");
    $("#mqc-toolbox .tab-pane").removeClass("active show");
  });

  // Hide toolbox when a modal is shown
  $(".modal").on("show.bs.modal", function (e) {
    if (toolboxOffcanvasDiv.classList.contains("show")) {
      toolboxOffcanvas.hide();
    }
  });

  // Listener to re-plot graphs if config loaded
  $(document).on("mqc_config_loaded", function (e) {
    $(".hc-plot:not(.not_rendered)").each(function () {
      let target = $(this).attr("id");
      renderPlot(target);
    });
  });

  // Initialize all modules
  initHighlights();
  initRename();
  initHideSamples();
  initSaveLoad();
  initExport();
  initAI();
  initFilters();
  initCitations();
  initHelp();
  initAICookies();
});
