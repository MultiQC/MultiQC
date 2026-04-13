////////////////////////////////////////////////
// MultiQC Report Citations and DOI Management
////////////////////////////////////////////////

// Make function available globally
window.initCitations = function () {
  // Find DOIs and modules
  var doi_list = { "10.1093/bioinformatics/btw354": "MultiQC" };
  $(".module-doi").each(function () {
    var module_id = $(this).closest(".mqc-module-section-first").find("h2").attr("id");
    doi_list[$(this).data("doi")] = module_id;
  });

  // Build for export toolbox page
  for (var doi in doi_list) {
    $("#mqc_citations_list").append(`
      <tr>
        <th>${doi_list[doi]}<td>
        <a href="https://dx.doi.org/${doi}" target="_blank">${doi}</a></td>
      </tr>
    `);
  }

  // Download DOIs
  $(".download-citations-btn").click(function (e) {
    e.preventDefault();
    var format = $(this).data("format");
    // Get BibTeX
    if (format == "bibtex") {
      var bibtex_string = "";
      // Kick off crossref api calls
      var ajax_promises = [];
      for (var doi in doi_list) {
        ajax_promises.push(
          $.get("https://api.crossref.org/works/" + doi + "/transform/application/x-bibtex", function (data) {
            bibtex_string += data + "\n";
          }),
        );
      }
      // Wait until all API calls are done
      $.when.apply(null, ajax_promises).then(function () {
        var blob = new Blob([bibtex_string], { type: "text/plain;charset=utf-8" });
        saveAs(blob, "multiqc_references.bib");
      });
    }
    // Download list of DOIs
    else {
      var doi_string = "";
      for (var doi in doi_list) {
        doi_string += doi + new Array(50 - doi.length).join(" ") + " # " + doi_list[doi] + "\n";
      }
      var blob = new Blob([doi_string], { type: "text/plain;charset=utf-8" });
      saveAs(blob, "multiqc_dois.txt");
    }
  });
};
