////////////////////////////////////////////////
// MultiQC Report Help Modal Functionality
////////////////////////////////////////////////

// Make function available globally
window.initHelp = function () {
  /////////////////////////
  // REGEX HELP MODAL
  /////////////////////////
  $(".regex_example_buttons button").click(function (e) {
    e.preventDefault();
    $(".regex_example_demo input").val($(this).data("example"));
    regex_example_test();
  });
  $(".regex_example_demo input").keyup(function (e) {
    regex_example_test();
  });

  function regex_example_test() {
    var re = $(".regex_example_demo input").val();
    console.log("Testing " + re);
    $(".regex_example_demo pre span").each(function () {
      // Remove any existing highlighting
      $(this).removeClass();
      $(this).find("hl").contents().unwrap();

      var text = $(this).text();
      var match = text.match(re);
      if (match) {
        console.log("Matches " + text);
        var matchStart = text.indexOf(match[0]);
        var beforeMatch = text.substring(0, matchStart);
        var matchText = match[0];
        var afterMatch = text.substring(matchStart + matchText.length);
        $(this).html(
          $("<span>").text(beforeMatch).addClass("text-muted").prop("outerHTML") +
            $("<hl>").text(matchText).addClass("mark text-success").prop("outerHTML") +
            $("<span>").text(afterMatch).addClass("text-muted").prop("outerHTML"),
        );
      } else {
        console.log("Matches " + text);
        $(this).addClass("text-muted");
      }
    });
  }
};
