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
      $(this).removeClass();
      if ($(this).text().match(re)) {
        console.log("Matches " + $(this).text());
        $(this).addClass("mark text-success");
      } else {
        console.log("Matches " + $(this).text());
        $(this).addClass("text-muted");
      }
    });
  }
};
