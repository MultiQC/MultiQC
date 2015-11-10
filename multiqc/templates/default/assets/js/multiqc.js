////////////////////////////////////////////////
// Base JS for MultiQC Reports
////////////////////////////////////////////////

$(function () {
  // Smooth scroll to top
  $("a[href='#top']").click(function(e) {
    $("html, body").animate({ scrollTop: 0 }, "slow");
  });
});

// From http://stackoverflow.com/questions/2901102/how-to-print-a-number-with-commas-as-thousands-separators-in-javascript
function numberWithCommas(x) {
  return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

// Helper config - is defined and object length > 0?
function notEmptyObj (obj){
  try{
    if(obj === undefined){ return false; }
    if(obj.length == 0){ return false; }
  } catch(e){ return false; }
  return true;
}


// Side nav expansion
$(function () {
  $('#side-nav-handle').click(function(e){
    $('.side-nav-sm, .side-nav-main').toggle();
    $('body').toggleClass('side-nav-mini');
    $(document).resize(); // send resize trigger for replotting
  });
});