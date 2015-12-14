////////////////////////////////////////////////
// Base JS for MultiQC Reports
////////////////////////////////////////////////

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
    $('.mainpage, .side-nav, .footer').toggleClass('hidden-nav');
    $('#side-nav-handle span').toggleClass('glyphicon-triangle-left glyphicon-triangle-right');
    // send resize trigger for replotting after css animation
    setTimeout(function(){ $(document).resize(); }, 510);
  });
});