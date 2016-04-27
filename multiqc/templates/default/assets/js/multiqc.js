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

$(function () {
  
  // Enable the bootstrap tooltip hovers
  $('[data-toggle="tooltip"]').tooltip();
  
  // Side nav expansion
  $('#side-nav-handle').click(function(e){
    $('.mainpage, .side-nav, .footer').toggleClass('hidden-nav');
    $('#side-nav-handle span').toggleClass('glyphicon-triangle-left glyphicon-triangle-right');
    // send resize trigger for replotting after css animation
    setTimeout(function(){ $(document).resize(); }, 510);
  });
  
  // Hide welcome alert if setting saved
  var hide_welcome = localStorage.getItem("mqc_hide_welcome");
  if(hide_welcome !== 'true'){
    $('#mqc_header_hr').slideUp();
    $('#mqc_welcome').slideDown();
  }
  $('#mqc_hide_welcome_btn').click(function(e){
    localStorage.setItem("mqc_hide_welcome", 'true');
  });
  $('#mqc_hide_welcome_btn, #mqc_welcome .close').click(function(e){
    $('#mqc_header_hr').show();
  });
  
});