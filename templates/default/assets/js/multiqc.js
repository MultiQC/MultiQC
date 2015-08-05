/* Javascript for MultiQC Default Template */

$('.sample_link').on('click', function(e){
	e.preventDefault();
	var thisid = $(this).text();
	updateId(thisid);
});

$('#next_button').on('click', function(e){
	e.preventDefault();
	var thisid = $('#sample_id').text();
	var regExp = /(\d+)$/;
	var numid = regExp.exec(thisid);
	var newId = parseInt(numid[0]) + 1;
	var name_trunk = thisid.substr(0, thisid.length - numid.length - 1);
	newId = name_trunk + newId;
	updateId(newId);
});

$('#prev_button').on('click', function(e){
	e.preventDefault();
	var thisid = $('#sample_id').text();
	var regExp = /(\d+)$/;
	var numid = regExp.exec(thisid);
	var newId = parseInt(numid[0]) - 1;
	var name_trunk = thisid.substr(0, thisid.length - numid.length - 1);
	newId = name_trunk + newId;
	updateId(newId);
});

function updateId(thisid){
	$('#sample_id').text(thisid);
	var imgbase = 'J.Nielson_14_02_fastq/'+thisid+'/fastqc/Images/';
	$('#quality_histogram .sample_img').attr('src', imgbase+'per_base_quality.png');
	$('#gc_plot .sample_img').attr('src', imgbase+'per_sequence_gc_content.png');
	$('#seq_content_plot .sample_img').attr('src', imgbase+'per_base_sequence_content.png');
}


// Make the project name float when we scroll down
$(window).bind('scroll', function() {
    var scrollY = window.pageYOffset;
	var ylimit = $('.sample_level').offset().top - 50;
    // They scrolled more than the original offset
	if (scrollY > ylimit) {
		$('.sample_header').addClass('fixed');
    } else {
		$('.sample_header').removeClass('fixed');
    }

});
