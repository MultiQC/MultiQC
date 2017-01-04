////////////////////////////////////////////////
// Static MatPlotLib Plots Javascript Code
////////////////////////////////////////////////

// On page load
$(function () {

    // Switch between counts and percentages in a bar plot
    $('.mqc_mplplot_bargraph_setcountspcnt button').click(function(e){
        e.preventDefault();
        if(!$(this).hasClass('active')){
            $(this).siblings('button.active').removeClass('active');
            $(this).addClass('active');
            var wrapper = $(this).closest('.mqc_mplplot_plotgroup');
            var current = '#'+wrapper.find('.mqc_mplplot:visible').attr('id');
            if (current.substr(current.length - 3) == '_pc') {
                var target = current.substr(0, current.length - 3);
            } else {
                var target = current + '_pc';
            }
            wrapper.find('.mqc_mplplot').hide();
            $(target).show();
        }
    });

    // Switch datasets in a bar plot
    $('.mqc_mplplot_bargraph_switchds button').click(function(e){
        e.preventDefault();
        if(!$(this).hasClass('active')){
            $(this).siblings('button.active').removeClass('active');
            $(this).addClass('active');
            var target = $(this).data('target');
            var wrapper = $(target).closest('.mqc_mplplot_plotgroup');
            if(wrapper.find('.mqc_mplplot_bargraph_setcountspcnt .pcnt').hasClass('active')){
                target += '_pc';
            }
            wrapper.find('.mqc_mplplot').hide();
            $(target).show();
        }
    });

});