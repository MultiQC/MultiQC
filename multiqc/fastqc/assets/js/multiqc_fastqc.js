// Javascript for the FastQC MultiQC Mod

// Per Base Sequence Content
function fastqc_seq_content_heatmap(data) {
    var num_samples = Object.keys(data).length;
    var ypos = 0;
    var max_bp = 0;
    var s_height = 15;
    var labels = [];
    // Convert the CSS percentage size into pixels
    var c_width = $("#fastqc_seq_heatmap").width();
    var c_height = num_samples * (s_height);
    $("#fastqc_seq_heatmap").prop("width", c_width);
    $("#fastqc_seq_heatmap").prop("height", c_height+1);
    var canvas = document.getElementById("fastqc_seq_heatmap");
    if (canvas.getContext) {
        var ctx = canvas.getContext("2d");
        ctx.strokeStyle = "#666666";
        // First, do labels and get max base pairs
        $.each(data, function(name, s){
            labels.push(name);
            $.each(s, function(bp, v){
                bp = parseInt(bp);
                if(bp > max_bp){
                    max_bp = bp;
                }
            });
        });
        $.each(data, function(name, s){
            var xpos = 0;
            var last_bp = 0;
            $.each(s, function(bp, v){
                bp = parseInt(bp);
                var this_width = (bp - last_bp) * (c_width / max_bp);
                last_bp = bp;
                var r = (v["T"] / 100)*255;
                var g = (v["A"] / 100)*255;
                var b = (v["C"] / 100)*255;
                ctx.fillStyle = chroma(r,g,b).css();
                // width+1 to avoid vertical white line gaps.
                ctx.fillRect (xpos, ypos, this_width+1, s_height);
                xpos += this_width;
            });
            // Draw a line under this row
            ctx.beginPath();
            ctx.moveTo(0, ypos);
            ctx.lineTo(c_width, ypos);
            ctx.stroke();
            ypos += s_height;
        });
        // Final line under row
        ctx.beginPath();
        ctx.moveTo(0, ypos);
        ctx.lineTo(c_width, ypos);
        ctx.stroke();
    }

    // Show the sequence base percentages on heatmap rollover
    // http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
    $("#fastqc_seq_heatmap").mousemove(function(e) {
        var pos = findPos(this);
        var x = e.pageX - pos.x;
        var y = e.pageY - pos.y;
        // Get label from y position
        var idx = Math.floor(y/s_height);
        var s_name = labels[idx];
        $('#fastqc_seq_heatmap_sname').text(s_name);
        // Get position from x pos
        var this_bp = Math.floor((x/c_width)*max_bp);
        $('#fastqc_seq_heatmap_key_pos').text(this_bp+' bp');
        // Get colour information
        var ctx = this.getContext("2d");
        var p = ctx.getImageData(x, y, 1, 1).data;
        var seq_t = (p[0]/255)*100;
        var seq_a = (p[1]/255)*100;
        var seq_c = (p[2]/255)*100;
        var seq_g = 100 - (seq_g + seq_a + seq_t);
        if (seq_g < 0){ seq_g = 0; }
        $("#fastqc_seq_heatmap_key_g").text(seq_g.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_a").text(seq_a.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_t").text(seq_t.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_c").text(seq_c.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_colourbar_g span").css("margin-left", seq_g);
        $("#fastqc_seq_heatmap_key_colourbar_a span").css("margin-left", seq_a);
        $("#fastqc_seq_heatmap_key_colourbar_t span").css("margin-left", seq_t);
        $("#fastqc_seq_heatmap_key_colourbar_c span").css("margin-left", seq_c);
    });
}

// Set up listeners etc on page load
$(function () {

  // Show the original plots - Sequence Quality
  $('#fastqc_qual_original img.original-plot').click(function(){
    $('#fastqc_quality_overlay').slideDown();
    $('#fastqc_qual_original').slideUp();
    $('#fastqc_quals_click_instr').text('Click to show original FastQC sequence quality plot.');
  });

  $('#fastqc_qual_orig_nextprev a').click(function(e){
    e.preventDefault();
    var name = $(this).attr('href').substr(1);
    fastqc_chg_original (name, fastqc_overlay_hist_data_names, '_per_base_quality.png', '#fastqc_qual_original', '#fastqc_quals_click_instr');
  });

  // Show the original plots - Sequence Content
  $("#fastqc_seq_heatmap").click(function(){
    var name = $('#fastqc_seq_heatmap_sname').text();
    fastqc_chg_original (name, fastqc_overlay_hist_data_names, '_per_base_sequence_content.png', '#fastqc_seq_original', '#fastqc_seq_heatmap_click_instr');
    $("#fastqc_seq_original").delay(100).slideDown();
    $(this).delay(100).slideUp();
    $('#fastqc_seq_heatmap_key').slideUp();
  });
  $("#fastqc_seq_original .original-plot").click(function(){
    $("#fastqc_seq_heatmap").slideDown();
    $('#fastqc_seq_original').slideUp();
    $('#fastqc_seq_heatmap_key').slideDown();
    $('#fastqc_seq_heatmap_click_instr').text('Click to show original FastQC sequence composition plot.');
  });
  $('#fastqc_seq_orig_nextprev a').click(function(e){
    e.preventDefault();
    var name = $(this).attr('href').substr(1);
    fastqc_chg_original (name, fastqc_overlay_hist_data_names, '_per_base_sequence_content.png', '#fastqc_seq_original', '#fastqc_seq_heatmap_click_instr', '#fastqc_seq_heatmap_sname');
  });

});

// Click event for sequence quality plot
function fastqc_chg_original (name, names, suffix, target, instr_target, name_target) {
    if (name_target == undefined){ name_target = target+" code"; }
    $(target+" img").attr('src', 'report_data/fastqc/'+name+suffix);
    $(name_target).text(name);
    var i = fastqc_overlay_hist_data_names.indexOf(name);
    var l = fastqc_overlay_hist_data_names.length;
    var n_i = i+1 < l ? i+1 : 0;
    var p_i = i-1 >= 0 ? i-1 : l - 1;
    var n = fastqc_overlay_hist_data_names[n_i];
    var p = fastqc_overlay_hist_data_names[p_i];
    $(target+" .prev_btn").attr('href', '#'+p);
    $(target+" .nxt_btn").attr('href', '#'+n);
    $(instr_target).text("Click plot to return to overview plot.");
}
