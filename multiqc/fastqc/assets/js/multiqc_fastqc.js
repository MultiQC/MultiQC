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
                var r = (v["G"] / 100)*255;
                var g = (v["A"] / 100)*255;
                var b = (v["T"] / 100)*255;
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
        var seq_g = (p[0]/255)*100;
        var seq_a = (p[1]/255)*100;
        var seq_t = (p[2]/255)*100;
        var seq_c = 100 - (seq_g + seq_a + seq_t);
        if (seq_c < 0){ seq_c = 0; }
        $("#fastqc_seq_heatmap_key_g").text(seq_g.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_a").text(seq_a.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_t").text(seq_t.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_c").text(seq_c.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_colourbar_g span").css("margin-left", seq_g);
        $("#fastqc_seq_heatmap_key_colourbar_a span").css("margin-left", seq_a);
        $("#fastqc_seq_heatmap_key_colourbar_t span").css("margin-left", seq_t);
        $("#fastqc_seq_heatmap_key_colourbar_c span").css("margin-left", seq_c);
    });

    // Show the original plots
    $("#fastqc_seq_heatmap").click(function(){
      var img = '<img style="max-width:800px; width:100%;" src="report_data/fastqc/'+$('#fastqc_seq_heatmap_sname').text()+'_per_base_sequence_content.png">';
      $("#fastqc_seq_original").html(img).delay(100).slideDown();
      $(this).delay(100).slideUp();
      $('#fastqc_seq_heatmap_key').slideUp();
      $('#fastqc_seq_heatmap_click_instr').text('Click to show the MultiQC sequence composition heatmap.');
    });
    $("#fastqc_seq_original").click(function(){
      $("#fastqc_seq_heatmap").slideDown();
      $(this).slideUp();
      $('#fastqc_seq_heatmap_key').slideDown();
      $('#fastqc_seq_heatmap_click_instr').text('Click to show original FastQC sequence composition plot.');
    });
}
