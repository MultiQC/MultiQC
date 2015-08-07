// Javascript for the FastQC MultiQC Mod

// Per Base Sequence Content
function fastqc_seq_content_heatmap(data) {
    var num_samples = Object.keys(data).length;
    var num_data = undefined;
    var ypos = 0;
    var xpos_start = 0;
    var labels_width = 0;
    var max_bp = 0;
    var s_height = 20;
    // Convert the CSS percentage size into pixels
    var c_width = $("#fastqc_seq_heatmap").width();
    var c_height = $("#fastqc_seq_heatmap").height();
    c_height = num_samples * (s_height);
    $("#fastqc_seq_heatmap").prop("width", c_width);
    $("#fastqc_seq_heatmap").prop("height", c_height+1);
    var canvas = document.getElementById("fastqc_seq_heatmap");
    if (canvas.getContext) {{
        var ctx = canvas.getContext("2d");
        ctx.font = "12px Arial";
        ctx.strokeStyle = "#333333";
        // First, do labels and get max base pairs
        $.each(data, function(name, s){{
            var txt_w = ctx.measureText(name).width;
            if(txt_w > labels_width) {{
                labels_width = txt_w;
            }}
            ctx.fillText(name, 0, ypos+14);
            ypos += s_height;
            $.each(s, function(bp, v){{
                bp = parseInt(bp);
                if(bp > max_bp){{
                    max_bp = bp;
                }}
            }});
        }});
        labels_width += 5;
        ypos = 0;
        xpos_start = labels_width;
        c_width -= labels_width;
        $.each(data, function(name, s){{
            if (num_data == undefined){{
                var s_width = c_width / Object.keys(s).length;
            }}
            var xpos = xpos_start;
            var last_bp = 0;
            $.each(s, function(bp, v){{
                bp = parseInt(bp);
                var this_width = (bp - last_bp) * (c_width / max_bp);
                last_bp = bp;
                var r = (v["G"] / 100)*255;
                var g = (v["A"] / 100)*255;
                var b = (v["T"] / 100)*255;
                ctx.fillStyle = chroma(r,g,b).css();
                ctx.fillRect (xpos, ypos, this_width, s_height-1);
                xpos += this_width;
            }});
            // Draw a line under this row
            ctx.beginPath();
            ctx.moveTo(0, ypos);
            ctx.lineTo(c_width+labels_width, ypos);
            ctx.stroke();
            ypos += s_height;
        }});
        // Final line under row
        ctx.beginPath();
        ctx.moveTo(0, ypos);
        ctx.lineTo(c_width+labels_width, ypos);
        ctx.stroke();
    }}

    // http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
    $("#fastqc_seq_heatmap").mousemove(function(e) {{
        var pos = findPos(this);
        var x = e.pageX - pos.x;
        var y = e.pageY - pos.y;
        var coord = "x=" + x + ", y=" + y;
        var c = this.getContext("2d");
        var p = c.getImageData(x, y, 1, 1).data;
        var seq_g = Math.round((p[0]/255)*100);
        var seq_a = Math.round((p[1]/255)*100);
        var seq_t = Math.round((p[2]/255)*100);
        var seq_c = 100 - (seq_g + seq_a + seq_t);
        $("#fastqc_seq_heatmap_key_g").text(seq_g+"%");
        $("#fastqc_seq_heatmap_key_a").text(seq_a+"%");
        $("#fastqc_seq_heatmap_key_t").text(seq_t+"%");
        $("#fastqc_seq_heatmap_key_c").text(seq_c+"%");
        $("#fastqc_seq_heatmap_key_colourbar_g span").css("margin-left", seq_g);
        $("#fastqc_seq_heatmap_key_colourbar_a span").css("margin-left", seq_a);
        $("#fastqc_seq_heatmap_key_colourbar_t span").css("margin-left", seq_t);
        $("#fastqc_seq_heatmap_key_colourbar_c span").css("margin-left", seq_c);
    }});
}
