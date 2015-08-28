// Javascript for the FastQC MultiQC Mod

///////////////
// Per Base Sequence Content
///////////////

// Global vars
s_height = 10;
num_samples = 0;
labels = [];
c_width = 0;
c_height = 0;
ypos = 0;
max_bp = 0;
// Function to plot heatmap
function fastqc_seq_content_heatmap() {
    num_samples = Object.keys(fastqc_seq_content_data).length;
    // Convert the CSS percentage size into pixels
    c_width = $("#fastqc_seq_heatmap").parent().width();
    c_height = $("#fastqc_seq_heatmap").parent().height();
    s_height = c_height / num_samples;
    // Resize the canvas properties
    $("#fastqc_seq_heatmap").prop("width", c_width);
    $("#fastqc_seq_heatmap").prop("height", c_height+1);
    var canvas = document.getElementById("fastqc_seq_heatmap");
    if (canvas.getContext) {
        var ctx = canvas.getContext("2d");
        ctx.strokeStyle = "#666666";
        // First, do labels and get max base pairs
        if(max_bp == 0){
            $.each(fastqc_seq_content_data, function(name, s){
                labels.push(name);
                $.each(s, function(bp, v){
                    bp = parseInt(bp);
                    if(bp > max_bp){
                        max_bp = bp;
                    }
                });
            });
        }
        ypos = 0;
        $.each(fastqc_seq_content_data, function(name, s){
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
        
        // Custom highlights
        $('#hc_col_filters li .f_text').each(function(){
            var f_text = $(this).text();
            var f_col = $(this).css('color');
            if(f_text == '[ all ]'){ return true; } // no initial colour, so highlighting all makes no sense
            $.each(labels, function(idx, label){
                if(label.indexOf(f_text) > -1){
                    var c_width = $("#fastqc_seq_heatmap").width();
                    var ypos = s_height * idx;
                    var canvas = document.getElementById("fastqc_seq_heatmap");
                    var ctx = canvas.getContext("2d");
                    var thisy = ypos + 2;
                    // Adjust height and position if we're next to another square of the same colour
                    // this avoids having a double thickness line
                    if(labels[idx-1] && labels[idx-1].indexOf(f_text) > -1){ thisy = ypos; }
                    ctx.beginPath();
                    ctx.lineWidth="2";
                    ctx.strokeStyle = f_col;
                    ctx.rect(1, thisy, c_width-2, s_height-2);
                    ctx.stroke();
                    ctx.closePath();
                }
            });
        });
    }
}

// Set up listeners etc on page load
$(function () {

    // Add the pass / warning / fails counts to the headings
    $.each(fastqc_passfails, function(k, v){
        var total = v['pass'] + v['warn'] + v['fail'];
        var p_bar = '<div class="progress fastqc_passfail_progress"> \
            <div class="progress-bar progress-bar-success" style="width: '+(v['pass']/total)*100+'%" title="'+v['pass']+'&nbsp;/&nbsp;'+total+' samples passed" data-toggle="tooltip">'+v['pass']+'</div> \
            <div class="progress-bar progress-bar-warning" style="width: '+(v['warn']/total)*100+'%" title="'+v['warn']+'&nbsp;/&nbsp;'+total+' samples with warnings" data-toggle="tooltip">'+v['warn']+'</div> \
            <div class="progress-bar progress-bar-danger" style="width: '+(v['fail']/total)*100+'%" title="'+v['fail']+'&nbsp;/&nbsp;'+total+' samples failed" data-toggle="tooltip">'+v['fail']+'</div> \
        </div>';
        $('#'+k).append(p_bar);
        $('#'+k).tooltip({selector: '[data-toggle="tooltip"]' });
    });

    // Show the overlay plots again (clicking the original)
    $('.original-plot').click(function(){
        $(this).closest('.showhide_orig').next('.fastqc-overlay-plot').slideDown();
        $(this).closest('.showhide_orig').slideUp();
        $(this).closest('.mqc-section').find('.instr').text('Click to show original FastQC plot.');
        highlight_fade_text($(this).closest('.mqc-section').find('.instr'));
        // Replot heatmap
        if ($(this).parents('#fastqc_seq').length) {
            fastqc_seq_content_heatmap();
        }
    });

    // prev / next buttons for original images
    $('.fastqc_prev_btn, .fastqc_nxt_btn').click(function(e){
        e.preventDefault();
        var name = $(this).attr('href').substr(1);
        var target = $(this).data('target');
        fastqc_chg_original (name, target);
    });
    
    /////////
    /// SEQ CONTENT HEATMAP LISTENERS
    /////////
    
    // Show the sequence base percentages on heatmap rollover
    // http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
    $("#fastqc_seq_heatmap").mousemove(function(e) {
        var pos = findPos(this);
        var x = e.pageX - pos.x;
        var y = e.pageY - pos.y;
        // Get label from y position
        var idx = Math.floor(y/s_height);
        var s_name = labels[idx];
        $('#fastqc_seq .s_name').text(s_name);
        var s_status = fastqc_s_statuses["fastqc_seq"][s_name];
        $("#fastqc_seq .s_status").text(s_status);
        if(s_status == 'pass'){ $("#fastqc_seq .s_status").removeClass().addClass('s_status label label-success'); }
        if(s_status == 'warn'){ $("#fastqc_seq .s_status").removeClass().addClass('s_status label label-warning'); }
        if(s_status == 'fail'){ $("#fastqc_seq .s_status").removeClass().addClass('s_status label label-danger'); }
        // Get position from x pos
        var this_bp = Math.floor((x/c_width)*max_bp);
        $('#fastqc_seq_heatmap_key_pos').text(this_bp+' bp');
        // Get colour information
        var ctx = this.getContext("2d");
        var p = ctx.getImageData(x, y, 1, 1).data;
        var seq_t = (p[0]/255)*100;
        var seq_a = (p[1]/255)*100;
        var seq_c = (p[2]/255)*100;
        var seq_g = 100 - (seq_t + seq_a + seq_c);
        if (seq_g < 0){ seq_g = 0; }
        $("#fastqc_seq_heatmap_key_t").text(seq_t.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_c").text(seq_c.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_a").text(seq_a.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_g").text(seq_g.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_colourbar_t span").css("margin-left", seq_t);
        $("#fastqc_seq_heatmap_key_colourbar_c span").css("margin-left", seq_c);
        $("#fastqc_seq_heatmap_key_colourbar_a span").css("margin-left", seq_a);
        $("#fastqc_seq_heatmap_key_colourbar_g span").css("margin-left", seq_g);
    });

    // Show the original plot on click (Sequence Content)
    $("#fastqc_seq_heatmap").click(function(){
        var name = $('#fastqc_seq .s_name').text();
        fastqc_chg_original (name, '#fastqc_seq');
        $("#fastqc_seq .showhide_orig").delay(100).slideDown();
        $("#fastqc_seq_heatmap_div").delay(100).slideUp();
    });
    
      
    // Highlight the custom heatmap
    $(document).on('mqc_highlights:reset', function(){
        fastqc_seq_content_heatmap();
    });
    
    // Seq content heatmap drag handle resized
    $('#fastqc_seq_plot').on('mqc_plotresize', function(){
        fastqc_seq_content_heatmap();
    });
    // Seq content - window resized
    $(window).resize(function() {
        fastqc_seq_content_heatmap();
    });

});

// Update a FastQC original plot
function fastqc_chg_original (name, target) {
    var suffix = $(target+" img.original-plot").data('fnsuffix');
    var names = fastqc_s_names[target.substr(1)];
    var statuses = fastqc_s_statuses[target.substr(1)]
    $(target+" img.original-plot").attr('src', 'report_data/fastqc/'+name+suffix);
    $(target+" .s_name").text(name);
    $(target+" .s_status").text(statuses[name]);
    if(statuses[name] == 'pass'){ $(target+" .s_status").removeClass().addClass('s_status label label-success'); }
    if(statuses[name] == 'warn'){ $(target+" .s_status").removeClass().addClass('s_status label label-warning'); }
    if(statuses[name] == 'fail'){ $(target+" .s_status").removeClass().addClass('s_status label label-danger'); }
    var i = names.indexOf(name);
    var l = names.length;
    var n_i = i+1 < l ? i+1 : 0;
    var p_i = i-1 >= 0 ? i-1 : l - 1;
    var n = names[n_i];
    var p = names[p_i];
    $(target+" .fastqc_prev_btn").attr('href', '#'+p);
    $(target+" .fastqc_nxt_btn").attr('href', '#'+n);
    if($(target).closest('.mqc-section').find('.instr').text() != "Click plot to return to overview plot."){
        $(target).closest('.mqc-section').find('.instr').text("Click plot to return to overview plot.");
        highlight_fade_text($(target).closest('.mqc-section').find('.instr'));
    }
}
