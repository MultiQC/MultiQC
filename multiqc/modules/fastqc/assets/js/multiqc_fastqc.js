// Javascript for the FastQC MultiQC Mod

///////////////
// Per Base Sequence Content
///////////////

// Global vars
s_height = 10;
num_samples = 0;
sample_names = {};
labels = [];
c_width = 0;
c_height = 0;
ypos = 0;
max_bp = 0;
highlight_regex_mode = false;
highlight_f_texts = [];
highlight_f_cols = [];
hidesamples_regex_mode = false;
hidesamples_f_texts = [];
// Function to plot heatmap
function fastqc_seq_content_heatmap() {
    
    // Get sample names, skipping hidden samples
    sample_names = [];
    $.each(Object.keys(fastqc_seq_content_data), function(i, name){
        var hide_sample = false;
        $.each(hidesamples_f_texts, function(idx, f_text){
            if((hidesamples_regex_mode && name.match(f_text))  || (!hidesamples_regex_mode && name.indexOf(f_text) > -1)){
                hide_sample = true;
            }
        });
        if(!hide_sample){ sample_names.push(name); }
    });
    num_samples = sample_names.length;
    if(num_samples == 0){
        $('#fastqc_seq').html('<p class="text-muted">No samples found.</p>');
        return;
    }
    
    // Convert the CSS percentage size into pixels
    c_width = $("#fastqc_seq_heatmap").parent().width();
    c_height = $("#fastqc_seq_heatmap").parent().height();
    s_height = c_height / num_samples;
    // Resize the canvas properties
    $("#fastqc_seq_heatmap").prop("width", c_width);
    $("#fastqc_seq_heatmap").prop("height", c_height+1);
    var canvas = document.getElementById("fastqc_seq_heatmap");
    if (canvas && canvas.getContext) {
        var ctx = canvas.getContext("2d");
        ctx.strokeStyle = "#666666";
        // First, do labels and get max base pairs
        max_bp = 0;
        labels = [];
        $.each(sample_names, function(idx, s_name){
            var s = fastqc_seq_content_data[s_name];
            labels.push(s_name);
            $.each(s, function(bp, v){
                bp = parseInt(bp);
                if(bp > max_bp){
                    max_bp = bp;
                }
            });
        });
        ypos = 0;
        $.each(sample_names, function(idx, s_name){
            var s = fastqc_seq_content_data[s_name]
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
        
        // Draw custom highlights
        $.each(highlight_f_texts, function(idx, f_text){
            var f_col = highlight_f_cols[idx];
            if(f_text == ''){ return true; } // no initial colour, so highlighting all makes no sense
            $.each(labels, function(idx, label){
                if((highlight_regex_mode && label.match(f_text)) || (!highlight_regex_mode && label.indexOf(f_text) > -1)){
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
    
    // Show the status next to the series name for original plots
    $('.mqc-section-fastqc .hc-plot-wrapper').on('mqc_original_chg_source', function(e, name){
        var pid = $(this).find('.hc-plot').attr('id');
        var s_name = get_orig_name(name, pid);
        var status = fastqc_s_statuses[pid][s_name];
        if (status === undefined) { status = '?'; }
        var label = $(this).find('.s_status');
        if(label.length == 0){
            $(this).find('h4').append(' <span class="label label-default s_status">status</span>');
            label = $(this).find('.s_status');
        }
        label.text(status);
        if(status == 'pass'){ label.removeClass().addClass('s_status label label-success'); }
        if(status == 'warn'){ label.removeClass().addClass('s_status label label-warning'); }
        if(status == 'fail'){ label.removeClass().addClass('s_status label label-danger'); }
        if(status == 'status'){ label.removeClass().addClass('s_status label label-default'); }
    });
    
    // Specific behaviour for adapter plot series click
    $('#mqc_fastqc_adapter_plot').on('mqc_original_series_click', function(e, name){
        var snames = name.split(" - ");
        hc_original_chg_source (snames[0], $(this).attr('id'));
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
        if(s_name === undefined){ return false; }
        $('#fastqc_seq .s_name').text(get_new_name(s_name, 'fastqc_seq'));
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
        hc_original_chg_source (name, 'fastqc_seq');
        $("#fastqc_seq .showhide_orig").delay(100).slideDown();
        $("#fastqc_seq_heatmap_div").delay(100).slideUp();
    });
    // Replot heatmap again (clicking the original)
    $('#fastqc_seq img.original-plot').click(function(){
        $("#fastqc_seq_heatmap_div").slideDown();
        fastqc_seq_content_heatmap();
    });
      
    // Highlight the custom heatmap
    $(document).on('mqc_highlights', function(e, f_texts, f_cols, regex_mode){
        highlight_regex_mode = regex_mode;
        highlight_f_texts = f_texts;
        highlight_f_cols = f_cols;
        fastqc_seq_content_heatmap();
    });
    
    // Hide samples the custom heatmap
    $(document).on('mqc_hidesamples', function(e, f_texts, regex_mode){
        hidesamples_regex_mode = regex_mode;
        hidesamples_f_texts = f_texts;
        fastqc_seq_content_heatmap();
    });
    
    // Rename samples on the custom heatmap
    $(document).on('mqc_renamesamples', function(e, f_texts, t_texts){
        $('#fastqc_seq .s_name').text(get_new_name($('#fastqc_seq .s_name').text(), 'fastqc_seq'));
        fastqc_seq_content_heatmap();
    });
    
    // Seq content heatmap drag handle resized
    $('#fastqc_seq').on('mqc_plotresize', function(){
        fastqc_seq_content_heatmap();
    });
    // Seq content - window resized
    $(window).resize(function() {
        fastqc_seq_content_heatmap();
    });

});


// Find the position of the mouse cursor over the canvas
// http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
function findPos(obj) {
  var curleft = 0, curtop = 0;
  if (obj.offsetParent) {
    do {
      curleft += obj.offsetLeft;
      curtop += obj.offsetTop;
    } while (obj = obj.offsetParent);
    return { x: curleft, y: curtop };
  }
  return undefined;
}