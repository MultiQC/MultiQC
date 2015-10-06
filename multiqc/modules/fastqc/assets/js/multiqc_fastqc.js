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
    var hidden_samples = 0;
    $.each(Object.keys(fastqc_seq_content_data), function(i, name){
        var hide_sample = false;
        $.each(hidesamples_f_texts, function(idx, f_text){
            if((hidesamples_regex_mode && name.match(f_text))  || (!hidesamples_regex_mode && name.indexOf(f_text) > -1)){
                hide_sample = true;
            }
        });
        if(!hide_sample){ sample_names.push(name); }
        else { hidden_samples += 1; }
    });
    num_samples = sample_names.length;
    $('#fastqc_seq_heatmap_div .samples-hidden-warning, #fastqc_seq_heatmap_div .fastqc-heatmap-no-samples').remove();
    $('#fastqc_seq_heatmap_div .hc-plot-wrapper').show();
    if(num_samples == 0){
        $('#fastqc_seq_heatmap_div .hc-plot-wrapper').hide();
        $('#fastqc_seq_heatmap_div').prepend('<p class="fastqc-heatmap-no-samples text-muted">No samples found.</p>');
    }
    if(hidden_samples > 0){
        $('#fastqc_seq_heatmap_div').prepend('<div class="samples-hidden-warning alert alert-danger"><span class="glyphicon glyphicon-info-sign"></span> <strong>Warning:</strong> '+hidden_samples+' samples hidden in toolbox. <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a></div>');
    }
    if(num_samples == 0){ return; }
    
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
    $.each(fastqc_passfails, function(k, vals){
        var pid = '#fastqc_'+k;
        var total = 0;
        var v = { 'pass': 0, 'warn': 0, 'fail': 0 }
        $.each(vals, function(s_name, status){
          total += 1;
          v[status] += 1;
        });
        var p_bar = '<div class="progress fastqc_passfail_progress"> \
            <div class="progress-bar progress-bar-success" style="width: '+(v['pass']/total)*100+'%" title="'+v['pass']+'&nbsp;/&nbsp;'+total+' samples passed" data-toggle="tooltip">'+v['pass']+'</div> \
            <div class="progress-bar progress-bar-warning" style="width: '+(v['warn']/total)*100+'%" title="'+v['warn']+'&nbsp;/&nbsp;'+total+' samples with warnings" data-toggle="tooltip">'+v['warn']+'</div> \
            <div class="progress-bar progress-bar-danger" style="width: '+(v['fail']/total)*100+'%" title="'+v['fail']+'&nbsp;/&nbsp;'+total+' samples failed" data-toggle="tooltip">'+v['fail']+'</div> \
        </div>';
        $(pid).append(p_bar);
        $(pid).tooltip({selector: '[data-toggle="tooltip"]' });
    });
    
    /////////
    /// SEQ CONTENT HEATMAP LISTENERS
    /////////
    
    // Seq Content heatmap mouse rollover
    $("#fastqc_seq_heatmap").mousemove(function(e) {
        
        // Replace the heading above the heatmap
        var pos = findPos(this);
        var x = e.pageX - pos.x;
        var y = e.pageY - pos.y;
        // Get label from y position
        var idx = Math.floor(y/s_height);
        var s_name = labels[idx];
        if(s_name === undefined){ return false; }
        var s_status = fastqc_passfails['sequence_content'][s_name];
        s_status_class = 'label-default';
        if(s_status == 'pass'){ s_status_class = 'label-success'; }
        if(s_status == 'warn'){ s_status_class = 'label-warning'; }
        if(s_status == 'fail'){ s_status_class = 'label-danger'; }
        $('#fastqc_sequence_content_plot .s_name').html(get_new_name(s_name, 'fastqc_seq') + ' <span class="label s_status '+s_status_class+'">'+s_status+'</span>');
        
        // Show the sequence base percentages on the bar plots below
        // http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
        var ctx = this.getContext("2d");
        var p = ctx.getImageData(x, y, 1, 1).data;
        var seq_t = (p[0]/255)*100;
        var seq_a = (p[1]/255)*100;
        var seq_c = (p[2]/255)*100;
        var seq_g = 100 - (seq_t + seq_a + seq_c);
        $("#fastqc_seq_heatmap_key_t span").text(seq_t.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_c span").text(seq_c.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_a span").text(seq_a.toFixed(0)+"%");
        $("#fastqc_seq_heatmap_key_g span").text(seq_g.toFixed(0)+"%");
        
        // Get base pair position from x pos
        var this_bp = Math.floor((x/c_width)*max_bp);
        $('#fastqc_seq_heatmap_key_pos').text(this_bp+' bp');
    });
    
    // Remove sample name again when mouse leaves
    $("#fastqc_seq_heatmap").mouseout(function(e) {
      $('#fastqc_sequence_content_plot .s_name').html('<em class="text-muted">rollover for sample name</em>');
      $('#fastqc_seq_heatmap_key_pos').text('-');
      $("#fastqc_seq_heatmap_key_t span").text('-');
      $("#fastqc_seq_heatmap_key_c span").text('-');
      $("#fastqc_seq_heatmap_key_a span").text('-');
      $("#fastqc_seq_heatmap_key_g span").text('-');
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
        $('#fastqc_sequence_content_plot .s_name').text(get_new_name($('#fastqc_seq .s_name').text(), 'fastqc_seq'));
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