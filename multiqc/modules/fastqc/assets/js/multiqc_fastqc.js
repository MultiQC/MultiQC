// Javascript for the FastQC MultiQC Mod

///////////////
// Per Base Sequence Content
///////////////

// Global vars
s_height = 10;
num_samples = 0;
sample_names = [];
sample_statuses = [];
labels = [];
c_width = 0;
c_height = 0;
ypos = 0;
max_bp = 0;

// Function to plot heatmap
function fastqc_seq_content_heatmap() {
    
    // Get sample names, rename and skip hidden samples
    sample_names = [];
    sample_statuses = [];
    var p_data = {};
    var hidden_samples = 0;
    $.each(fastqc_seq_content_data, function(s_name, data){
        // rename sample names
        var t_status = fastqc_passfails['sequence_content'][s_name];
        $.each(window.mqc_rename_f_texts, function(idx, f_text){
            if(window.mqc_rename_regex_mode){
                var re = new RegExp(f_text,'g');
                s_name = s_name.replace(re, window.mqc_rename_t_texts[idx]);
            } else {
                s_name = s_name.replace(f_text, window.mqc_rename_t_texts[idx]);
            }
        });
        sample_statuses[s_name] = t_status;
        p_data[s_name] = JSON.parse(JSON.stringify(data)); // clone data
        
        var hide_sample = false;
        $.each(window.mqc_hide_f_texts, function(idx, f_text){
            if((window.mqc_hide_regex_mode && s_name.match(f_text))  || (!window.mqc_hide_regex_mode && s_name.indexOf(f_text) > -1)){
                hide_sample = true;
            }
        });
        if(!hide_sample){ sample_names.push(s_name); }
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
        $('#fastqc_seq_heatmap_div').prepend('<div class="samples-hidden-warning alert alert-warning"> \
            <span class="glyphicon glyphicon-info-sign"></span> \
            <strong>Warning:</strong> '+hidden_samples+' samples hidden in toolbox. \
            <a href="#mqc_hidesamples" class="alert-link" onclick="mqc_toolbox_openclose(\'#mqc_hidesamples\', true); return false;">See toolbox.</a>\
        </div>');
    }
    if(num_samples == 0){ return; }
    
    // Convert the CSS percentage size into pixels
    c_width = $('#fastqc_seq_heatmap').parent().width() - 5; // -5 for status bar
    c_height = $('#fastqc_seq_heatmap').parent().height() - 2; // -2 for bottom line padding
    s_height = c_height / num_samples;
    // Resize the canvas properties
    $('#fastqc_seq_heatmap').prop({
        'width': c_width,
        'height': c_height+1
    });
    var canvas = document.getElementById('fastqc_seq_heatmap');
    if (canvas && canvas.getContext) {
        var ctx = canvas.getContext('2d');
        ctx.strokeStyle = '#666666';
        // First, do labels and get max base pairs
        max_bp = 0;
        labels = [];
        $.each(sample_names, function(idx, s_name){
            var s = p_data[s_name];
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
          
            // Add a 5px wide bar indicating either status or Highlight
            var status = sample_statuses[s_name];
            var s_col = '#999999';
            if(status == 'pass'){ s_col = '#5cb85c'; }
            if(status == 'warn'){ s_col = '#f0ad4e'; }
            if(status == 'fail'){ s_col = '#d9534f'; }
            // Override status colour with highlights
            $.each(window.mqc_highlight_f_texts, function(idx, f_text){
                if((window.mqc_highlight_regex_mode && s_name.match(f_text)) || (!window.mqc_highlight_regex_mode && s_name.indexOf(f_text) > -1)){
                  s_col = window.mqc_highlight_f_cols[idx];
                }
            });
            ctx.fillStyle = s_col;
            ctx.fillRect (0, ypos+1, 5, s_height-2);
            
            // Plot the squares for the heatmap
            var s = p_data[s_name]
            var xpos = 6;
            var last_bp = 0;
            $.each(s, function(bp, v){
                bp = parseInt(bp);
                var this_width = (bp - last_bp) * (c_width / max_bp);
                last_bp = bp;
                var r = (v['T'] / 100)*255;
                var g = (v['A'] / 100)*255;
                var b = (v['C'] / 100)*255;
                ctx.fillStyle = chroma(r,g,b).css();
                // width+1 to avoid vertical white line gaps.
                ctx.fillRect (xpos, ypos, this_width+1, s_height);
                xpos += this_width;
            });
            // Draw a line under this row
            ctx.beginPath();
            ctx.moveTo(6, ypos);
            ctx.lineTo(c_width, ypos);
            ctx.stroke();
            ypos += s_height;
        });
        // Final line under row
        ctx.beginPath();
        ctx.moveTo(6, ypos);
        ctx.lineTo(c_width, ypos);
        ctx.stroke();
    }
}

// Set up listeners etc on page load
$(function () {

    // Add the pass / warning / fails counts to each of the FastQC submodule headings
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
    $('#fastqc_seq_heatmap').mousemove(function(e) {
        
        // Replace the heading above the heatmap
        var pos = findPos(this);
        var x = e.pageX - pos.x;
        var y = e.pageY - pos.y;
        // Get label from y position
        var idx = Math.floor(y/s_height);
        var s_name = sample_names[idx];
        if(s_name === undefined){ return false; }
        var s_status = sample_statuses[s_name];
        var s_status_class = 'label-default';
        if(s_status == 'pass'){ s_status_class = 'label-success'; }
        if(s_status == 'warn'){ s_status_class = 'label-warning'; }
        if(s_status == 'fail'){ s_status_class = 'label-danger'; }
        $('#fastqc_sequence_content_plot .s_name').html(s_name + ' <span class="label s_status '+s_status_class+'">'+s_status+'</span>');
        
        // Show the sequence base percentages on the bar plots below
        // http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
        var ctx = this.getContext('2d');
        var p = ctx.getImageData(x, y, 1, 1).data;
        var seq_t = (p[0]/255)*100;
        var seq_a = (p[1]/255)*100;
        var seq_c = (p[2]/255)*100;
        var seq_g = 100 - (seq_t + seq_a + seq_c);
        $('#fastqc_seq_heatmap_key_t span').text(seq_t.toFixed(0)+'%');
        $('#fastqc_seq_heatmap_key_c span').text(seq_c.toFixed(0)+'%');
        $('#fastqc_seq_heatmap_key_a span').text(seq_a.toFixed(0)+'%');
        $('#fastqc_seq_heatmap_key_g span').text(seq_g.toFixed(0)+'%');
        
        // Get base pair position from x pos
        var this_bp = Math.floor((x/c_width)*max_bp);
        $('#fastqc_seq_heatmap_key_pos').text(this_bp+' bp');
    });
    
    // Remove sample name again when mouse leaves
    $('#fastqc_seq_heatmap').mouseout(function(e) {
      $('#fastqc_sequence_content_plot .s_name').html('<em class="text-muted">rollover for sample name</em>');
      $('#fastqc_seq_heatmap_key_pos').text('-');
      $('#fastqc_seq_heatmap_key_t span').text('-');
      $('#fastqc_seq_heatmap_key_c span').text('-');
      $('#fastqc_seq_heatmap_key_a span').text('-');
      $('#fastqc_seq_heatmap_key_g span').text('-');
    });
      
    // Highlight the custom heatmap
    $(document).on('mqc_highlights mqc_hidesamples mqc_renamesamples mqc_plotresize', function(e){
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