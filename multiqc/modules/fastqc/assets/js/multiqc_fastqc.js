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
current_single_plot = undefined;

// Function to plot heatmap
function fastqc_seq_content_heatmap() {
    
    // Get sample names, rename and skip hidden samples
    sample_names = [];
    sample_statuses = [];
    var p_data = {};
    var hidden_samples = 0;
    $.each(fastqc_seq_content_data, function(s_name, data){
        // rename sample names
        var t_status = fastqc_passfails['per_base_sequence_content'][s_name];
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
        for (i = 0; i < window.mqc_hide_f_texts.length; i++) {
            var f_text = window.mqc_hide_f_texts[i];
            if(window.mqc_hide_regex_mode){
                if(s_name.match(f_text)){ hide_sample = true; }
            } else {
                if(s_name.indexOf(f_text) > -1){ hide_sample = true; }
            }
        }
        if(window.mqc_hide_mode == 'show'){
            hide_sample = !hide_sample;
        }
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
                // Very old versions of FastQC give counts instead of percentages
                if(v['t'] > 100){
                    var t = v['t'] + v['a'] + v['c'] + v['g'];
                    v['t'] = (v['t']/t)*100;
                    v['a'] = (v['a']/t)*100;
                    v['c'] = (v['c']/t)*100;
                    v['g'] = (v['g']/t)*100;
                }
                var r = (v['t'] / 100)*255;
                var g = (v['a'] / 100)*255;
                var b = (v['c'] / 100)*255;
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
            <div class="progress-bar progress-bar-success" style="width: '+(v['pass']/total)*100+'%" title="'+v['pass']+'&nbsp;/&nbsp;'+total+' samples passed">'+v['pass']+'</div> \
            <div class="progress-bar progress-bar-warning" style="width: '+(v['warn']/total)*100+'%" title="'+v['warn']+'&nbsp;/&nbsp;'+total+' samples with warnings">'+v['warn']+'</div> \
            <div class="progress-bar progress-bar-danger" style="width: '+(v['fail']/total)*100+'%" title="'+v['fail']+'&nbsp;/&nbsp;'+total+' samples failed">'+v['fail']+'</div> \
        </div>';
        $(pid).append(p_bar);
    });
    
    // Create popovers on click
    $('.mqc-section-fastqc .fastqc_passfail_progress .progress-bar').mouseover(function(){
        // Does this element already have a popover?
        if ($(this).attr('data-original-title')) { return false; }
        // Create it
        var pid = $(this).closest('h3').attr('id');
        var k = pid.substr(7);
        var vals = fastqc_passfails[k];
        var passes = $(this).hasClass('progress-bar-success') ? true : false;
        var warns = $(this).hasClass('progress-bar-warning') ? true : false;
        var fails = $(this).hasClass('progress-bar-danger') ? true : false;
        var pclass = '';
        if(passes) pclass = 'success';
        if(warns) pclass = 'warning';
        if(fails) pclass = 'danger';
        var samples = Array();
        $.each(vals, function(s_name, status){
            if(status == 'pass' && passes) samples.push(s_name);
            else if(status == 'warn' && warns) samples.push(s_name);
            else if(status == 'fail' && fails) samples.push(s_name);
        });
        $($(this)).popover({
            title: $(this).attr('title'),
            content: samples.sort().join('<br>'),
            html: true,
            trigger: 'hover click focus',
            placement: 'bottom auto',
            template: '<div class="popover popover-'+pclass+'" role="tooltip"> \
                <div class="arrow"></div>\
                <h3 class="popover-title"></h3>\
                <div class="fastqc-popover-intro">\
                    Click bar to fix in place <br>\
                    <a href="#" class="fastqc-status-highlight"><span class="glyphicon glyphicon-pushpin"></span> Highlight these samples</a><br>\
                    <a href="#" class="fastqc-status-hideothers"><span class="glyphicon glyphicon-eye-close"></span> Show only these samples</a>\
                </div>\
                <div class="popover-content"></div>\
            </div>'
        }).popover('show');
    });
    
    // Listener for Status higlight click
    $('.mqc-section-fastqc .fastqc_passfail_progress').on('click', '.fastqc-status-highlight', function(e){
        e.preventDefault();
        // Get sample names and highlight colour
        var samples = $(this).parent().parent().find('.popover-content').html().split('<br>');
        var f_col = mqc_colours[mqc_colours_idx];
        // Add sample names to the toolbox
        for (i = 0; i < samples.length; i++) {
            var f_text = samples[i];
            $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'"/><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
        }
        // Apply highlights and open toolbox
        apply_mqc_highlights();
        mqc_toolbox_openclose('#mqc_cols', true);
        // Update next highlight colour
        mqc_colours_idx += 1;
        if(mqc_colours_idx >= mqc_colours.length){ mqc_colours_idx = 0; }
        $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
        // Hide the popover
        $(this).closest('.popover').popover('hide');
    });
    
    // Listener for Status hide others click
    $('.mqc-section-fastqc .fastqc_passfail_progress').on('click', '.fastqc-status-hideothers', function(e){
        e.preventDefault();
        // Get sample names
        var samples = $(this).parent().parent().find('.popover-content').html().split('<br>');
        // Check if we're already hiding anything, remove after confirm if so
        if($('#mqc_hidesamples_filters li').length > 0){
            if(!confirm($('#mqc_hidesamples_filters li').length+' Hide filters already exist - discard?')){
                return false;
            } else {
                $('#mqc_hidesamples_filters').empty();
            }
        }
        // Set to "show only" and disable regex
        $('.mqc_hidesamples_showhide[value="show"]').prop("checked", true);
        $('#mqc_hidesamples .mqc_regex_mode .re_mode').removeClass('on').addClass('off').text('off');
        // Add sample names to the toolbox
        for (i = 0; i < samples.length; i++) {
            var f_text = samples[i];
            $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
        }
        // Apply highlights and open toolbox
        apply_mqc_hidesamples();
        mqc_toolbox_openclose('#mqc_hidesamples', true);
        // Hide the popover
        $(this).closest('.popover').popover('hide');
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
        $('#fastqc_per_base_sequence_content_plot .s_name').html(s_name + ' <span class="label s_status '+s_status_class+'">'+s_status+'</span>');
        
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
      $('#fastqc_per_base_sequence_content_plot .s_name').html('<em class="text-muted">rollover for sample name</em>');
      $('#fastqc_seq_heatmap_key_pos').text('-');
      $('#fastqc_seq_heatmap_key_t span').text('-');
      $('#fastqc_seq_heatmap_key_c span').text('-');
      $('#fastqc_seq_heatmap_key_a span').text('-');
      $('#fastqc_seq_heatmap_key_g span').text('-');
    });
    
    // Click sample
    $('#fastqc_seq_heatmap').click(function(e) {
      e.preventDefault();
      // Get label from y position
      var pos = findPos(this);
      var x = e.pageX - pos.x;
      var y = e.pageY - pos.y;
      var idx = Math.floor(y/s_height);
      var s_name = sample_names[idx];
      if(s_name !== undefined){
        plot_single_seqcontent(s_name);
      }
    });
    $('#mqc-module-section-fastqc').on('click', '#fastqc_sequence_content_single_prev', function(e){
        e.preventDefault();
        var idx = sample_names.indexOf(current_single_plot) - 1;
        if(idx < 0){ idx = sample_names.length - 1; }
        plot_single_seqcontent(sample_names[idx]);
    });
    $('#mqc-module-section-fastqc').on('click', '#fastqc_sequence_content_single_next', function(e){
        e.preventDefault();
        var idx = sample_names.indexOf(current_single_plot) + 1;
        if(idx == sample_names.length){ idx = 0; }
        plot_single_seqcontent(sample_names[idx]);
    });
    $('#mqc-module-section-fastqc').on('click', '#fastqc_sequence_content_single_back', function(e){
        e.preventDefault();
        $('#fastqc_per_base_sequence_content_plot').slideDown();
        $('#fastqc_sequence_content_single_wrapper').slideUp(function(){
          $(this).remove();
        });
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

function plot_single_seqcontent(s_name){
  current_single_plot = s_name;
  var data = fastqc_seq_content_data[s_name];
  var plot_data = [
    {'name': '% T', 'data':[]},
    {'name': '% C', 'data':[]},
    {'name': '% A', 'data':[]},
    {'name': '% G', 'data':[]}
  ];
  for (var d in data){
    var base = data[d]['base'].toString().split('-');
    base = parseFloat(base[0]);
    plot_data[0]['data'].push([base, data[d]['t']]);
    plot_data[1]['data'].push([base, data[d]['c']]);
    plot_data[2]['data'].push([base, data[d]['a']]);
    plot_data[3]['data'].push([base, data[d]['g']]);
  }
  
  // Create plot div if it doesn't exist, and hide overview
  if($('#fastqc_sequence_content_single_wrapper').length == 0) {
    $('#fastqc_per_base_sequence_content_plot').slideUp();
    var newplot = '<div id="fastqc_sequence_content_single_wrapper"> \
    <div id="fastqc_sequence_content_single_controls"><div class="btn-group"> \
      <button class="btn btn-default btn-sm" id="fastqc_sequence_content_single_prev">&laquo; Prev</button> \
      <button class="btn btn-default btn-sm" id="fastqc_sequence_content_single_next">Next &raquo;</button> \
    </div> <button class="btn btn-default btn-sm" id="fastqc_sequence_content_single_back">Back to overview heatmap</button></div>\
    <div class="hc-plot-wrapper"><div id="fastqc_sequence_content_single" class="hc-plot hc-line-plot"><small>loading..</small></div></div></div>';
    $(newplot).insertAfter('#fastqc_per_base_sequence_content_plot').hide().slideDown();
  }

  $('#fastqc_sequence_content_single').highcharts({
    chart: {
      type: 'line',
      zoomType: 'x'
    },
    colors: ['#dc0000', '#0000dc', '#00dc00', '#404040'],
    title: {
      text: 'Per Base Sequence Content: '+s_name,
      x: 30 // fudge to center over plot area rather than whole plot
    },
    xAxis: {
      title: { text: 'Position (bp)' },
      allowDecimals: false,
    },
    yAxis: {
      title: { text: '% Reads' },
      max: 100,
      min: 0,
    },
    legend: {
      floating: true,
      layout: 'vertical',
      align: 'right',
      verticalAlign: 'top',
      y: 40
    },
    tooltip: {
      backgroundColor: '#FFFFFF',
      borderColor: '#CCCCCC',
      formatter: function () {
        var texts = [];
        var bars = [];
        $.each(this.points, function () {
          texts.push('<span style="display: inline-block; border-left: 3px solid '+this.color+'; padding-left:5px; margin-bottom: 2px;"></div>' + this.y.toFixed(1) + this.series.name + '</span>');
          bars.push('<div class="progress-bar" style="width:'+this.y+'%; float:left; font-size:8px; line-height:12px; padding:0; background-color:'+this.color+';">'+this.series.name.replace('%','').trim()+'</div>');
        });
        return'<p style="font-weight:bold; text-decoration: underline;">Position: ' + this.x + ' bp</p>\
            <p>'+texts.join('<br>')+'</p><div class="progress" style="height: 12px; width: 150px; margin:0;">'+bars.join('')+'</div>';
      },
			useHTML: true,
      crosshairs: true,
      shared: true,
    },
    plotOptions: {
      series: {
        animation: false,
        marker: { enabled: false },
      }
    },
    series: plot_data
  });
}


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