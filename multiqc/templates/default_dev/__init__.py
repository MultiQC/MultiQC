"""
=============
 default_dev
=============

This theme is visually identical to 'default', however it keeps
the source CSS and JavaScript files separate to the final report
file and imports them like a regular website.

This is primarily to help with development - it's great to have
everything included in one file when using MultiQC (it makes it
easier to share reports), however it can get messy when trying
to debug the report file output.

Note that this template uses a couple of special variables below
to get MultiQC to create a subdirectory for the output, plus
copy required files.

"""
import os

template_parent = 'default'

base64_plots = False

template_dir = os.path.dirname(__file__)
base_fn = 'base.html'

output_subdir = 'multiqc_report'
copy_files = ['assets']