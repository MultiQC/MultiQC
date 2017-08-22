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
from multiqc.utils import config

template_parent = 'default'

base64_plots = False

template_dir = os.path.dirname(__file__)
base_fn = 'base.html'

output_subdir = 'multiqc_report'
copy_files = ['assets']

# This has already been done in the main script, do now if it was false
if not config.export_plots:
    tmp_dir = config.data_tmp_dir.rstrip('multiqc_data')
    config.plots_tmp_dir = os.path.join(tmp_dir, 'multiqc_plots')
    config.plots_dir = config.plots_tmp_dir
    if not os.path.exists(config.plots_dir):
        os.makedirs(config.plots_dir)
config.export_plots = True
