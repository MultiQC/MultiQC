"""
=============
 sections
=============

This theme is very similar to 'default', except that
the output from each module is hidden by default.
Clicking the module names in the side-navigation
(on the left) hide the visible module and show the
module that was clicked on.

This approach is useful when reports are very large
and take a long time for the browser to render. By only
displaying one piece of the report at a time, the report
is easier to use.

"""
import os
from multiqc.utils import config

template_parent = 'default'
template_dir = os.path.dirname(__file__)
base_fn = 'base.html'

config.collapse_tables = False
