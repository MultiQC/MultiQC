"""
=========
 default
=========

The main MultiQC report template, lovingly known to its admirers
simply as "default"

Note, this is where most of the MultiQC report interactive functionality
is based and will be developed. Unless you want to do some really radical
changes, you probably don't want to replace this theme. Instead, you can
create a child theme that starts with 'default' and then overwrites
certain files.

For more information about creating child themes, see the docs:
docs/templates.md

"""
import os
from .plots import bar, line, scatter, heatmap, violin, table
from ...utils import config

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

bargraph = bar.plot
linegraph = line.plot
scatter = scatter.plot
heatmap = heatmap.plot
violin = violin.plot
table = table.plot


__all__ = ["bargraph", "linegraph", "scatter", "heatmap", "violin", "table"]

if config.development:
    copy_files = ["assets"]
