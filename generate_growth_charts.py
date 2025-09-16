#!/usr/bin/env python

# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "plotly",
#     "pydriller",
# ]
# ///

import datetime
import fnmatch
import json
import plotly.graph_objects as go
from pydriller import Repository, ModificationType
import re

# Find when each new module was added to MultiQC
modules = {}
mods_plot_x = []
mods_plot_y = []

committers = {}
committers_plot_x = []
committers_plot_y = []
for commit in Repository(".").traverse_commits():
    # Count new modules
    for modification in commit.modified_files:
        if modification.change_type == ModificationType.ADD:
            if fnmatch.fnmatch(modification.new_path, "multiqc/*"):
                mod_match = re.match(
                    r"multiqc/modules/([^/\.]+)", modification.new_path
                )
                if mod_match:
                    mod = mod_match.group(1)
                    if mod not in modules:
                        modules[mod] = str(commit.committer_date)
                        # Plotting data points
                        mods_plot_x.append(commit.committer_date)
                        mods_plot_y.append(len(modules))
    # Count new committers
    if commit.committer.name not in committers:
        committers[commit.committer.name] = str(commit.committer_date)
        # Plotting data points
        committers_plot_x.append(commit.committer_date)
        committers_plot_y.append(len(committers))

print(json.dumps(modules, indent=4))
print(json.dumps(committers, indent=4))

# Make a graph
fig = go.Figure()
fig.add_trace(
    go.Scatter(x=mods_plot_x, y=mods_plot_y, fill="tozeroy")
)  # fill down to xaxis
fig.update_layout(
    title="MultiQC modules over time",
    xaxis_title="Date",
    yaxis_title="Number of modules",
    font=dict(family="Open Sans", size=18, color="#ffffff"),
    plot_bgcolor="rgba(0, 0, 0, 0)",
    paper_bgcolor="rgba(0, 0, 0, 0)",
)
fig.show()
fig.write_image("modules_over_time.svg")

# Add BOSC 2017
bosc_datetime = datetime.datetime.strptime("2017-07-22", "%Y-%m-%d")
fig.add_trace(go.Scatter(x=[bosc_datetime, bosc_datetime], y=[0, 50], mode="lines"))
fig.show()
fig.write_image("modules_over_time_bosc_label.svg")

# Committers over time
fig = go.Figure()
fig.add_trace(
    go.Scatter(x=committers_plot_x, y=committers_plot_y, fill="tozeroy")
)  # fill down to xaxis
fig.update_layout(
    title="MultiQC code contributors over time",
    xaxis_title="Date",
    yaxis_title="Number of contributors",
    font=dict(family="Open Sans", size=18, color="#ffffff"),
    plot_bgcolor="rgba(0, 0, 0, 0)",
    paper_bgcolor="rgba(0, 0, 0, 0)",
)
fig.show()
fig.write_image("code_contributors_over_time.svg")