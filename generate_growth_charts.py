#!/usr/bin/env python

# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "plotly",
#     "pydriller",
# ]
# ///

import fnmatch
import json
import plotly.graph_objects as go
from pydriller import Repository, ModificationType
import re

def extract_coauthors(commit_msg):
    """Extract co-authors from commit message"""
    coauthors = []
    # Look for Co-authored-by: Name <email>
    coauthor_pattern = r'Co-authored-by:\s*([^<\n\r]+?)(?:\s*<[^>]+>)?\s*$'
    matches = re.findall(coauthor_pattern, commit_msg, re.MULTILINE | re.IGNORECASE)
    
    for match in matches:
        # Clean up the name
        name = match.strip()
        # Skip bot accounts, empty names, and single character names
        if (name and len(name) > 2 and 
            not any(bot in name.lower() for bot in ['bot', 'github-actions', 'multiqc bot']) and
            not name.startswith('Co-authored-by')):
            coauthors.append(name)
    
    return coauthors

# Find when each new module was added to MultiQC
modules = {}
mods_plot_x = []
mods_plot_y = []

# Track both committers and co-authors
contributors = {}
contributors_plot_x = []
contributors_plot_y = []

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
    
    # Count new contributors (main committer)
    if commit.committer.name not in contributors:
        contributors[commit.committer.name] = str(commit.committer_date)
        # Plotting data points
        contributors_plot_x.append(commit.committer_date)
        contributors_plot_y.append(len(contributors))
    
    # Count co-authors from commit message
    coauthors = extract_coauthors(commit.msg)
    for coauthor in coauthors:
        if coauthor not in contributors:
            contributors[coauthor] = str(commit.committer_date)
            # Plotting data points
            contributors_plot_x.append(commit.committer_date)
            contributors_plot_y.append(len(contributors))

print(f"Total modules found: {len(modules)}")
print(f"Total contributors found (including co-authors): {len(contributors)}")

# Modules over time
fig = go.Figure()
fig.add_trace(
    go.Scatter(x=mods_plot_x, y=mods_plot_y, fill="tozeroy", name="Modules")
)
fig.update_layout(
    title="MultiQC modules over time",
    xaxis_title="Date",
    yaxis_title="Number of modules",
    font=dict(family="Open Sans", size=18, color="#ffffff"),
    plot_bgcolor="rgba(0, 0, 0, 0)",
    paper_bgcolor="rgba(0, 0, 0, 0)",
)
fig.write_image("modules_over_time.svg")

# Contributors over time (including co-authors from squash merges)
fig = go.Figure()
fig.add_trace(
    go.Scatter(x=contributors_plot_x, y=contributors_plot_y, fill="tozeroy", name="Contributors")
)
fig.update_layout(
    title="MultiQC code contributors over time",
    xaxis_title="Date",
    yaxis_title="Number of contributors",
    font=dict(family="Open Sans", size=18, color="#ffffff"),
    plot_bgcolor="rgba(0, 0, 0, 0)",
    paper_bgcolor="rgba(0, 0, 0, 0)",
)
fig.write_image("contributors_over_time.svg")

print(f"\nGenerated charts:")
print(f"- modules_over_time.svg ({len(modules)} modules)")
print(f"- contributors_over_time.svg ({len(contributors)} contributors, including co-authors from squash merges)")