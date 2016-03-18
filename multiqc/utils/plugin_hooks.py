#!/usr/bin/env python

""" MultiQC plugin hooks. Enables MultiQC plugins
to run their own custom subroutines at predefined
trigger points during MultiQC execution. """

import pkg_resources

# Load the hooks
hook_functions = {}
for entry_point in pkg_resources.iter_entry_points('multiqc.hooks.v1'):
  nicename = str(entry_point).split('=')[0].strip()
  try:
    hook_functions[nicename].append(entry_point.load())
  except KeyError:
    hook_functions[nicename] = [entry_point.load()]

# Function to run the hooks
def mqc_trigger (trigger):
  for hook in hook_functions.get(trigger, []):
    hook()
