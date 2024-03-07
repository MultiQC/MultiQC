""" MultiQC plugin hooks. Enables MultiQC plugins
to run their own custom subroutines at predefined
trigger points during MultiQC execution. """

from importlib_metadata import entry_points

# Load the hooks
hook_functions = {}
for entry_point in entry_points(group="multiqc.hooks.v1"):
    try:
        hook_functions[entry_point.name].append(entry_point.load())
    except KeyError:
        hook_functions[entry_point.name] = [entry_point.load()]


# Function to run the hooks
def mqc_trigger(trigger):
    for hook in hook_functions.get(trigger, []):
        hook()
