import logging

from multiqc.plots import table_object
from .violin import plot as plot_violin
from multiqc.utils import config

logger = logging.getLogger(__name__)


# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def plot(dt: table_object.DataTable) -> str:
    s_names = set()
    for d in dt.data:
        for s_name in d.keys():
            s_names.add(s_name)

    if len(s_names) >= config.max_table_rows and dt.pconfig.get("no_beeswarm") is not True:
        logger.debug(f"Plotting violin instead of table, {len(s_names)} samples")
        return plot_violin(dt)
    else:
        # Plot a table with a switch to a violin plot
        return plot_violin(dt, show_table_by_default=True)
