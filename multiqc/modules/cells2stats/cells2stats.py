import json
import logging

from multiqc.base_module import BaseMultiqcModule

from .cells2stats_bar_plots import plot_cell_segmentation, plot_barcoding, plot_cell_assignment, plot_controls
from .cells2stats_tables import tabulate_wells, tabulate_batches

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="cells2stats",
            anchor="cells2stats",
            href="https://docs.elembio.io/docs/cells2stats/introduction/",
            info="Generate output files and statistics from Element Biosciences Teton Assay",
            doi="",
        )

        self.c2s_run_data = dict()

        observed_run_names = set()

        for f in self.find_log_files("cells2stats/run"):
            data = json.loads(f["f"])
            if data["RunName"] in observed_run_names:
                unique_run_name = f"{data['RunName']} {data['AnalysisID'][-5:]}"
            else:
                unique_run_name = data["RunName"]
                observed_run_names.add(data["RunName"])
            self.c2s_run_data[unique_run_name] = data
            analysis_version = data.get("AnalysisVersion", None)
            if analysis_version is not None:
                self.add_software_version(analysis_version, sample=unique_run_name)
            self.add_data_source(f=f, s_name=unique_run_name, module="cells2stats")

        for plotting_function in [
            tabulate_wells,
            tabulate_batches,
            plot_cell_segmentation,
            plot_barcoding,
            plot_cell_assignment,
            plot_controls,
        ]:
            plot_html, plot_name, anchor, description, helptext, plot_content = plotting_function(self.c2s_run_data)
            self.add_section(
                anchor=anchor,
                name=plot_name,
                helptext=helptext,
                plot=plot_html,
                description=description,
            )
            self.write_data_file(plot_content, anchor)
