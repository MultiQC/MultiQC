"""
Super Special-Case MultiQC module to produce report section on MultiQC performance
"""

import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table
from multiqc import report, config
from multiqc.plots.plotly.bar import BarPlotConfig
from multiqc.plots.table_object import TableConfig

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        info = """- this analysis is about the MultiQC run itself, profiling the time spent
            on different parts of the MultiQC execution. It is designed to help
            developers optimise how they run MultiQC, to get the most efficient
            and fastest configuration possible. For more information, see the
            <a href="https://multiqc.info/docs/#optimising-run-time" target="_blank">MultiQC documentation</a>"""
        super(MultiqcModule, self).__init__(
            name="Run time " + ("and memory " if config.profile_memory else "") + "profiling",
            anchor="multiqc_runtime",
            info=info,
        )

        self.alert = ""
        if config.profile_memory:
            self.alert = (
                "<div class='alert alert-info'>Note that memory profiling can slow down run times of each module</div>"
            )
        else:
            self.alert = "Note: to enable memory profiling, run MultiQC with <code>--profile-memory</code>. Note that it can skew the run time of each module."

        log.info("Running profiling module")
        self.module_table()
        if config.profile_memory:
            self.module_memory_section()
        self.module_times_section()
        self.search_pattern_times_section()
        self.file_search_counts_section()

    def module_table(self):
        """
        Table with time and memory usage per module
        """
        headers = {
            "run_time": {
                "title": "Run time",
                "description": "Time spent running the module",
                "suffix": "s",
                "format": "{:.2f}",
                "scale": "Oranges",
            },
            "peak_mem": {
                "title": "Peak memory",
                "description": "Peak memory usage during module execution",
                "suffix": " MB",
                "format": "{:.2f}",
                "scale": "Greys",
            },
            "mem_change": {
                "title": "Memory change",
                "description": "Change in memory usage during module execution",
                "suffix": " MB",
                "format": "{:.2f}",
                "scale": "Blues",
            },
        }

        table_data: Dict[str, Dict[str, float]] = {}
        for key in report.runtimes.mods:
            table_data[key] = {
                "run_time": report.runtimes.mods[key],
            }

        if config.profile_memory:
            for key in report.peak_memory_bytes_per_module:
                table_data[key].update(
                    {
                        "peak_mem": report.peak_memory_bytes_per_module[key] / 1024 / 1024,
                        "mem_change": report.diff_memory_bytes_per_module[key] / 1024 / 1024,
                    }
                )

        self.add_section(
            name="Per module",
            anchor="per_module_benchmark",
            description=self.alert,
            plot=table.plot(
                table_data,
                headers,
                pconfig=TableConfig(
                    id="per_module_benchmark_table",
                    title="Module run times" + " and memory usage" if config.profile_memory else "",
                    col1_header="Module",
                ),
            ),
        )

    def file_search_counts_section(self):
        """Count of all files iterated through by MultiQC, by category"""

        file_search_counts: Dict[str, int] = {k: len(paths) for k, paths in report.file_search_stats.items()}

        pdata: Dict[str, Dict] = dict()
        pcats: Dict[str, Dict] = dict()
        for key in sorted(
            file_search_counts.keys(),
            key=lambda k: file_search_counts[k],
            reverse=True,
        ):
            if "skipped_" in key:
                s_name = f"Skipped: {key.replace('skipped_', '').replace('_', ' ').capitalize()}"
                pcats[key] = {"name": key, "color": "#999999"}
            else:
                s_name = key
                pcats[key] = {"name": key, "color": "#7cb5ec"}
            pdata[s_name] = {key: file_search_counts[key]}

        self.add_section(
            name="Files searched counts",
            anchor="multiqc_runtime_files_searched",
            description="""
                Number of files searched by MultiQC, categorised by what happened to them.
                **Total file searches: {}**.
            """.format(sum(file_search_counts.values())),
            helptext="""
                Note that only files are considered in this plot - skipped directories are not shown.

                Some search patterns do not discard files after they match (`shared: true`), so it is possible
                that some files may be double-counted in this plot.

                * `Skipped: No match` - File was searched, but didn't match any search patterns
                * `Skipped: Ignore pattern` - File matched a MultiQC ignore pattern (see `-x` / `--ignore` / `config.fn_ignore_paths`)
                * `Skipped: Filesize limit` - File was skipped because it was too large (see `config.log_filesize_limit`)
                * `Skipped: Symlinks` - File was a symlink and skipped (see `config.ignore_symlinks`)
                * `Skipped: Not a file` - File could not be read (eg. was a unix pipe or something)
            """,
            plot=bargraph.plot(
                pdata,
                pcats,
                BarPlotConfig(
                    id="multiqc_runtime_files_searched_plot",
                    title="MultiQC: Files searched",
                    ylab="Number of files",
                    use_legend=False,
                    cpswitch=False,
                ),
            ),
        )

    def search_pattern_times_section(self):
        """Section with a bar plot showing the time spent on each search pattern"""

        pdata: Dict[str, Dict] = dict()
        for key in sorted(report.runtimes.sp.keys(), key=lambda k: report.runtimes.sp[k], reverse=True):
            pdata[key] = {"Run time": report.runtimes.sp[key]}

        pconfig = {
            "id": "multiqc_runtime_search_patterns_plot",
            "title": "MultiQC: Time per search pattern key",
            "ylab": "Run time",
            "use_legend": False,
            "cpswitch": False,
            "suffix": "s",
        }

        self.add_section(
            name="Search patterns run times",
            anchor="multiqc_runtime_search_patterns",
            description="""
                Time spent running each search pattern to find files for MultiQC modules.
                **Total file search time: {:.2f} seconds**.
            """.format(report.runtimes.total_sp),
            helptext="""
                **NOTE: Usually, MultiQC run time is fairly insignificant - in the order of seconds.
                Unless you are running MultiQC on many thousands of analysis files, optimising this process
                will have limited practical benefit.**

                MultiQC works by recursively looking through all files found in the analysis directories.
                After skipping any that are too big / binary file types etc., it uses the search patterns
                defined in `multiqc/search_patterns.yaml`.
                These work by matching either file names or file contents. Generally speaking, matching
                filenames is super fast and matching file contents is slower.

                Please see the [MultiQC Documentation](https://multiqc.info/docs/#optimising-run-time)
                for information on how to optimise MultiQC to speed this process up.
                The plot below shows which search keys are running and how long each has taken to run in
                total. This should help to guide you to where optimisation is most worthwhile.
            """,
            plot=bargraph.plot(pdata, None, pconfig),
        )

    def module_times_section(self):
        """Section with a bar plot showing the time spent on each search pattern"""

        pdata = dict()
        for key in report.runtimes.mods:
            pdata[key] = {"Time": report.runtimes.mods[key]}

        pconfig = {
            "id": "multiqc_runtime_modules_plot",
            "title": "MultiQC: Time per module",
            "ylab": "Run time",
            "use_legend": False,
            "cpswitch": False,
            "suffix": "s",
        }

        description = f"""
            Time spent running each module.
            **Total modules run time: {report.runtimes.total_mods:.2f} seconds**.
            <br><br>{self.alert}
        """

        self.add_section(
            name="Per module run times",
            anchor="multiqc_runtime_modules",
            description=description,
            plot=bargraph.plot(pdata, None, pconfig),
        )

    def module_memory_section(self):
        """
        Section with a bar plot showing the memory usage of each module
        """
        pdata: Dict[str, Dict] = {}
        for key in report.peak_memory_bytes_per_module:
            pdata[key] = {"Peak memory": report.peak_memory_bytes_per_module[key] / 1024 / 1024}
        for key in report.diff_memory_bytes_per_module:
            pdata[key]["Memory change"] = report.diff_memory_bytes_per_module[key] / 1024 / 1024

        pconfig = {
            "id": "multiqc_runtime_memory_plot",
            "title": "MultiQC: Memory usage per module",
            "ylab": "Memory",
            "suffix": " MB",
            "stacking": "overlay",
        }
        self.add_section(
            name="Per module memory usage",
            anchor="multiqc_runtime_memory",
            description="Memory usage per each module. The <span style='color: #7cb5ec'>blue</span> "
            "bar indicates how much more memory MultiQC occupies after finishing running the module, which roughly should"
            "correspond to the size of the parsed data, which is loaded into memory. "
            "The <span style='color: #888888'>grey</span> bar shows the peak memory usage during the module "
            "execution - some memory could be cleaned after module is finished."
            f"<br><br>{self.alert}",
            plot=bargraph.plot(
                pdata,
                {
                    "Peak memory": {"name": "Peak memory", "color": "#999999"},
                    "Memory change": {"name": "Memory change", "color": "#7cb5ec"},
                },
                pconfig=pconfig,
            ),
        )
