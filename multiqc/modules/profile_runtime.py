""" Super Special-Case MultiQC module to produce report section on MultiQC run time """

import logging

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph
from multiqc.utils import report

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Run Time",
            anchor="multiqc_runtime",
            info="""
                This analysis is about the MultiQC run itself, profiling the time spent
                on different parts of the MultiQC execution. It is designed to help
                developers optimise how they run MultiQC, to get the most efficient
                and fastest configuration possible. For more information, see the
                <a href="https://multiqc.info/docs/#optimising-run-time" target="_blank">MultiQC documentation</a>.
            """,
        )

        log.info("Running run time profiling module")

        self.file_search_stats_section()

        self.search_pattern_times_section()

    def file_search_stats_section(self):
        """Count of all files iterated through by MultiQC, by category"""

        pdata = dict()
        pcats = dict()
        for key in sorted(report.file_search_stats, key=report.file_search_stats.get, reverse=True):
            if "skipped_" in key:
                s_name = f"Skipped: {key.replace('skipped_', '').replace('_', ' ').capitalize()}"
                pcats[key] = {"name": key, "color": "#999999"}
            else:
                s_name = key
                pcats[key] = {"name": key, "color": "#7cb5ec"}
            pdata[s_name] = {key: report.file_search_stats[key]}

        pconfig = {
            "id": "multiqc_runtime_files_searched_plot",
            "title": "MultiQC: Files searched",
            "ylab": "Number of files",
            "use_legend": False,
            "cpswitch": False,
        }

        self.add_section(
            name="Files searched",
            anchor="multiqc_runtime_files_searched",
            description="""
                Number of files searched by MultiQC, categorised by what happened to them.
                **Total file searches: {}**.
            """.format(sum(report.file_search_stats.values())),
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
            plot=bargraph.plot(pdata, pcats, pconfig),
        )

    def search_pattern_times_section(self):
        """Section with a bar plot showing the time spent on each search pattern"""

        pdata = dict()
        for key in sorted(report.runtimes["sp"], key=report.runtimes["sp"].get, reverse=True):
            pdata[key] = {"time": report.runtimes["sp"][key]}

        pconfig = {
            "id": "multiqc_runtime_search_patterns_plot",
            "title": "MultiQC: Time per search pattern key",
            "ylab": "Time (seconds)",
            "use_legend": False,
            "cpswitch": False,
        }

        self.add_section(
            name="Search patterns",
            anchor="multiqc_runtime_search_patterns",
            description="""
                Time spent running each search pattern to find files for MultiQC modules.
                **Total file search time: {:.2f} seconds**.
            """.format(report.runtimes["total_sp"]),
            helptext="""
                **NOTE: Usually, MultiQC run time is fairly insignificant - in the order of seconds.
                Unless you are running MultiQC on many thousands of analysis files, optimising this process
                will have limited practical benefit.**

                MultiQC works by recursively looking through all files found in the analysis directories.
                After skipping any that are too big / binary file types etc, it uses the search patterns
                defined in `multiqc/utils/search_patterns.yaml`.
                These work by matching either file names or file contents. Generally speaking, matching
                filenames is super fast and matching file contents is slower.

                Please see the [MultiQC Documentation](https://multiqc.info/docs/#optimising-run-time)
                for information on how to optimise MultiQC to speed this process up.
                The plot below shows which search keys are running and how long each has taken to run in
                total. This should help to guide you to where optimisation is most worthwhile.
            """,
            plot=bargraph.plot(pdata, None, pconfig),
        )
