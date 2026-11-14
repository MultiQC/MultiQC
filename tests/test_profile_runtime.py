from multiqc import report
from multiqc.core.special_case_modules.profile_runtime import MultiqcModule


def test_profile_runtime_bar_plots_use_single_category():
    """One category per search pattern / module made the plot data grow quadratically."""
    n = 30
    report.runtimes.sp = {f"mod{i}/pattern": 0.1 * i for i in range(n)}
    report.runtimes.mods = {f"mod{i}": 0.1 * i for i in range(n)}
    report.runtimes.total_sp = sum(report.runtimes.sp.values())
    report.runtimes.total_mods = sum(report.runtimes.mods.values())
    report.file_search_stats = {}

    MultiqcModule()
    for plot_id in ("multiqc_runtime_search_patterns_plot", "multiqc_runtime_modules_plot"):
        plot = report.plot_by_id[plot_id]
        assert len(plot.datasets[0].cats) == 1
