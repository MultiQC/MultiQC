"""Inline-fixture unit tests for the dupblaster module.

The fixtures are dupblaster 0.3.0 output for synthetic query-grouped SAM files with Illumina read names:
two libraries on three flowcell lanes, one library with `--sequencing-duplicate-detection off`, and one
single-end library. The 0.2.0 `--stats` fixture is the same first library in that release's column order.
"""

from typing import Dict, List

import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.dupblaster import MultiqcModule
from multiqc.plots.table_object import InputRow
from multiqc.types import ColumnKey, SampleGroup

METRICS_HEADER = (
    "sample\tlibrary\tdupblaster_version\ttotal_templates\tduplicate_templates\tfrac_duplicates\tmapped_pairs"
    "\tunmapped_pairs\tduplicate_pairs\traw_sequencing_duplicate_pairs\tcorrected_sequencing_duplicate_pairs"
    "\tlibrary_duplicate_pairs\tfrac_duplicate_pairs\tfrac_sequencing_duplicate_pairs\testimated_library_size"
    "\tmapped_orphans\tduplicate_orphans\tunmapped_orphans\tunmated_templates\n"
)
TWO_LIBRARY_METRICS = (
    METRICS_HEADER
    + "sample1\tlibA\t0.3.0\t59693\t18893\t0.321098\t57735\t800\t18556\t3701\t3439\t15117\t0.321399\t0.059565"
    "\t78415\t1158\t337\t0\t0\n"
    "sample1\tlibB\t0.3.0\t37640\t12140\t0.327237\t36400\t500\t11925\t2333\t2110\t9815\t0.327610\t0.057967"
    "\t47819\t740\t215\t0\t0\n"
)
DETECTION_OFF_METRICS = (
    METRICS_HEADER
    + "sample2\tlibC\t0.3.0\t37249\t6649\t0.181748\t35924\t600\t6541\t\t\t\t0.182079\t\t86271\t725\t108\t0\t0\n"
)
SINGLE_END_METRICS = (
    METRICS_HEADER + "sample3\tlibD\t0.3.0\t42024\t11430\t0.275927\t0\t0\t0\t\t\t\t0.000000\t\t\t41424\t11430\t600\t0\n"
)
STATS_0_2_0 = (
    "sample\tlibrary\tdupblaster_version\ttotal_templates\tduplicate_templates\tfrac_duplicates\tmapped_pairs"
    "\tduplicate_pairs\tmapped_orphans\tduplicate_orphans\tunmapped_orphans\tunmapped_pairs\tunmated_templates"
    "\testimated_library_size\n"
    "sample1\tlibA\t0.2.0\t59693\t18893\t0.321098\t57735\t18556\t1158\t337\t0\t800\t0\t64780\n"
)

UNITS_HEADER = (
    "sample\tlibrary\tsequencing_unit\ttemplates\ttiles\tsequencing_duplicate_pairs\tfrac_sequencing_duplicate_pairs\n"
)
TWO_LIBRARY_UNITS = (
    UNITS_HEADER + "sample1\tlibA\tHAAAADSXY:1\t18662\t936\t1159\t0.062105\n"
    "sample1\tlibA\tHAAAADSXY:2\t39073\t936\t2542\t0.065058\n"
    "sample1\tlibB\tHAAAADSXY:2\t22491\t936\t1450\t0.064470\n"
)
SINGLE_END_UNITS = UNITS_HEADER + "sample3\tlibD\t\t0\t0\t0\t0.000000\n"

SAMPLED = (
    "sample\tlibrary\tcategory\ttotal\tunique\tduplicates\tfrac_duplicates\twindow_total\twindow_unique"
    "\twindow_duplicates\twindow_frac_duplicates\n"
    "sample2\tlibC\tpairs\t10000\t9104\t896\t0.089600\t10000\t9104\t896\t0.089600\n"
    "sample2\tlibC\tpairs\t20000\t17524\t2476\t0.123800\t10000\t8420\t1580\t0.158000\n"
    "sample2\tlibC\tpairs\t35924\t29383\t6541\t0.182079\t15924\t11859\t4065\t0.255275\n"
)

SPECTRUM = (
    "sample\tlibrary\tcategory\tn_observations\tn_molecules\n"
    "sample3\tlibD\tsingle_end\t1\t23033\n"
    "sample3\tlibD\tsingle_end\t2\t5303\n"
    "sample3\tlibD\tsingle_end\t9\t2\n"
    "sample3\tlibD\tsingle_end\t10\t1\n"
    "sample3\tlibD\tsingle_end\t412\t1\n"
    "sample3\tlibD\tsingle_end\t929\t1\n"
    "sample3\tlibD\tsingle_end\t948\t1\n"
)


@pytest.fixture
def run_dupblaster(tmp_path):
    """Factory: write each `filename: content` pair to a temp dir and run the dupblaster module on it."""
    # Module tests live outside tests/, so tests/conftest.py's autouse reset does not reach them.
    config.reset()

    def _run(files: Dict[str, str]) -> MultiqcModule:
        for filename, content in files.items():
            (tmp_path / filename).write_text(content)
        report.reset()
        config.preserve_module_raw_data = True
        report.analysis_files = [tmp_path]
        report.search_files(["dupblaster"])
        return MultiqcModule()

    yield _run
    config.reset()


def _general_stats(sample: str) -> Dict[ColumnKey, object]:
    values: Dict[ColumnKey, object] = {}
    for rows_by_group in report.general_stats_data.values():
        for row in rows_by_group.get(SampleGroup(sample), []):
            values.update(row.data)
    return values


def _general_stats_header(key: str) -> Dict:
    for headers in report.general_stats_headers.values():
        if ColumnKey(key) in headers:
            return headers[ColumnKey(key)]
    raise AssertionError(f"{key!r} not in general stats headers")


def _categories(sample: str) -> Dict[str, float]:
    dataset = report.plot_by_id["dupblaster_categories"].datasets[0]
    index = dataset.samples.index(sample)
    return {cat.name: cat.data[index] for cat in dataset.cats if cat.data[index] == cat.data[index]}


def _unit_rows(sample: str) -> List[InputRow]:
    for dataset in report.plot_by_id["dupblaster_sequencing_units_table"].datasets:
        for section in dataset.dt.section_by_id.values():
            if SampleGroup(sample) in section.rows_by_sgroup:
                return section.rows_by_sgroup[SampleGroup(sample)]
    raise AssertionError(f"{sample!r} not in the sequencing units table")


def test_each_library_of_a_file_gets_its_own_row(run_dupblaster):
    module = run_dupblaster({"sample1.dupblaster.duplicate-metrics.tsv": TWO_LIBRARY_METRICS})

    assert set(module.saved_raw_data["multiqc_dupblaster"]) == {"sample1 (libA)", "sample1 (libB)"}
    stats = _general_stats("sample1 (libA)")
    assert stats[ColumnKey("frac_duplicates")] == pytest.approx(0.321098)
    assert stats[ColumnKey("frac_sequencing_duplicate_pairs")] == pytest.approx(0.059565)
    assert _general_stats_header("frac_duplicates")["modify"](0.321098) == pytest.approx(32.1098)
    assert stats[ColumnKey("estimated_library_size")] == 78415
    assert not _general_stats_header("frac_duplicates").get("hidden")
    assert _general_stats_header("frac_sequencing_duplicate_pairs")["hidden"] is True
    assert [version for _, version in module.versions["dupblaster"]] == ["0.3.0"]


def test_categories_split_duplicates_and_add_up_to_total_templates(run_dupblaster):
    module = run_dupblaster({"sample1.dupblaster.duplicate-metrics.tsv": TWO_LIBRARY_METRICS})

    assert _categories("sample1 (libA)") == {
        "Unique pairs": 57735 - 18556,
        "Unique unpaired": 1158 - 337,
        "Library duplicate pairs": 15117,
        "Sequencing duplicate pairs": 3439,
        "Duplicate unpaired": 337,
        "Unmapped": 800,
    }
    for s_name, d in module.saved_raw_data["multiqc_dupblaster"].items():
        assert sum(_categories(s_name).values()) == d["total_templates"]
    assert module.saved_raw_data["multiqc_dupblaster"]["sample1 (libA)"][
        "frac_duplicate_pairs_sequencing"
    ] == pytest.approx(3439 / 18556)


def test_library_without_the_sequencing_split_shows_unsplit_duplicates(run_dupblaster):
    run_dupblaster({"sample2.dupblaster.duplicate-metrics.tsv": DETECTION_OFF_METRICS})

    assert _categories("sample2")["Duplicate pairs (not split)"] == 6541
    assert "Sequencing duplicate pairs" not in _categories("sample2")
    stats = _general_stats("sample2")
    assert ColumnKey("frac_sequencing_duplicate_pairs") not in stats
    assert stats[ColumnKey("estimated_library_size")] == 86271


def test_single_end_library_has_no_pair_rate(run_dupblaster):
    module = run_dupblaster({"sample3.dupblaster.duplicate-metrics.tsv": SINGLE_END_METRICS})

    d = module.saved_raw_data["multiqc_dupblaster"]["sample3"]
    assert d["frac_duplicate_pairs"] is None
    assert d["frac_duplicates"] == pytest.approx(0.275927)
    assert _categories("sample3")["Unique unpaired"] == 41424 - 11430
    assert _categories("sample3")["Duplicate unpaired"] == 11430
    assert ColumnKey("estimated_library_size") not in _general_stats("sample3")


def test_reads_the_0_2_0_stats_file(run_dupblaster):
    module = run_dupblaster({"sample1.stats.tsv": STATS_0_2_0})

    d = module.saved_raw_data["multiqc_dupblaster"]["sample1"]
    assert d["estimated_library_size"] == 64780
    assert "frac_duplicate_pairs_sequencing" not in d
    assert _categories("sample1")["Duplicate pairs (not split)"] == 18556


def test_lanes_nest_under_a_library_summary_row(run_dupblaster):
    run_dupblaster({"sample1.dupblaster.sequencing-units.tsv": TWO_LIBRARY_UNITS})

    summary, *lanes = _unit_rows("sample1 (libA)")
    assert summary.sample == "sample1 (libA)"
    assert summary.data[ColumnKey("units")].raw == 2
    assert summary.data[ColumnKey("templates")].raw == 18662 + 39073
    assert summary.data[ColumnKey("frac_sequencing_duplicate_pairs")].raw == pytest.approx(
        (1159 + 2542) / (18662 + 39073)
    )
    assert [row.sample for row in lanes] == ["sample1 (libA) HAAAADSXY:1", "sample1 (libA) HAAAADSXY:2"]
    assert [row.sample for row in _unit_rows("sample1 (libB)")] == ["sample1 (libB) HAAAADSXY:2"]


def test_placeholder_unit_row_of_a_library_without_pairs_is_skipped(run_dupblaster):
    with pytest.raises(ModuleNoSamplesFound):
        run_dupblaster({"sample3.dupblaster.sequencing-units.tsv": SINGLE_END_UNITS})


def test_duplication_sampled_draws_window_and_cumulative_rates(run_dupblaster):
    run_dupblaster({"sample2.dupblaster.duplication-sampled.tsv": SAMPLED})

    window, cumulative = report.plot_by_id["dupblaster_duplication_sampled"].datasets
    assert dict(next(line for line in window.lines if line.name == "sample2").pairs) == pytest.approx(
        {10000: 8.96, 20000: 15.8, 35924: 25.5275}
    )
    assert dict(next(line for line in cumulative.lines if line.name == "sample2").pairs) == pytest.approx(
        {10000: 8.96, 20000: 12.38, 35924: 18.2079}
    )


def test_duplication_spectrum_bins_molecules_and_weights_templates(run_dupblaster):
    run_dupblaster({"sample3.dupblaster.duplication-spectrum.tsv": SPECTRUM})

    molecules, templates = report.plot_by_id["dupblaster_duplication_spectrum"].datasets
    by_bin = {cat.name: cat.data[0] for cat in molecules.cats}
    assert {k: by_bin[k] for k in ("1x", "2x", "9x", "10-49x", "100-499x", "500-999x")} == {
        "1x": 23033,
        "2x": 5303,
        "9x": 2,
        "10-49x": 1,
        "100-499x": 1,
        "500-999x": 2,
    }
    weighted = {cat.name: cat.data[0] for cat in templates.cats}
    assert weighted["500-999x"] == 929 + 948
    assert weighted["2x"] == 2 * 5303


def test_empty_sample_column_falls_back_to_the_file_name(run_dupblaster):
    metrics = DETECTION_OFF_METRICS.replace("\nsample2\t", "\n\t")
    module = run_dupblaster({"run7.dupblaster.duplicate-metrics.tsv": metrics})

    assert set(module.saved_raw_data["multiqc_dupblaster"]) == {"run7"}


def test_use_filename_as_sample_name(run_dupblaster):
    config.use_filename_as_sample_name = ["dupblaster"]
    module = run_dupblaster(
        {
            "lane3.dupblaster.duplicate-metrics.tsv": TWO_LIBRARY_METRICS,
            "lane3.dupblaster.sequencing-units.tsv": TWO_LIBRARY_UNITS,
        }
    )

    assert set(module.saved_raw_data["multiqc_dupblaster"]) == {"lane3 (libA)", "lane3 (libB)"}
    assert _unit_rows("lane3 (libA)")[0].sample == "lane3 (libA)"


def test_row_with_the_wrong_number_of_fields_is_an_error(run_dupblaster):
    with pytest.raises(ValueError, match="line 2 has 18 fields, expected 19"):
        run_dupblaster({"bad.tsv": METRICS_HEADER + "sample1\tlibA\t0.3.0" + "\t1" * 15 + "\n"})


def test_header_only_files_find_no_samples(run_dupblaster):
    with pytest.raises(ModuleNoSamplesFound):
        run_dupblaster({"empty.dupblaster.duplicate-metrics.tsv": METRICS_HEADER})
