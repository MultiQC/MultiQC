"""Tests for the OrthoFinder module: parsers and integration."""

import pytest

from multiqc import report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.orthofinder.orthofinder import (
    MultiqcModule,
    parse_number,
)


def _loaded(path):
    """Mimic the dict find_log_files yields, so parsers can be tested directly."""
    return {"f": path.read_text(), "fn": path.name, "root": str(path.parent), "s_name": path.stem}


class TestParseNumber:
    def test_integral_text_gives_int(self):
        assert parse_number("27628") == 27628
        assert isinstance(parse_number("27628"), int)

    def test_decimal_text_gives_float(self):
        assert parse_number("88.5") == pytest.approx(88.5)

    def test_surrounding_whitespace_is_tolerated(self):
        assert parse_number("  42  ") == 42

    def test_non_numeric_raises(self):
        # A date or filename reaching this function means the format changed, and
        # that should surface rather than become a zero.
        with pytest.raises(ValueError):
            parse_number("2026-01-15")


class TestParsePerSpecies:
    def _parse(self, fixtures_dir):
        report.reset()
        m = MultiqcModule.__new__(MultiqcModule)
        path = fixtures_dir / "run_alpha" / "Comparative_Genomics_Statistics" / "Statistics_PerSpecies.tsv"
        return m._parse_per_species(_loaded(path))

    def test_one_entry_per_species(self, fixtures_dir):
        parsed = self._parse(fixtures_dir)
        assert set(parsed) == {"species_a", "species_b", "species_c"}

    def test_values_are_read_per_column(self, fixtures_dir):
        parsed = self._parse(fixtures_dir)
        assert parsed["species_a"]["genes"] == 20000
        assert parsed["species_b"]["genes"] == 18000
        assert parsed["species_c"]["genes"] == 22000
        assert parsed["species_c"]["percent_genes_in_orthogroups"] == pytest.approx(90.0)

    def test_stops_before_the_distribution_block(self, fixtures_dir):
        """The histogram below the blank line must not be read as more metrics."""
        parsed = self._parse(fixtures_dir)
        for data in parsed.values():
            assert len(data) == 10
            assert not any(k.startswith("'") for k in data)

    def test_a_metric_label_below_the_blank_line_is_not_read(self, tmp_path):
        """Parsing stops at the blank line, not merely on unrecognised labels.

        The real distribution block happens to use labels no metric shares, so a
        parser that only filtered by label would look correct on real files while
        still being wrong. This pins the stop itself.
        """
        m = MultiqcModule.__new__(MultiqcModule)
        path = tmp_path / "Statistics_PerSpecies.tsv"
        path.write_text("\tsp_a\nNumber of genes\t10\n\nNumber of genes\t999\n")
        assert m._parse_per_species(_loaded(path)) == {"sp_a": {"genes": 10}}

    def test_empty_file_yields_nothing(self, tmp_path):
        m = MultiqcModule.__new__(MultiqcModule)
        path = tmp_path / "Statistics_PerSpecies.tsv"
        path.write_text("")
        assert m._parse_per_species(_loaded(path)) == {}

    def test_header_without_species_yields_nothing(self, tmp_path):
        m = MultiqcModule.__new__(MultiqcModule)
        path = tmp_path / "Statistics_PerSpecies.tsv"
        path.write_text("\t\t\nNumber of genes\t1\t2\n")
        assert m._parse_per_species(_loaded(path)) == {}

    def test_crlf_line_endings_are_handled(self, tmp_path):
        """OrthoFinder output is CRLF in the wild, and a stray \\r must not survive."""
        m = MultiqcModule.__new__(MultiqcModule)
        path = tmp_path / "Statistics_PerSpecies.tsv"
        path.write_bytes(
            b"\tsp_a\tsp_b\r\nNumber of genes\t10\t20\r\n\r\nNumber of genes per-species in orthogroup\r\n'0\t1\t2\r\n"
        )
        parsed = m._parse_per_species(_loaded(path))
        assert parsed == {"sp_a": {"genes": 10}, "sp_b": {"genes": 20}}


class TestParseOverall:
    def _parse(self, fixtures_dir):
        m = MultiqcModule.__new__(MultiqcModule)
        path = fixtures_dir / "run_alpha" / "Comparative_Genomics_Statistics" / "Statistics_Overall.tsv"
        return m._parse_overall(_loaded(path))

    def test_reads_the_run_level_numbers(self, fixtures_dir):
        data = self._parse(fixtures_dir)
        assert data["species"] == 3
        assert data["orthogroups"] == 9000
        assert data["single_copy_orthogroups"] == 900
        assert data["mean_orthogroup_size"] == pytest.approx(6.0)

    def test_filename_and_date_rows_are_skipped(self, fixtures_dir):
        """Rows whose value is a date or a filename are not run statistics."""
        data = self._parse(fixtures_dir)
        assert "Date" not in data
        assert all(isinstance(v, (int, float)) for v in data.values())
        assert len(data) == 18

    def test_a_metric_label_below_the_blank_line_is_not_read(self, tmp_path):
        """As above: pin the stop at the blank line, not the label allowlist."""
        m = MultiqcModule.__new__(MultiqcModule)
        path = tmp_path / "Statistics_Overall.tsv"
        path.write_text("Number of species\t3\n\nNumber of species\t999\n")
        assert m._parse_overall(_loaded(path)) == {"species": 3}

    def test_stops_before_the_distribution_table(self, fixtures_dir):
        data = self._parse(fixtures_dir)
        assert "Percentage of genes" not in data
        assert "Average number of genes per-species in orthogroup" not in data


class TestRunName:
    def test_named_after_the_run_dir_not_the_stats_dir(self, tmp_path):
        """Every run has a Comparative_Genomics_Statistics dir, so it names nothing."""
        m = MultiqcModule.__new__(MultiqcModule)
        m.clean_s_name = lambda name, f: name
        stats = tmp_path / "Results_Jan01" / "Comparative_Genomics_Statistics"
        stats.mkdir(parents=True)
        path = stats / "Statistics_Overall.tsv"
        path.write_text("Number of species\t3\n")
        assert m._run_name(_loaded(path)) == "Results_Jan01"

    def test_relative_root_resolves_against_cwd(self, tmp_path, monkeypatch):
        """`multiqc .` from inside the run dir hands the module a relative root.

        Without resolving it, the stats dir has no parent to fall back on and every
        run collapses to the same name.
        """
        m = MultiqcModule.__new__(MultiqcModule)
        m.clean_s_name = lambda name, f: name
        run_dir = tmp_path / "Results_Feb02"
        stats = run_dir / "Comparative_Genomics_Statistics"
        stats.mkdir(parents=True)
        (stats / "Statistics_Overall.tsv").write_text("Number of species\t3\n")
        monkeypatch.chdir(run_dir)
        f = {"f": "", "fn": "Statistics_Overall.tsv", "root": "Comparative_Genomics_Statistics", "s_name": "x"}
        assert m._run_name(f) == "Results_Feb02"

    def test_other_directories_are_used_as_is(self, tmp_path):
        m = MultiqcModule.__new__(MultiqcModule)
        m.clean_s_name = lambda name, f: name
        other = tmp_path / "my_run"
        other.mkdir()
        path = other / "Statistics_Overall.tsv"
        path.write_text("Number of species\t3\n")
        assert m._run_name(_loaded(path)) == "my_run"


class TestModuleIntegration:
    def test_module_runs_on_fixture_data(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir)]
        report.search_files(["orthofinder"])

        m = MultiqcModule()

        assert set(m.per_species_data) == {"species_a", "species_b", "species_c"}
        assert m.per_species_data["species_a"]["genes_in_orthogroups"] == 18000
        assert len(m.overall_data) == 1
        run_name = next(iter(m.overall_data))
        assert m.overall_data[run_name]["orthogroups"] == 9000

    def test_no_samples_raises(self, tmp_path):
        report.reset()
        report.analysis_files = [str(tmp_path)]
        report.search_files(["orthofinder"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()
