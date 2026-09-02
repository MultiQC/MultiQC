import pytest

from multiqc.modules.cramino.cramino import (
    build_cramino_stats,
    filter_karyotype_contigs,
    parse_cramino_lines,
)

ALIGNED_WITH_KARYOTYPE = """File name\tsample.bam
Number of alignments\t1000
% from total alignments\t96.00
Number of reads\t950
Yield [Gb]\t1.50
Mean coverage\t5.00
Yield [Gb] (>25kb)\t0.20
N50\t600
N75\t500
Median length\t550.00
Mean length\t700.00

Median identity\t99.40
Mean identity\t99.00
Modal identity\t100.0



# Normalized read count per chromosome

chr1\t0.90
chrX\t0.80
chrUn_KI270333v1\t42.00

Path\t/work/sample.bam
Creation time\tNA
"""

UNALIGNED_UBAM = """File name\tsample.ubam.bam
Number of alignments\t500
% from total alignments\t100.00
Number of reads\t500
Yield [Gb]\t0.01
Mean coverage\t0.00
Yield [Gb] (>25kb)\t0.00
N50\t400
N75\t350
Median length\t380.00
Mean length\t390.50

Path\t/work/sample.ubam.bam
Creation time\tNA
"""

ESTIMATED_IDENTITY = """File name\tsample.bam
Number of alignments\t200
% from total alignments\t99.00
Number of reads\t198
Yield [Gb]\t0.02
Mean coverage\t0.10
Yield [Gb] (>25kb)\t0.00
N50\t300
N75\t250
Median length\t280.00
Mean length\t290.00

Median est. identity\t97.00
Mean est. identity\t96.50
Modal est. identity\t98.0

Path\t/work/sample.bam
Creation time\tNA
"""

NO_CONTIGS_FOUND = """File name\tsample.bam
Number of alignments\t10
% from total alignments\t100.00
Number of reads\t10
Yield [Gb]\t0.0001
Mean coverage\t0.00
Yield [Gb] (>25kb)\t0.00
N50\t100
N75\t90
Median length\t95.00
Mean length\t96.00



# Warning - no contigs found in BAM file!


Path\t/work/sample.bam
Creation time\tNA
"""

INVALID_FORMAT = """This is not cramino output
just some random text
"""


class TestParseCraminoLines:
    def test_aligned_with_karyotype(self):
        summary, karyotype = parse_cramino_lines(ALIGNED_WITH_KARYOTYPE.splitlines())
        assert summary["Number of alignments"] == "1000"
        assert summary["Median identity"] == "99.40"
        assert summary["Path"] == "/work/sample.bam"
        assert summary["Creation time"] == "NA"
        assert karyotype == {"chr1": 0.90, "chrX": 0.80, "chrUn_KI270333v1": 42.00}

    def test_unaligned_ubam_has_no_identity_or_karyotype(self):
        summary, karyotype = parse_cramino_lines(UNALIGNED_UBAM.splitlines())
        assert "Median identity" not in summary
        assert karyotype == {}

    def test_estimated_identity_labels(self):
        summary, karyotype = parse_cramino_lines(ESTIMATED_IDENTITY.splitlines())
        assert summary["Median est. identity"] == "97.00"
        assert "Median identity" not in summary

    def test_no_contigs_found_warning(self):
        summary, karyotype = parse_cramino_lines(NO_CONTIGS_FOUND.splitlines())
        assert karyotype == {}
        assert summary["Number of alignments"] == "10"

    def test_empty_file(self):
        summary, karyotype = parse_cramino_lines([])
        assert summary == {}
        assert karyotype == {}

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError):
            parse_cramino_lines(INVALID_FORMAT.splitlines())


class TestBuildCraminoStats:
    def test_full_aligned_stats(self):
        summary, _ = parse_cramino_lines(ALIGNED_WITH_KARYOTYPE.splitlines())
        stats = build_cramino_stats(summary)
        assert stats["num_alignments"] == 1000
        assert stats["num_reads"] == 950
        assert stats["yield_gb"] == 1.50
        assert stats["n50"] == 600
        assert stats["median_identity"] == 99.40
        assert stats["mean_identity"] == 99.00
        assert stats["modal_identity"] == 100.0

    def test_ubam_stats_have_no_identity_keys(self):
        summary, _ = parse_cramino_lines(UNALIGNED_UBAM.splitlines())
        stats = build_cramino_stats(summary)
        assert "median_identity" not in stats
        assert "mean_identity" not in stats
        assert "modal_identity" not in stats
        assert stats["num_reads"] == 500

    def test_estimated_identity_mapped_to_canonical_keys(self):
        summary, _ = parse_cramino_lines(ESTIMATED_IDENTITY.splitlines())
        stats = build_cramino_stats(summary)
        assert stats["median_identity"] == 97.00
        assert stats["mean_identity"] == 96.50
        assert stats["modal_identity"] == 98.0

    def test_missing_required_field_raises(self):
        with pytest.raises(KeyError):
            build_cramino_stats({"Number of alignments": "10"})


class TestFilterKaryotypeContigs:
    KARYOTYPE = {
        "chr1": 0.9,
        "chrX": 0.8,
        "chrM": 71.0,
        "chrUn_KI270333v1": 42.0,
        "chr1_KI270706v1_random": 1.1,
        "HLA-A*01:01:01:01": 2.0,
    }

    def test_default_no_filtering(self):
        result = filter_karyotype_contigs(self.KARYOTYPE, include_contigs=[], exclude_contigs=[])
        assert result == self.KARYOTYPE

    def test_exclude_patterns(self):
        result = filter_karyotype_contigs(
            self.KARYOTYPE,
            include_contigs=[],
            exclude_contigs=["chrM", "chrUn*", "*_random", "HLA*"],
        )
        assert result == {"chr1": 0.9, "chrX": 0.8}

    def test_include_patterns(self):
        result = filter_karyotype_contigs(self.KARYOTYPE, include_contigs=["chr*"], exclude_contigs=[])
        assert set(result.keys()) == {"chr1", "chrX", "chrM", "chr1_KI270706v1_random", "chrUn_KI270333v1"}

    def test_exclude_takes_precedence_over_include(self):
        result = filter_karyotype_contigs(
            self.KARYOTYPE,
            include_contigs=["chr*"],
            exclude_contigs=["chrM"],
        )
        assert "chrM" not in result
        assert "chr1" in result
