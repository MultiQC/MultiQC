"""Tests for the seqkit stats module"""

import pytest

from multiqc.modules.seqkit.stats import parse_stats_report


# Sample seqkit stats output with all columns (--all --tabular)
SAMPLE_STATS_ALL = """file	format	type	num_seqs	sum_len	min_len	avg_len	max_len	Q1	Q2	Q3	sum_gap	N50	N50_num	Q20(%)	Q30(%)	AvgQual	GC(%)	sum_n
sample1.fq.gz	FASTQ	DNA	1000000	100000000	100	100.0	100	100	100	100	0	100	1	98.5	86.2	25.64	31.85	1000
sample2.fq.gz	FASTQ	DNA	2000000	200000000	100	100.0	100	100	100	100	0	100	1	99.1	97.3	29.80	34.54	500
"""

# Sample seqkit stats output with basic columns only
SAMPLE_STATS_BASIC = """file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
reads.fastq	FASTQ	DNA	50000	5000000	80	100.0	150
assembly.fasta	FASTA	DNA	1000	10000000	500	10000.0	50000
"""

# Sample with Windows-style paths
SAMPLE_STATS_WINDOWS_PATH = """file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
C:\\data\\reads.fq.gz	FASTQ	DNA	1000	100000	100	100.0	100
"""

# Sample with stdin input
SAMPLE_STATS_STDIN = """file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
-	FASTQ	DNA	5000	500000	100	100.0	100
"""

# Empty file (just header)
SAMPLE_STATS_EMPTY = """file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
"""

# Invalid file (missing required columns)
SAMPLE_STATS_INVALID = """name	count	length
sample1	1000	100000
"""

# Space-separated format (default seqkit stats output without --tabular)
SAMPLE_STATS_SPACE_SEPARATED = """file                format  type  num_seqs   sum_len      min_len  avg_len  max_len
sample1.fq.gz       FASTQ   DNA   1000000    100000000    100      100.0    100
sample2.fq.gz       FASTQ   DNA   2000000    200000000    100      100.0    100
"""


class TestParseStatsReport:
    """Tests for parse_stats_report function"""

    def test_parse_all_columns(self):
        """Test parsing output with all columns from --all flag"""
        result = parse_stats_report(SAMPLE_STATS_ALL)

        assert len(result) == 2
        assert "sample1" in result
        assert "sample2" in result

        # Check sample1 data
        s1 = result["sample1"]
        assert s1["file"] == "sample1.fq.gz"
        assert s1["format"] == "FASTQ"
        assert s1["type"] == "DNA"
        assert s1["num_seqs"] == 1000000
        assert s1["sum_len"] == 100000000
        assert s1["min_len"] == 100
        assert s1["avg_len"] == 100.0
        assert s1["max_len"] == 100
        assert s1["Q1"] == 100
        assert s1["Q2"] == 100
        assert s1["Q3"] == 100
        assert s1["sum_gap"] == 0
        assert s1["N50"] == 100
        assert s1["N50_num"] == 1
        assert s1["Q20_pct"] == 98.5
        assert s1["Q30_pct"] == 86.2
        assert s1["AvgQual"] == 25.64
        assert s1["GC_pct"] == 31.85
        assert s1["sum_n"] == 1000

    def test_parse_basic_columns(self):
        """Test parsing output with basic columns only"""
        result = parse_stats_report(SAMPLE_STATS_BASIC)

        assert len(result) == 2
        assert "reads" in result
        assert "assembly" in result

        # Check reads data
        reads = result["reads"]
        assert reads["format"] == "FASTQ"
        assert reads["num_seqs"] == 50000
        assert reads["avg_len"] == 100.0

        # Check assembly data
        assembly = result["assembly"]
        assert assembly["format"] == "FASTA"
        assert assembly["num_seqs"] == 1000
        assert assembly["avg_len"] == 10000.0

        # These columns should not be present in basic output
        assert "Q20_pct" not in reads
        assert "N50" not in reads

    def test_parse_windows_path(self):
        """Test parsing handles Windows-style paths correctly"""
        result = parse_stats_report(SAMPLE_STATS_WINDOWS_PATH)

        assert len(result) == 1
        # os.path.basename should handle Windows paths
        assert "reads" in result

    def test_parse_stdin_input(self):
        """Test parsing handles stdin input (file = '-')"""
        result = parse_stats_report(SAMPLE_STATS_STDIN)

        assert len(result) == 1
        assert "stdin" in result
        assert result["stdin"]["num_seqs"] == 5000

    def test_parse_empty_file(self):
        """Test parsing empty file returns empty dict"""
        result = parse_stats_report(SAMPLE_STATS_EMPTY)
        assert result == {}

    def test_parse_invalid_format(self):
        """Test parsing invalid format returns empty dict"""
        result = parse_stats_report(SAMPLE_STATS_INVALID)
        assert result == {}

    def test_parse_single_line(self):
        """Test parsing file with only header line"""
        result = parse_stats_report("file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len")
        assert result == {}

    def test_sample_name_extension_stripping(self):
        """Test that various sequence file extensions are properly stripped"""
        test_cases = [
            ("sample.fq.gz", "sample"),
            ("sample.fastq.gz", "sample"),
            ("sample.fq", "sample"),
            ("sample.fastq", "sample"),
            ("sample.fa.gz", "sample"),
            ("sample.fasta.gz", "sample"),
            ("sample.fa", "sample"),
            ("sample.fasta", "sample"),
            ("sample.txt", "sample.txt"),  # Unknown extension not stripped
        ]

        for filename, expected_name in test_cases:
            data = f"file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n{filename}\tFASTQ\tDNA\t100\t10000\t100\t100.0\t100"
            result = parse_stats_report(data)
            assert expected_name in result, f"Expected {expected_name} for {filename}"

    def test_parse_space_separated(self):
        """Test parsing space-separated output (default seqkit stats format)"""
        result = parse_stats_report(SAMPLE_STATS_SPACE_SEPARATED)

        assert len(result) == 2
        assert "sample1" in result
        assert "sample2" in result

        # Check sample1 data
        s1 = result["sample1"]
        assert s1["format"] == "FASTQ"
        assert s1["type"] == "DNA"
        assert s1["num_seqs"] == 1000000
        assert s1["sum_len"] == 100000000
        assert s1["avg_len"] == 100.0
