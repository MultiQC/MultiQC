HEADER = (
    "chrom\tpos\tref\tgenotype\tfilters\tA\tC\tG\tT\tN\tconsensus_read\tconsensus_insert\tconsensus_call"
    "\tconsensus_qual\ta\tc\tg\tt\tn\n"
)
ROW = "chr1\t{pos}\tA\tA/T\tPASS\t1\t0\t0\t1\t0\t{read}\tchr1:1-100\tT\t40\t1\t0\t0\t1\t0\n"
REVIEW = HEADER + ROW.format(pos=10, read="r1") + ROW.format(pos=10, read="r2") + ROW.format(pos=20, read="r1")


def test_review_summary(run_fgumi):
    module = run_fgumi({"S1.txt": REVIEW})
    assert module.saved_raw_data["multiqc_fgumi_review"]["S1"] == {"sites": 2, "consensus_reads": 2, "observations": 3}
