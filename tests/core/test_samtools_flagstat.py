#!/usr/bin/env python3

import unittest
import os
import math


def slurp_file(fname):
    with open(os.path.dirname(__file__) + "/../data/modules/samtools/flagstat/" + fname) as fh:
        return fh.read()


# From samtools 1.3
rep1 = slurp_file("small.samtools13.flagstat.log.txt")

# Same BAM file in samools 1.2
rep2 = slurp_file("small.samtools12.flagstat.log.txt")


class T(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def test_rep1(self):
        """Test that parsing rep1 produces expected results"""
        from multiqc.modules.samtools.flagstat import parse_single_report

        res1 = parse_single_report(rep1)

        # I expect 13 + 13 + 3 + 3 + 1 things reported in total
        self.assertEqual(len(res1), 13 + 13 + 3 + 3 + 1)

        self.assertEqual((res1["total_passed"], res1["total_failed"]), (5414, 0))

        self.assertEqual(res1["flagstat_total"], 5414)

        self.assertEqual(res1["mapped_passed_pct"], 98.82)

        # I expect mapped_failed_pct to be float('nan')
        self.assertTrue(math.isnan(res1["mapped_failed_pct"]))

    def test_rep2(self):
        """I expect rep2 to give identical results to rep1."""
        from multiqc.modules.samtools.flagstat import parse_single_report

        res1 = parse_single_report(rep1)
        res2 = parse_single_report(rep2)

        # But since nan != nan we have to strip these out.
        nans = [k for k, v in res1.items() if math.isnan(v)]
        for k in nans:
            del res1[k]
            del res2[k]

        self.assertEqual(res1, res2)


if __name__ == "__main__":
    unittest.main()
