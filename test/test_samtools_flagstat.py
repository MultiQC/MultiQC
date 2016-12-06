#!/usr/bin/env python3

import unittest
import sys, os, math

# This line allows the tests to run if you just naively run this script.
# But the preferred way is to use run_tests.sh
sys.path.insert(0,'.')

from multiqc.modules.samtools.flagstat import parse_single_report

# From samtools 1.3
rep1 = """
==> small.bam.flagstat <==
5414 + 0 in total (QC-passed reads + QC-failed reads)
13 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
5350 + 0 mapped (98.82% : N/A)
5401 + 0 paired in sequencing
2709 + 0 read1
2692 + 0 read2
5011 + 0 properly paired (92.78% : N/A)
5273 + 0 with itself and mate mapped
64 + 0 singletons (1.18% : N/A)
206 + 0 with mate mapped to a different chr
81 + 0 with mate mapped to a different chr (mapQ>=5)
"""

# Same BAM file in samools 1.2
rep2 = """
==> small.bam.flagstat2 <==
5414 + 0 in total (QC-passed reads + QC-failed reads)
13 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
5350 + 0 mapped (98.82%:-nan%)
5401 + 0 paired in sequencing
2709 + 0 read1
2692 + 0 read2
5011 + 0 properly paired (92.78%:-nan%)
5273 + 0 with itself and mate mapped
64 + 0 singletons (1.18%:-nan%)
206 + 0 with mate mapped to a different chr
81 + 0 with mate mapped to a different chr (mapQ>=5)
"""

class T(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None

    def test_rep1(self):
        """Test that parsing rep1 produces expected results
        """
        res1 = parse_single_report(rep1)

        #I expect 13 + 13 + 3 + 3 + 1 things reported in total
        self.assertEqual(len(res1), 13 + 13 + 3 + 3 + 1 )

        self.assertEqual( (res1['total_passed'], res1['total_failed']),
                                (5414,                 0) )

        self.assertEqual(res1['flagstat_total'], 5414)

        self.assertEqual(res1['mapped_passed_pct'], 98.82)

        #I expect mapped_failed_pct to be float('nan')
        self.assertTrue(math.isnan(res1['mapped_failed_pct']))

    def test_rep2(self):
        """I expect rep2 to give identical results to rep1.
        """
        res1 = parse_single_report(rep1)
        res2 = parse_single_report(rep2)

        # But since nan != nan we have to strip these out.
        nans = [ k for k, v in res1.items() if math.isnan(v) ]
        for k in nans:
            del(res1[k])
            del(res2[k])

        self.assertEqual(res1, res2)
