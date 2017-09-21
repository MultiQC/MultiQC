#!/usr/bin/env python
"""MultiQC module to parse the output from SNPsplit"""
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(name='SNPsplit',
                                            anchor='SNPsplit',
                                            target='SNPsplit',
                                            href='https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/',
                                            info='A tool to determine allele-specific alignments from high-throughput sequencing experiments that have been aligned to N-masked genomes')

        n = self.parse_SNPsplit()
        if n > 0:
            log.info("Found {} SNPsplit outputs".format(n))
        else:
            log.debug("Could not find any SNPsplit outputs in {}".format(config.analysis_dir))

    def parse_SNPsplit(self):
        d = dict()
        for f in self.find_log_files('SNPsplit'):
            parsed = self.parseSNPsplitOutput(f['f'], f['fn'])
            if len(parsed) == 2:
                if parsed[0] in d:
                    log.warning("Replacing duplicate sample {}".format(parsed[0]))
                d[parsed[0]] = parsed[1]
                self.add_data_source(f, section='SNPsplit')

        if len(d) > 0:
            # Allele-tagging report, the order matches that in the text file
            cats = OrderedDict()
            cats['skipped'] = dict(name='Unaligned reads')
            cats['unassignable'] = dict(name='Reads that could not be assigned to a genome')
            cats['genome1'] = dict(name='Reads assigned to genome 1')
            cats['genome2'] = dict(name='Reads assigned to genome 2')
            # Subtract ambiguous and C->T from unassignable
            CT = False
            for k, v in d.items():
                if 'unassignable' and 'ambiguous' in v:
                    v['unassignable'] -= v['ambiguous']
                if 'unassignable' and 'unassignableCT' in v:
                    v['unassignable'] -= v['unassignableCT']
                    CT = True
            if CT:
                cats['unassignableCT'] = dict(name='Reads unassigned due to a C->T SNP')
            cats['ambiguous'] = dict(name='Reads did not match either genome')
            cats['conflictingSNPs'] = dict(name='Reads containing conflicting SNPs')

            self.add_section(name="Allele-tagging report",
                             description="Per-sample metrics of how many reads were assigned to each genome.",
                             helptext="Possible reasons for not assigning reads to a genome are that they aren't aligned ('Unaligned reads'), they don't overlap a SNP ('Reads that could not be assigned to a genome'), they overlap informative SNPs, but don't contain the expected nucleotide for either genome ('Reads did not match either genome'), or the reads overlapped multiple informative SNPs, but there was a conflict in support for assignment to one genome over the other between the SNPs ('Reads contain conflicting SNPs'). Additionally, for bisulfite sequencing data reads overlapping C->T SNPs may be unassigned ('Reads unassigned due to a C->T SNP').",
                             anchor="SNPsplit",
                             plot=bargraph.plot(d, cats))

            # Allele-specific sorting report
            cats = OrderedDict()
            cats['unassignableReads'] = dict(name='Reads that could not be assigned to a genome')
            cats['genome1Reads'] = dict(name='Reads assigned to genome 1')
            cats['genome2Reads'] = dict(name='Reads assigned to genome 2')
            cats['G1UA'] = dict(name='PE Reads with one mate assigned to genome 1 and the other unassignable')
            cats['G2UA'] = dict(name='PE Reads with one mate assigned to genome 2 and the other unassignable')
            cats['G1G2'] = dict(name='PE Reads with mates assigned to different genomes')
            cats['conflictingReads'] = dict(name='Reads containing conflicting SNPs')
            self.add_section(name="Allele-specific sorting",
                             description="Per-sample metrics of how reads and pairs of reads were sorted into each genome.",
                             helptext="Note that metrics here may differ from those in the allele-tagging report. This occurs when paired-end reads are used, since 'tagging' only one read in a pair as arising from one genome can suffice in both reads being sorted there. Similarly, if two reads in a pair are tagged as arising from different genomes then the pair becomes unassignable.",
                             anchor="SNPsplit",
                             plot=bargraph.plot(d, cats))

        return len(d)

    def parseSNPsplitOutput(self, f, fname):
        l = []
        d = dict()
        first = True

        try:
            for line in f.splitlines():
                # Parse the sample name
                match = re.match(r"Input file:\W+'([^']+)'", line)
                if match:
                    if len(l) == 0:
                        l.append(match.group(1))
                    continue

                # Allele-tagging report
                match = re.match(r"Reads were unaligned and hence skipped: (\d+)", line)
                if match:
                    d["skipped"] = int(match.group(1))
                    continue
                match = re.match(r"(\d+) reads were unassignable", line)
                if match:
                    d["unassignable"] = int(match.group(1))
                    continue
                match = re.match(r"(\d+) reads were specific for genome 1", line)
                if match:
                    d["genome1"] = int(match.group(1))
                    continue
                match = re.match(r"(\d+) reads were specific for genome 2", line)
                if match:
                    d["genome2"] = int(match.group(1))
                    continue
                # Generally not present
                match = re.match(r"(\d+) reads that were unassignable contained C>T SNPs", line)
                if match:
                    d["unassignableCT"] = int(match.group(1))
                    continue
                match = re.match(r"(\d+) reads did not contain one of the expected bases", line)
                if match:
                    d["ambiguous"] = int(match.group(1))
                    continue
                match = re.match(r"(\d+) contained conflicting allele-specific SNPs", line)
                if match:
                    d["conflictingSNPs"] = int(match.group(1))
                    continue

                # Allele-specific sorting report
                if line.startswith("Reads were unassignable"):
                    d["unassignableReads"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Read pairs were unassignable (UA/UA):"):
                    d["unassignableReads"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Reads were specific for genome 1"):
                    d["genome1Reads"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Read pairs were specific for genome 1 (G1/G1)"):
                    d["genome1Reads"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Reads were specific for genome 2"):
                    d["genome2Reads"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Read pairs were specific for genome 2 (G2/G2)"):
                    d["genome2Reads"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Read pairs were a mix of G1 and UA"):
                    d["G1UA"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Read pairs were a mix of G2 and UA"):
                    d["G2UA"] = int(line.split("\t")[-1].split()[0])
                    continue
                if line.startswith("Read pairs were a mix of G1 and G2"):
                    d["G1G2"] = int(line.split("\t")[-1].split()[0])
                    continue
                match = re.match(r"Reads contained conflicting SNP information:\W+(\d+)", line)
                if match:
                    d["conflictingReads"] = int(match.group(1))
                    continue
                if line.startswith("Read pairs contained conflicting SNP information"):
                    d["conflictingReads"] = int(line.split("\t")[-1].split()[0])
                    continue
        except:
            log.warning("{} was initially flagged as the output from SNPsplit, but that seems to not be the case. Skipping...".format(fname))
            return []
        l.append(d)

        return l
