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
                for k, v in parsed[1].items():
                    log.info("{}: {}".format(k, v))

        if len(d) > 0:
            # Allele-tagging report
            cats = OrderedDict()
            #cats['processed'] = dict(name='Total alignments processed')
            cats['skipped'] = dict(name='Unaligned reads')
            cats['unassignable'] = dict(name='Reads that could not be assigned to a genome')
            cats['genome1'] = dict(name='Reads assigned to genome 1')
            cats['genome2'] = dict(name='Reads assigned to genome 2')
            cats['conflicting'] = dict(name='Reads containing conflicting SNPs')
            cats['ambiguous'] = dict(name='Reads did not match either genome')
            # Subtract ambiguous from unassignable
            for k, v in d.items():
                if 'unassignable' and 'ambiguous' in v:
                    v['unassignable'] -= v['ambiguous']
            self.add_section(name="Allele-tagging report",
                             anchor="SNPsplit",
                             plot=bargraph.plot(d, cats))

            # SNP-coverage report

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
                match = re.match(r"Processed (\d+) read alignments in total", line)
                if match:
                    d["processed"] = int(match.group(1))
                    continue
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

                # SNP coverage report
                match = re.match(r"SNPs stored in total:\W+(\d+)", line) 
                if match:
                    d["totalStoredSNPs"] = int(match.group(1))
                    continue
                match = re.match(r"N-containing reads:\W+(\d+)", line) 
                if match:
                    d["Nreads"] = int(match.group(1))
                    continue
                match = re.match(r"non-N:\W+(\d+)", line) 
                if match:
                    d["nonNreads"] = int(match.group(1))
                    continue
                match = re.match(r"total:\W+(\d+)", line) 
                if match:
                    d["totalReads"] = int(match.group(1))
                    continue
                if line.startswith("Reads had a deletion of the N-masked position (and were thus dropped)"):
                    d["deletedNReads"] = int(line.split("\t")[1].split()[0])
                    continue
                match = re.match(r"Of which had multiple deletions of N-masked positions within the same read:	(\d+)", line)
                if match:
                    d["multiDeletedNReads"] = int(match.group(1))
                    continue

                # Of valid N-containing reads
                match = re.match(r"N was present in the list of known SNPs:\W+(\d+)", line)
                if match:
                    d["NIsKnownSNP"] = int(match.group(1))
                    continue
                # Generally not present
                match = re.match(r"Positions were skipped since they involved C>T SNPs:\W+(\d+)", line)
                if match:
                    d["skippedCTSNPs"] = int(match.group(1))
                    continue
                match = re.match(r"N was not present in the list of SNPs:\W+(\d+)", line)
                if match:
                    d["NIsUnknownSNP"] = int(match.group(1))
                    continue

                # Allele-specific paired-end sorting report
                match = re.match(r"Read pairs/singletons processed in total:\W+(\d+)", line)
                if match:
                    d["totalReadProcessed"] = int(match.group(1))
                    continue
                if line.startswith("Reads were unassignable (not overlapping SNPs):		"):
                    d["unassignableReads"] = int(line.split("\t")[2].split()[0])
                    continue
                match = re.match(r"Reads were specific for genome 1:\W+(\d+)", line)
                if match:
                    d["genome1Reads"] = int(match.group(1))
                    continue
                match = re.match(r"Reads were specific for genome 2:\W+(\d+)", line)
                if match:
                    d["genome2Reads"] = int(match.group(1))
                    continue
                match = re.match(r"Reads contained conflicting SNP information:\W+(\d+)", line)
                if match:
                    d["conflictingReads"] = int(match.group(1))
                    continue
        except:
            log.warning("{} was initially flagged as the output from SNPsplit, but that seems to not be the case. Skipping...".format(fname))
            return []

        # check the number of entries in d (see SNPsplit source code)
        log.info(len(d.keys()))
        l.append(d)

        return l
