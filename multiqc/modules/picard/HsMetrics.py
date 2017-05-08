#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard HsMetrics """

from collections import OrderedDict, defaultdict
import logging
import os
import re

from multiqc import config
from multiqc.plots import table, linegraph

# Initialise the logger
log = logging.getLogger(__name__)

DESC = {
        'BAIT_SET': 'The name of the bait set used in the hybrid selection.',
        'BAIT_TERRITORY': 'The number of bases which have one or more baits on top of them.',
        'BAIT_DESIGN_EFFICIENCY': 'Target terrirtoy / bait territory. 1 == perfectly efficient, 0.5 = half of baited bases are not target.',
        'TOTAL_READS': 'The total number of reads in the SAM or BAM file examine.',
        'PF_READS': 'The number of reads that pass the vendor\'s filter.',
        'PF_UNIQUE_READS': 'The number of PF reads that are not marked as duplicates.',
        'PCT_PF_READS': 'PF reads / total reads. The percent of reads passing filter.',
        'PCT_PF_UQ_READS': 'PF Unique Reads / Total Reads.',
        'PF_UQ_READS_ALIGNED': 'The number of PF unique reads that are aligned with mapping score > 0 to the reference genome.',
        'PCT_PF_UQ_READS_ALIGNED': 'PF Reads Aligned / PF Reads.',
        'PF_BASES_ALIGNED': 'The number of PF unique bases that are aligned with mapping score > 0 to the reference genome.',
        'PF_UQ_BASES_ALIGNED': 'The number of bases in the PF aligned reads that are mapped to a reference base. Accounts for clipping and gaps.',
        'ON_BAIT_BASES': 'The number of PF aligned bases that mapped to a baited region of the genome.',
        'NEAR_BAIT_BASES': 'The number of PF aligned bases that mapped to within a fixed interval of a baited region, but not on a baited region.',
        'OFF_BAIT_BASES': 'The number of PF aligned bases that mapped to neither on or near a bait.',
        'ON_TARGET_BASES': 'The number of PF aligned bases that mapped to a targeted region of the genome.',
        'PCT_SELECTED_BASES': 'On+Near Bait Bases / PF Bases Aligned.',
        'PCT_OFF_BAIT': 'The percentage of aligned PF bases that mapped neither on or near a bait.',
        'ON_BAIT_VS_SELECTED': 'The percentage of on+near bait bases that are on as opposed to near.',
        'MEAN_BAIT_COVERAGE': 'The mean coverage of all baits in the experiment.',
        'MEAN_TARGET_COVERAGE': 'The mean coverage of targets.',
        'MEDIAN_TARGET_COVERAGE': 'The median coverage of targets.',
        'PCT_USABLE_BASES_ON_BAIT': 'The number of aligned, de-duped, on-bait bases out of the PF bases available.',
        'PCT_USABLE_BASES_ON_TARGET': 'The number of aligned, de-duped, on-target bases out of the PF bases available.',
        'FOLD_ENRICHMENT': 'The fold by which the baited region has been amplified above genomic background.',
        'ZERO_CVG_TARGETS_PCT': 'The fraction of targets that did not reach coverage=1 over any base.',
        'PCT_EXC_DUPE': 'The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.',
        'PCT_EXC_MAPQ': 'The fraction of aligned bases that were filtered out because they were in reads with low mapping quality.',
        'PCT_EXC_BASEQ': 'The fraction of aligned bases that were filtered out because they were of low base quality.',
        'PCT_EXC_OVERLAP': 'The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.',
        'PCT_EXC_OFF_TARGET': 'The fraction of aligned bases that were filtered out because they did not align over a target base.',
        'OLD_80_BASE_PENALTY': 'The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.',
        'AT_DROPOUT': 'A measure of how undercovered <= 50% GC regions are relative to the mean.',
        'GC_DROPOUT': 'A measure of how undercovered >= 50% GC regions are relative to the mean.',
        'HET_SNP_SENSITIVITY': 'The theoretical HET SNP sensitivity.',
        'HET_SNP_Q': 'The Phred Scaled Q Score of the theoretical HET SNP sensitivity.',
        'FOLD_80_BASE_PENALTY': 'The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.',
        'TARGET_TERRITORY': 'The unique number of target bases in the experiment where target is usually exons etc.'
        }


def parse_reports(self):
    """ Find Picard HsMetrics reports and parse their data """

    # Set up vars
    self.picard_HsMetrics_data = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files('picard/hsmetrics', filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None
        for l in f['f']:
            # New log starting
            if 'picard.analysis.directed.CalculateHsMetrics' in l or \
               'picard.analysis.directed.CollectHsMetrics' in l and 'INPUT' in l:
                s_name = None
                keys = None

                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
                    parsed_data[s_name] = dict()

            if s_name is not None:
                if 'picard.analysis.directed.HsMetrics' in l and '## METRICS CLASS' in l:
                    keys = f['f'].readline().strip("\n").split("\t")
                elif keys:
                    vals = l.strip("\n").split("\t")
                    if len(vals) == len(keys):
                        j = 'NA'
                        if keys[0] == 'BAIT_SET':
                            j = vals[0]
                        parsed_data[s_name][j] = dict()
                        for i, k in enumerate(keys):
                            try:
                                parsed_data[s_name][j][k] = float(vals[i])
                            except ValueError:
                                parsed_data[s_name][j][k] = vals[i]
                    else:
                        s_name = None
                        keys = None

        # Remove empty dictionaries
        for s_name in list(parsed_data.keys()):
            for j in parsed_data[s_name].keys():
                if len(parsed_data[s_name][j]) == 0:
                    parsed_data[s_name].pop(j, None)
            if len(parsed_data[s_name]) == 0:
                parsed_data.pop(s_name, None)

        # Manipulate sample names if multiple baits found
        for s_name in parsed_data.keys():
            for j in parsed_data[s_name].keys():
                this_s_name = s_name
                if(len(parsed_data[s_name]) > 1):
                    this_s_name = "{}: {}".format(s_name, j)
                if this_s_name in self.picard_HsMetrics_data:
                    log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], this_s_name))
                self.add_data_source(f, this_s_name, section='HsMetrics')
                self.picard_HsMetrics_data[this_s_name] = parsed_data[s_name][j]


    # Filter to strip out ignored sample names
    self.picard_HsMetrics_data = self.ignore_samples(self.picard_HsMetrics_data)

    if len(self.picard_HsMetrics_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_HsMetrics_data, 'multiqc_picard_HsMetrics')

        # Add to general stats table
        # Swap question marks with -1
        data = self.picard_HsMetrics_data
        for s_name in data:
            if data[s_name]['FOLD_ENRICHMENT'] == '?':
                data[s_name]['FOLD_ENRICHMENT'] = -1

        self.general_stats_headers['FOLD_ENRICHMENT'] = {
            'title': 'Fold Enrichment',
            'min': 0,
            'format': '{:,.0f}',
            'scale': 'Blues',
        }
        try:
            covs = config.picard_config['general_stats_target_coverage']
            assert type(covs) == list
            assert len(covs) > 0
            covs = [str(i) for i in covs]
            log.debug("Custom Picard coverage thresholds: {}".format(", ".join([i for i in covs])))
        except (AttributeError, TypeError, AssertionError):
            covs = ['30']
        for c in covs:
            self.general_stats_headers['PCT_TARGET_BASES_{}X'.format(c)] = {
                'id': 'picard_target_bases_{}X'.format(c),
                'title': 'Target Bases {}X'.format(c),
                'description': 'Percent of target bases with coverage &ge; {}X'.format(c),
                'max': 100,
                'min': 0,
                'suffix': '%',
                'format': '{:,.0f}',
                'scale': 'RdYlGn',
                'modify': lambda x: self.multiply_hundred(x)
            }
        for s_name in data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( data[s_name] )
        data_table = _clean_table(data)
        self.add_section (
                name = 'HSMetrics',
                anchor = 'picard_hsmetrics',
                plot = table.plot(data_table, _get_headers(data_table))
        )
        tbases = _add_target_bases(data)
        self.add_section (
            name = tbases['name'],
            anchor = tbases['anchor'],
            description = tbases['description'],
            plot = tbases['plot']
        )
        hs_pen = _add_hs_penalty(data)
        if hs_pen is not None:
            self.add_section (
                name = hs_pen['name'],
                anchor = hs_pen['anchor'],
                description = hs_pen['description'],
                plot = hs_pen['plot']
            )

    # Return the number of detected samples to the parent module
    return len(self.picard_HsMetrics_data)

def _clean_table(data):
    data_table = defaultdict(dict)
    for s in data:
        for h in data[s]:
            if h.startswith("PCT"):
                continue
            if h.startswith("HS"):
                continue
            if h == "GENOME_SIZE":
                continue
            data_table[s][h] = data[s][h]
    return data_table

def _get_headers(data):
    header = {}
    for s in data:
        for h in data[s]:
            try:
                float(data[s][h])
            except:
                continue
            if h not in header:
                this = {
                'title': h.replace("_", " ").lower().capitalize(),
                'description' : DESC[h] if h in DESC else None,
                'scale': 'RdYlGn',
                'min': 0,
                'namespace': 'HsMetrics'
                }
                if h.find("READS") > -1:
                    this.update({
                        'shared_key': "read_count",
                        'modify': lambda x: x * config.read_count_multiplier
                        })
                    this['title'] = "{} {}".format(config.read_count_prefix, this['title'])
                elif h.find("BASES") > -1:
                    this.update({'shared_key': 'bases_count',
                                 'modify': lambda x: x / 1000000.0})
                    this['title'] = "Mb {0}".format(this['title'])
                if h in ["BAIT_TERRITORY", "TOTAL_READS", "TARGET_TERRITORY", "AT_DROPOUT", "GC_DROPOUT"]:
                    this.update({'hidden': True})
                header.update({h : this})
    return OrderedDict(sorted(header.items(), key=lambda t: t[1]['title']))

def _add_target_bases(data):
    subtitle = "The percentage of all target bases with at least <code>x</code> fold coverage."
    data_clean = defaultdict(dict)
    for s in data:
        for h in data[s]:
            if h.startswith("PCT_TARGET"):
                data_clean[s][int(h.replace("PCT_TARGET_BASES_", "")[:-1])] = data[s][h] * 100.0

    pconfig = { 'id': 'picard_percentage_target_bases',
                'title': 'Percentage of target bases',
                'xlab': 'Fold Coverage',
                'ylab': 'Pct of bases',
                'ymax': 100,
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',}
    return {
        'name': 'Target Region Coverage',
        'anchor': 'picard_hsmetrics_target_bases',
        'description': subtitle,
        'plot' : linegraph.plot(data_clean, pconfig)
    }

def _add_hs_penalty(data):
    subtitle = "The \"hybrid selection penalty\" incurred to get 80% of target bases to a given coverage. Can be used with the formula <code>required_aligned_bases = bait_size_bp * desired_coverage * hs_penalty</code>."
    data_clean = defaultdict(dict)
    any_non_zero = False
    for s in data:
        for h in data[s]:
            if h.startswith("HS_PENALTY"):
                data_clean[s][(h.replace("HS_PENALTY_", " ")[:-1])] = data[s][h]
                if data[s][h] > 0:
                    any_non_zero = True

    pconfig = { 'id': 'picard_hybrid_selection_penalty',
                'title': 'Hybrid Selection Penalty',
                'xlab': 'Fold Coverage',
                'ylab': 'Pct of bases',
                'ymax': 100,
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',}

    if any_non_zero:
        return {
            'name': 'HS penalty',
            'anchor': 'picard_hsmetrics_hs_penalty',
            'description': subtitle,
            'plot': linegraph.plot(data_clean, pconfig)
        }
