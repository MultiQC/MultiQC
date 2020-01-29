#!/usr/bin/env python

""" MultiQC module to do indexhopping checking """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import json
import os
import os.path
import itertools
import pprint

from multiqc.modules.base_module import BaseMultiqcModule

from pkg_resources import get_distribution
from multiqc.utils import report, util_functions, config
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger('multiqc')
#__version__ = get_distribution('indexhopping').version

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Check that these hooks haven't been disabled in the config file
        if getattr(config, 'disable_indexhopping', False) is True:
            log.debug("Skipping indexhopping as specified in config file")
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name = 'Indexhopping',
            anchor = 'indexhopping',
            info = 'PF clusters showing unexpected combinations of barcodes used in SampleSheet. Number of PF clusters and percentage of all PF clusters for a given lane.'
        )
        self.f_count = 0

        for f in self.find_log_files('bcl2fastq'):
            self.single_run(os.path.join(f['root'],f['fn']))


    def get_all_tags(self, stats):
        ret = dict()
        skip_re = re.compile('^N+$')
        for lane_stats in stats['ConversionResults']:
            ln = lane_stats['LaneNumber']
            ret[ln] = dict()
            ret[ln]['total'] = lane_stats['TotalClustersPF']
            ret[ln]['left'] = dict()
            ret[ln]['right'] = dict()
            ret[ln]['problem'] = dict()
            ret[ln]['problem']['left'] = []
            ret[ln]['problem']['right'] = []
            g = {0:'left',1:'right'}
            for demux in lane_stats['DemuxResults']:
                # skip samples with no indexes
                if not 'IndexMetrics' in demux:
                    continue
                for c in demux['IndexMetrics']:
                    index_seq = c['IndexSequence']
                    ind_tags = []
                    # skip single indexes
                    if index_seq.find('+') is -1:
                        continue
                    ind_tags.extend(index_seq.split('+'))
                    skip = 0
                    for i in 0,1:
                        if skip_re.match(ind_tags[i]):
                            skip = 1
                    if skip:
                        continue
                    for i in 0,1:
                        if ind_tags[i] in ret[ln][g[i]]:
                            # problem
                            ret[ln][g[i]][ind_tags[i]].append(demux['SampleId'])
                            if ind_tags[i] not in ret[ln]['problem'][g[i]]:
                                ret[ln]['problem'][g[i]].append(ind_tags[i])
                        else:
                            ret[ln][g[i]][ind_tags[i]] = [ demux['SampleId'] ]
        return ret


    def dual_index_lanes(self, run):
        ret = dict()
        for info in run['ReadInfosForLanes']:
            lane = info['LaneNumber']
            cnt = 0
            for rinfo in info['ReadInfos']:
                # some samples have index true with 0 cycles
                if rinfo['IsIndexedRead'] and rinfo['NumCycles'] > 0:
                    cnt += 1
            if cnt > 1:
                ret[lane] = True
                log.info("Dual:" + str(lane))
        return ret


    def has_same_tag(self, f, s):
        for i in f:
            if i in s:
                return 1
        return 0

    
    def get_sample_cluster_count(self,data):
        ret = dict()
        for counts in data['ConversionResults']:
            for dm in counts['DemuxResults']:
                sample = dm['SampleName']
                ret[sample] = dict()
                ret[sample]['clusters'] = dm['NumberReads']
                ret[sample]['indexhopping'] = 0
        return ret


    def single_run(self, f):
        with open(f) as json_file:
            try:
                run = json.load(json_file)
            except ValueError:
                log.warn('Could not parse file as json: {}'.format(f))
                return
        #try:
            #run = json.load(file(f))
        #except ValueError:
            #log.warn('Could not parse file as json: {}'.format(f))
            #return
        # first go through all tags and put them to a sample
        # then check through all missed combinations if match two different samples, hit!
        tags = self.get_all_tags(run)
        table_data = dict()
        lanes = dict()
        dual_index = self.dual_index_lanes(run)
        file_data = OrderedDict()
        self.perc = dict()
        self.f_count += 1
        run_id = run['RunId']
        file_data['file'] = f
        sample_count = self.get_sample_cluster_count(run)
        for ub in run['UnknownBarcodes']:
            lane = ub['Lane']
            if not lane in dual_index.keys():
                continue
            lanes[lane] = 0;
            ln = 'Lane ' + str(lane)
            file_data[ln] = []
            # for some reason there can be samples on lanes with no index, so skip lanes with no index
            if not tags[lane]['left']:
                continue
            for i in sorted(ub['Barcodes'], key=ub['Barcodes'].get, reverse=True):
                unmatch = i.split('+')
                tag1 = tag2 = None
                if unmatch[0] in tags[lane]['left']:
                    tag1 = tags[lane]['left'][unmatch[0]]
                if unmatch[1] in tags[lane]['right']:
                    tag2 = tags[lane]['right'][unmatch[1]]
                if tag1 is None or tag2 is None:
                    continue
                if self.has_same_tag(tag1,tag2):
                    continue
                lanes[lane] += ub['Barcodes'][i]
                for j in itertools.chain(tag1, tag2):
                    sample_count[j]['indexhopping'] += ub['Barcodes'][i]
                t1 = ",".join(tag1)
                t2 = ",".join(tag2)
                file_data[ln].append([ t1 + ' + ' + t2, i, ub['Barcodes'][i] ])
            self.perc[lane] = dict()
            self.perc[lane]['perc'] = float(lanes[lane]) / tags[lane]['total'] * 100
            run_id_lane = run_id + ' - ' + str(lane)
            table_data[run_id_lane] = dict()
            table_data[run_id_lane]['count'] = lanes[lane]
            table_data[run_id_lane]['perc'] = float(lanes[lane]) / tags[lane]['total'] * 100

        # Add section for report
        if not table_data:
            self.add_section (
                #name = 'Indexhopping',
                anchor = 'indexhopping-skip',
                description = 'No dual index samples found, cannot calculate indexhopping.'
            )
            return

        # warn if there is elements in tags['problem']['left'] or tags['problem']['right']
        desc_lanestats = 'Statistics about each lane' + self.get_problems_value(tags)
        log.info(desc_lanestats)
        self.add_section (
            name = 'Lane Statistics',
            anchor = 'indexhopping-lanestats',
            description = desc_lanestats,
            plot = self.lane_stats_table(table_data)
        )
        self.write_data_file(file_data, 'multiqc_indexhopping_' + str(self.f_count), data_format='json')
        cats = OrderedDict()
        cats["perc"] = { 'name': 'Percentage', 'color' : '#6CCAB0' } # green 7EDDA7
        self.add_section (
            name = 'Indexhopping ' + f,
            anchor = 'indexhopping-lane',
            #description = 'Indexhopping shows how many unmatched index pairs are in multiple samples.',
            plot = bargraph.plot(
                self.perc,
                cats,
                {
                    'id': 'indexhopping-lane',
                    'title': 'Indexhopping per lane %',
                    'tt_decimals': 5,
                    'cpswitch': False,
                    'tt_percentages': False
                }
            )
        )

        self.add_section (
            anchor = 'indexhopping-sample',
            plot = bargraph.plot(
                sample_count,
                { 'clusters':{ 'name':'Total clusters' }, 'indexhopping':{ 'name':'Indexhopping' } },
                {
                    'id': 'indexhopping-sample',
                    'cpswitch_c_active': False,
                    'title': 'Indexhopping per sample'
                }
            )
        )


    def get_problems_value(self, tag):
        # Lane x has same tags in left index on samples
        ret = '.'
        add = ''
        for lane in tag:
            if not tag[lane]['problem']['left'] and not tag[lane]['problem']['right']:
                continue
            if tag[lane]['problem']['left']:
                add += ' Lane ' + str(lane) + ' has multiple samples with same tag in i5 index: ' + ', '.join(tag[lane]['problem']['left']) + '.'
            if tag[lane]['problem']['right']:
                add += ' Lane ' + str(lane) + ' has multiple samples with same tag in i7 index: ' + ', '.join(tag[lane]['problem']['right']) + '.'
        if add:
            ret += "\n\nPROBLEM!" + add
        return ret

    def lane_stats_table(self, table_data):
        headers = OrderedDict()
        headers['count'] = {
            'title': 'Indexhopping PF Clusters',
            'description': 'indexhopping PF clusters.',
            'scale' : False,
            'format': '{:,.0f}',
            #'shared_key': 'yield'
        }
        headers['perc'] = {
            'title': 'Indexhopping %',
            'description': 'Percentage of indexhopping PF clusters.',
            'scale' : False,
            'format': '{:,.5f}',
            #'shared_key': 'yield'
        }
        table_config = {
            'namespace': 'indexhopping',
            'id': 'indexhopping-lane-stats-table',
            'table_title': 'Indexhopping Lane Statistics',
            'col1_header': 'Run ID - Lane',
            'no_beeswarm': True
        }
        return table.plot(table_data, headers, table_config)
