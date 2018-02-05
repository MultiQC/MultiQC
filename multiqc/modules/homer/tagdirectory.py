#!/usr/bin/env python

""" MultiQC module to parse output from HOMER tagdirectory """

import logging
import math
import os
import re
from multiqc.plots import bargraph, linegraph, scatter, table
from collections import OrderedDict

# Initialise the logger
log = logging.getLogger(__name__)


class TagDirReportMixin():


    def homer_tagdirectory(self):
        """ Find HOMER tagdirectory logs and parse their data """
        self.parse_gc_content()
        self.parse_re_dist()
        self.parse_tagLength_dist()
        self.parse_tagInfo_data()
        self.parse_FreqDistribution_data()
        self.homer_stats_table_tagInfo()

        return sum([len(v) for v in self.tagdir_data.values()])


    def parse_gc_content(self):
        """parses and plots GC content and genome GC content files"""
        # Find and parse GC content:
        for f in self.find_log_files('homer/GCcontent', filehandles=True):
            # Get the s_name from the parent directory
            s_name = os.path.basename(f['root'])
            s_name = self.clean_s_name(s_name, f['root'])
            parsed_data = self.parse_twoCol_file(f)
            if parsed_data is not None:
                if s_name in self.tagdir_data['GCcontent']:
                    log.debug("Duplicate GCcontent sample log found! Overwriting: {}".format(s_name))

                self.add_data_source(f, s_name, section='GCcontent')
                self.tagdir_data['GCcontent'][s_name] = parsed_data

        ## get esimated genome content distribution:
        for f in self.find_log_files('homer/genomeGCcontent', filehandles=True):
            parsed_data = self.parse_twoCol_file(f)
            if parsed_data is not None:
                if s_name + "_genome" in self.tagdir_data['GCcontent']:
                    log.debug("Duplicate genome GCcontent sample log found! Overwriting: {}".format(s_name+ "_genome"))

                self.add_data_source(f, s_name + "_genome", section='GCcontent')
                self.tagdir_data['GCcontent'][s_name + "_genome"] = parsed_data

        self.tagdir_data['GCcontent'] = self.ignore_samples(self.tagdir_data['GCcontent'])

        if len(self.tagdir_data['GCcontent']) > 0:
            self.add_section (
                name = 'Per Sequence GC Content',
                anchor = 'homer_per_sequence_gc_content',
                description = 'This plot shows the distribution of GC content.',
                plot = self.GCcontent_plot()
            )


    def parse_re_dist(self):
        """parses and plots restriction distribution files"""
        # Find and parse homer restriction distribution reports
        for f in self.find_log_files('homer/RestrictionDistribution', filehandles=True):
            s_name = os.path.basename(f['root'])
            s_name = self.clean_s_name(s_name, f['root'])
            parsed_data = self.parse_restriction_dist(f)
            if parsed_data is not None:
                if s_name in self.tagdir_data['restriction']:
                    log.debug("Duplicate Restriction Distribution sample log found! Overwriting: {}".format(s_name))

                self.add_data_source(f, s_name, section='restriction')
                self.tagdir_data['restriction'][s_name] = parsed_data
                self.tagdir_data['restriction_norm'][s_name] = self.normalize(parsed_data)

        self.tagdir_data['restriction'] = self.ignore_samples(self.tagdir_data['restriction'])
        self.tagdir_data['restriction_norm'] = self.ignore_samples(self.tagdir_data['restriction_norm'])

        if len(self.tagdir_data['restriction']) > 0:
            self.add_section (
                name = 'Restriction Site Tag Dist',
                anchor = 'homer-restrictionDist',
                description = 'This plot shows the distribution of tags around restriction enzyme cut sites.',
                helptext = '''Hi-C data often involves the digestion of DNA using a restriction enzyme.
                    A good quality control for the experiment is the centering of reads around the
                    restriction enzyme cut site.''',
                plot = self.restriction_dist_chart()
            )


    def parse_tagLength_dist(self):
        """parses and plots tag length distribution files"""
        # Find and parse homer tag length distribution reports
        for f in self.find_log_files('homer/LengthDistribution', filehandles=True):
            s_name = os.path.basename(f['root'])
            s_name = self.clean_s_name(s_name, f['root'])
            parsed_data = self.parse_length_dist(f)
            if parsed_data is not None:
                if s_name in self.tagdir_data['length']:
                    log.debug("Duplicate Length Distribution sample log found! Overwriting: {}".format(s_name))
                self.add_data_source(f, s_name, section='length')
                self.tagdir_data['length'][s_name] = parsed_data

        self.tagdir_data['length'] = self.ignore_samples(self.tagdir_data['length'])
        if len(self.tagdir_data['length']) > 0:
            self.add_section (
                name = 'Tag Length Distribution',
                anchor = 'homer-tagLength',
                description = 'This plot shows the distribution of tag length.',
                helptext = 'This is a good quality control for tag length inputed into Homer.',
                plot = self.length_dist_chart()
            )


    def parse_tagInfo_data(self):
        """parses and plots taginfo files"""
        # Find and parse homer taginfo reports
        for f in self.find_log_files('homer/tagInfo', filehandles=True):
            s_name = os.path.basename(f['root'])
            s_name = self.clean_s_name(s_name, f['root'])
            parsed_data = self.parse_tag_info_chrs(f)
            if parsed_data is not None:
                if s_name in self.tagdir_data['taginfo_total']:
                    log.debug("Duplicate tag info sample log found! Overwriting: {}".format(s_name))
                self.add_data_source(f, s_name, section='taginfo')
                self.tagdir_data['taginfo_total'][s_name] = parsed_data[0]
                self.tagdir_data['taginfo_total_norm'][s_name] = self.normalize(parsed_data[0])
                self.tagdir_data['taginfo_uniq'][s_name] = parsed_data[1]
                self.tagdir_data['taginfo_uniq_norm'][s_name] = self.normalize(parsed_data[1])


        for f in self.find_log_files('homer/tagInfo', filehandles=True):
            s_name = os.path.basename(f['root'])
            s_name = self.clean_s_name(s_name, f['root'])
            ## collected tag_info data for general stats table and store under 'header'
            parsed_data = self.parse_tag_info(f)
            if parsed_data is not None:
                self.tagdir_data['header'][s_name] = parsed_data

        self.tagdir_data['taginfo_total'] = self.ignore_samples(self.tagdir_data['taginfo_total'])
        self.tagdir_data['taginfo_total_norm'] = self.ignore_samples(self.tagdir_data['taginfo_total_norm'])
        self.tagdir_data['taginfo_uniq'] = self.ignore_samples(self.tagdir_data['taginfo_uniq'])
        self.tagdir_data['taginfo_uniq_norm'] = self.ignore_samples(self.tagdir_data['taginfo_uniq_norm'])

        if len(self.tagdir_data['taginfo_total']) > 0:
            self.add_section (
                name = 'Chromosomal Coverage',
                anchor = 'homer-tagInfo',
                description = 'This plot shows the distribution of tags along chromosomes.',
                helptext = '''This is a good quality control for tag distribution and
                        could be a good indication of large duplications or deletions.''',
                plot = self.tag_info_chart()
            )


    def parse_FreqDistribution_data(self):
        """parses and plots taginfo files"""
        # Find and parse homer tag FreqDistribution_1000 reports
        for f in self.find_log_files('homer/FreqDistribution', filehandles=True):
            s_name = os.path.basename(f['root'])
            s_name = self.clean_s_name(s_name, f['root'])
            parsed_data = self.parse_FreqDist(f)
            if parsed_data is not None:
                if s_name in self.tagdir_data['FreqDistribution']:
                    log.debug("Duplicate Freq Distribution sample log found! Overwriting: {}".format(s_name))

                self.add_data_source(f, s_name, section='FreqDistribution')
                self.tagdir_data['FreqDistribution'][s_name] = parsed_data


        self.tagdir_data['FreqDistribution'] = self.ignore_samples(self.tagdir_data['FreqDistribution'])

        if len(self.tagdir_data['FreqDistribution']) > 0:
            self.add_section (
                name = 'Frequency Distribution',
                anchor = 'homer-FreqDistribution',
                description = 'This plot shows the distribution of distance between PE tags.',
                helptext = '''It is expected the the frequency of PE tags decays with
                    increasing distance between the PE tags. This plot gives an idea
                     of the proportion of short-range versus long-range interactions.''',
                plot = self.FreqDist_chart()
            )



    def homer_stats_table_tagInfo(self):
        """ Add core HOMER stats to the general stats table from tagInfo file"""

        if len(self.tagdir_data['header']) == 0:
            return None

        headers = OrderedDict()
        headers['UniqPositions'] = {
            'title': 'Uniq Pos',
            'description': 'Numer of Unique Di-Tags Passed Through HOMER',
            'format': '{:,.0f}',
            'modify': lambda x: x * 0.000001,
            'suffix': "M"

        }
        headers['TotalPositions'] = {
            'title': 'Total Pos',
            'description': 'Numer of Total Di-Tags Passed Through HOMER',
            'format': '{:,.0f}',
            'modify': lambda x: x * 0.000001,
            'suffix': "M"
        }
        headers['fragmentLengthEstimate'] = {
            'title': 'fragment Length',
            'description': 'Estimate of Fragnment Length',
            'format': '{:,.0f}'
        }
        headers['peakSizeEstimate'] = {
            'title': 'Peak Size',
            'description': 'Estimate of Peak Size',
            'format': '{:,.0f}'
        }
        headers['tagsPerBP'] = {
            'title': 'tagsPerBP',
            'description': 'average tags Per basepair',
            'format': '{:,.3f}',
        }
        headers['TagsPerPosition'] = {
            'title': 'averageTagsPerPosition',
            'description': 'Average Tags Per Position',
            'format': '{:,.2f}'
        }
        headers['averageTagLength'] = {
            'title': 'TagLength',
            'description': 'Average Tag Length',
            'format': '{:,.0f}'
        }
        headers['averageFragmentGCcontent'] = {
            'title': 'GCcontent',
            'description': 'Average Fragment GC content',
            'max': 1,
            'min': 0,
            'format': '{:,.2f}'
        }
        self.general_stats_addcols(self.tagdir_data['header'], headers, 'HOMER')


    def homer_stats_table_interChr(self):
        """ Add core HOMER stats to the general stats table from FrequencyDistribution file"""

        headers = OrderedDict()
        headers['InterChr'] = {
            'title': 'InterChr',
            'description': 'Fraction of Reads forming inter chromosomal interactions',
            'format': '{:,.4f}'
        }
        self.general_stats_addcols(self.tagdir_data['FreqDistribution'], headers, 'Homer-InterChr')


    def normalize(self, mydict, target=100):
       raw = sum(mydict.values())
       factor = target/raw
       return { key:value*factor for key,value in mydict.items() }


    def parse_twoCol_file(self, f):
        """ Parse HOMER tagdirectory GCcontent file. """
        parsed_data = dict()
        firstline = True
        for l in f['f']:
            if firstline:    #skip first line
                firstline = False
                continue
            s = l.split("\t")
            if len(s) > 1:
                k = float(s[0].strip())
                v = float(s[2].strip())
                parsed_data[k] = v

        return parsed_data



    def parse_restriction_dist(self, f):
        """ Parse HOMER tagdirectory petagRestrictionDistribution file. """
        parsed_data = dict()
        firstline = True
        for l in f['f']:
            if firstline:    #skip first line
                firstline = False
                continue
            s = l.split("\t")
            if len(s) > 1:
                nuc = float(s[0].strip())
                v1 = float(s[1].strip())
                v2 = float(s[2].strip())
                v = v1 + v2
                #parsed_data.update({nuc:v1})
                #parsed_data.update({nuc:v2})
                parsed_data.update({nuc:v})
        return parsed_data


    def parse_length_dist(self, f):
        """ Parse HOMER tagdirectory tagLengthDistribution file. """
        parsed_data = dict()
        firstline = True
        for l in f['f']:
            if firstline:    #skip first line
                firstline = False
                continue
            s = l.split("\t")
            if len(s) > 1:
                k = float(s[0].strip())
                v = float(s[1].strip())
                parsed_data[k] = v

        return parsed_data


    def parse_tag_info(self, f):
        """ Parse HOMER tagdirectory taginfo.txt file to extract statistics in the first 11 lines. """
        # General Stats Table
        tag_info = dict()
        for l in f['f']:
            s = l.split("=")
            if len(s) > 1:
                if s[0].strip() == 'genome':
                    ss = s[1].split("\t")
                    if len(ss) > 2:
                        tag_info['genome'] = ss[0].strip()
                        try:
                            tag_info['UniqPositions'] = float(ss[1].strip())
                            tag_info['TotalPositions'] = float(ss[2].strip())
                        except:
                            tag_info['UniqPositions'] = ss[1].strip()
                            tag_info['TotalPositions'] = ss[2].strip()
                try:
                    tag_info[s[0].strip()] = float(s[1].strip())
                except ValueError:
                    tag_info[s[0].strip()] = s[1].strip()
        return tag_info




    def parse_tag_info_chrs(self, f, convChr=True):
        """ Parse HOMER tagdirectory taginfo.txt file to extract chromosome coverage. """
        parsed_data_total = OrderedDict()
        parsed_data_uniq = OrderedDict()
        remove = ["hap", "random", "chrUn", "cmd", "EBV", "GL", "NT_"]
        for l in f['f']:
            s = l.split("\t")
            key = s[0].strip()
            # skip header
            if '=' in l or len(s) != 3:
                continue
            if convChr:
                if any(x in key for x in remove):
                    continue
            try:
                vT = float(s[1].strip())
                vU = float(s[2].strip())
            except ValueError:
                continue

            parsed_data_total[key] = vT
            parsed_data_uniq[key] = vU

        return [parsed_data_total, parsed_data_uniq]


    def parse_FreqDist(self, f):
        """ Parse HOMER tagdirectory petag.FreqDistribution_1000 file. """
        parsed_data = dict()
        firstline = True
        for l in f['f']:
            if firstline:
                firstline = False
                continue
            else:
                s = l.split("\t")
                if len(s) > 1:
                    k = s[0].strip()
                    if k.startswith("More than "):
                        k = re.sub("More than ", "", k)

                    k = float(k)
                    v = float(s[1].strip())
                    parsed_data[k] = v
        return parsed_data

    def parse_FreqDist_interChr(self, f):
        """ Parse HOMER tagdirectory petag.FreqDistribution_1000 file to get inter-chromosomal interactions. """
        parsed_data = dict()
        firstline = True
        for l in f['f']:
            if firstline:
                firstline = False
                interChr = float(re.sub("\)", "", l.split(":")[1]))
            else:
                break
        parsed_data['interChr'] = interChr
        return parsed_data


    def restriction_dist_chart (self):

        """ Make the petagRestrictionDistribution plot """

        pconfig = {
            'id': 'petagRestrictionDistribution',
            'title': 'Restriction Distribution',
            'ylab': 'Reads',
            'xlab': 'Distance from cut site (bp)',
            'data_labels': [
                        {'name': 'Number of Tags'},
                        {'name': 'Percenatge'}
                            ]
       }
        datasets = [
            self.tagdir_data['restriction'],
            self.tagdir_data['restriction_norm']
        ]

        return linegraph.plot(datasets, pconfig)




    def length_dist_chart (self):

        """ Make the tagLengthDistribution plot """

        pconfig = {
            'id': 'tagLengthDistribution',
            'cpswitch': True,
            'title': 'Tag Length Distribution',
            'ylab': 'Fraction of Tags',
            'xlab': 'Tag Length (bp)'
        }
        return linegraph.plot(self.tagdir_data['length'], pconfig)



    def GCcontent_plot (self):
        """ Create the HTML for the Homer GC content plot """

        pconfig = {
            'id': 'homer-tag-directory-gc-content',
            'title': 'Homer: Tag Directory Per Sequence GC Content',
            'smooth_points': 200,
            'smooth_points_sumcounts': False,
            'ylab': 'Normalized Count',
            'xlab': '% GC',
            'ymin': 0,
            'xmax': 1,
            'xmin': 0,
            'yDecimals': True,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}'
        }
        return linegraph.plot(self.tagdir_data['GCcontent'], pconfig)


    def tag_info_chart (self):

        """ Make the taginfo.txt plot """

        ## TODO: human chrs on hg19. How will this work with GRCh genome or other, non human, genomes?
        # nice if they are ordered by size
        ucsc = ["chr" + str(i) for i in range(1,23)].append([ "chrX", "chrY", "chrM"])
        ensembl = list(range(1,23)).append([ "X", "Y", "MT"])
        pconfig = {
            'id': 'tagInfo',
            'title': 'Homer: Tag Info Distribution',
            'ylab': 'Tags',
            'cpswitch_counts_label': 'Number of Tags'
        }

        ## check if chromosomes starts with "chr" (UCSC) or "#" (ensembl)
        sample1 = next(iter(self.tagdir_data['taginfo_total']))
        chrFormat = next(iter(self.tagdir_data['taginfo_total'][sample1]))

        if ("chr" in chrFormat):
            chrs = ucsc
        else:
            chrs = ensembl

        return bargraph.plot(self.tagdir_data['taginfo_total'], chrs, pconfig)


    def FreqDist_chart (self):
        """ Make the petag.FreqDistribution_1000 plot """
        # Take a log of the data before plotting so that we can
        # reduce the number of points to plot evenly
        pdata = {}
        for idx, s_name in enumerate(self.tagdir_data['FreqDistribution']):
            pdata[s_name] = {}
            for x, y in self.tagdir_data['FreqDistribution'][s_name].items():
                try:
                    pdata[s_name][math.log(float(x))] = y
                except ValueError:
                    pass
        pconfig = {
            'id': 'FreqDistribution',
            'title': 'Frequency Distribution',
            'ylab': 'Fraction of Reads',
            'xlab': 'Log10(Distance between regions)',
            'data_labels': ['Reads', 'Percent'],
            'smooth_points': 500,
            'smooth_points_sumcounts': False,
            'yLog' : True
        }
        return linegraph.plot(pdata, pconfig)
