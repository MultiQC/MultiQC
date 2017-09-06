#!/usr/bin/env python

""" MultiQC module to parse output from HOMER tagdirectory """

import logging
import os
import re
from multiqc.plots import bargraph, linegraph, scatter, table
from collections import OrderedDict

# Initialise the logger
log = logging.getLogger(__name__)


class TagDirReportMixin():


    def homer_tagdirectory(self):
        """ Find HOMER tagdirectory logs and parse their data """



        # Find and parse GC content:
        for f in self.find_log_files('homer/GCcontent', filehandles=True):
            # Get the s_name from the parent directory
            s_name = os.path.basename(f['root'])
            s_name = self.clean_s_name(s_name, f['root'])
            parsed_data = self.parse_GCcontent(f)
            if parsed_data is not None:
                if s_name in self.tagdir_data['GCcontent']:
                    log.debug("Duplicate GCcontent sample log found! Overwriting: {}".format(s_name))

                self.add_data_source(f, s_name, section='GCcontent')
                self.tagdir_data['GCcontent'][s_name] = parsed_data

        ## get esimated genome content distribution:
        for f in self.find_log_files('homer/genomeGCcontent', filehandles=True):
            parsed_data = self.parse_GCcontent(f)
            if parsed_data is not None:
                if s_name + "_genome" in self.tagdir_data['GCcontent']:
                    log.debug("Duplicate genome GCcontent sample log found! Overwriting: {}".format(s_name+ "_genome"))

                self.add_data_source(f, s_name + "_genome", section='GCcontent')
                self.tagdir_data['GCcontent'][s_name + "_genome"] = parsed_data


        ## plot GC content:  
        self.tagdir_data = self.ignore_samples(self.tagdir_data) 
        description = '<p>This plot shows the distribution of GC content.</p>'
        helptext = '<p>This is a good quality control for GC bias</p>'

        self.add_section (
            name = 'Per Sequence GC Content',
            anchor = 'homer_per_sequence_gc_content',
            description = description,
            plot = self.GCcontent_plot()
        )



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
                
        self.tagdir_data = self.ignore_samples(self.tagdir_data)

        description = '<p>This plot shows the distribution of tags around restriction enzyme cut sites.</p>'
        helptext = '<p>Hi-C data often involves the digestion of DNA using a restriction enzyme. A good quality control for the experiment is the centering of reads around the restriction enzyme cut site.</p>'
        
        self.add_section (
            name = 'PE Tag Distribution Around Restriction Sites',
            anchor = 'homer-restrictionDist',
            description = description,
            helptext = helptext,
            plot = self.restriction_dist_chart()
        )
        
        

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
        
        self.tagdir_data = self.ignore_samples(self.tagdir_data)
        description = '<p>This plot shows the distribution of tag length.</p>'
        helptext = '<p>This is a good quality control for tag length inputed into Homer.</p>'

        self.add_section (
        name = 'Tag Length Distribution',
        anchor = 'homer-tagLength',
        description = description,
        helptext = helptext,
        plot = self.length_dist_chart()
        )
        

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


        self.tagdir_data = self.ignore_samples(self.tagdir_data)
        description = '<p>This plot shows the distribution of tags along chromosomes.</p>'
        helptext = '<p>This is a good quality control for tag distribution and could be a good indication of large duplications or deletions.</p>'

        self.add_section (
        name = 'Tag Info Chromosomal Coverage',
        anchor = 'homer-tagInfo',
        description = description,
        helptext = helptext,
        plot = self.tag_info_chart()
        )


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


        
        self.tagdir_data = self.ignore_samples(self.tagdir_data)
        description = '<p>This plot shows the distribution of distance between PE tags.</p>'
        helptext = '<p>It is expected the the frequency of PE tags decays with increasing distance between the PE tags. This plot gives an idea of the proportion of short-range versus long-range interactions.</p>'

        self.add_section (
        name = 'Frequency Distribution',
        anchor = 'homer-FreqDistribution',
        description = description,
        helptext = helptext,
        plot = self.FreqDist_chart()
        )


        # Basic Stats Table
        self.homer_stats_table_tagInfo()

        if len(self.tagdir_data) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.tagdir_data, 'multiqc_homer_tagdirectory')
        return len(self.tagdir_data)




    def homer_stats_table_tagInfo(self):
        """ Add core HOMER stats to the general stats table from tagInfo file"""

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


    def reduce_points_average(self, mydict, iterations = 1):
        """
        reduces the numbers of points in dictionary by averging every 2 consecutive points
        """
        ## should just reduce points by skipping some, or better to average?
        counter = 0
        myOrderedDict = OrderedDict(sorted(mydict.items()))
        keys = myOrderedDict.keys()
        new_keys = [(a + b) / 2.0 for a, b in zip(keys[::2], keys[1::2])]
        values = myOrderedDict.values()
        new_values = [(a + b) / 2.0 for a, b in zip(values[::2], values[1::2])]
        newDict = dict(zip(new_keys, new_values))
        counter = counter + 1
        
        while counter < iterations:
            self.reduce_points_average(newDict, iterations - 1)
            
        return newDict


    def reduce_points_skip(self, mydict, iterations = 1):
        """
        reduces the numbers of points in dictionary by skipping every other point
        """
        counter = 0
        myOrderedDict = OrderedDict(sorted(mydict.items()))
        keys = myOrderedDict.keys()
        new_keys = keys[::2]
        values = myOrderedDict.values()
        new_values = values[::2]
        newDict = dict(zip(new_keys, new_values))
        counter = counter + 1
        
        while counter < iterations:
            self.reduce_points_skip(newDict, iterations - 1)
            
        return newDict


    def normalize(self, mydict, target=100):
       raw = sum(mydict.values())
       factor = target/raw
       return {key:value*factor for key,value in mydict.iteritems()}


    def parse_GCcontent(self, f):
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
        counter = 0
        for l in f['f']:
            if counter == 1:  
                s = l.split("\t")
                tag_info['UniqPositions'] = float(s[1].strip())
                tag_info['TotalPositions'] = float(s[2].strip())
            if counter == 4:
                s = l.split("\t")
                tag_info['fragmentLengthEstimate'] = float(s[0].strip().split("=")[1])
            if counter == 5:
                s = l.split("\t")
                tag_info['peakSizeEstimate'] = float(s[0].strip().split("=")[1])
            if counter == 6:
                s = l.split("\t")
                tag_info['tagsPerBP'] = float(s[0].strip().split("=")[1])
            if counter == 7:
                s = l.split("\t")
                tag_info['averageTagsPerPosition'] = float(s[0].strip().split("=")[1])          
            if counter == 8:
                s = l.split("\t")
                tag_info['averageTagLength'] = float(s[0].strip().split("=")[1])
            if counter == 9:
                s = l.split("\t")
                tag_info['gsizeEstimate'] = float(s[0].strip().split("=")[1])
            if counter == 10:
                s = l.split("\t")
                tag_info['averageFragmentGCcontent'] = float(s[0].strip().split("=")[1])
            if counter == 11:
                break
            counter = counter + 1
        return tag_info




    def parse_tag_info_chrs(self, f, convChr=True):
        """ Parse HOMER tagdirectory taginfo.txt file to extract chromosome coverage. """
        parsed_data_total = dict()
        parsed_data_uniq = dict()
        remove = ["hap", "random", "chrUn", "cmd"]
        ## skip first 11 lines
        counter = 0
        for l in f['f']:
            ## skip first 11 lines
            if counter < 11:
                counter = counter + 1
                continue    
            s = l.split("\t")

            key = s[0].strip()

            if len(s) > 1:
                if convChr:
                    if any(x in key for x in remove):
                        continue

                vT = float(s[1].strip())
                vU = float(s[2].strip())

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
        datasets = [
            self.tagdir_data['length'],
        ]

        return linegraph.plot(datasets, pconfig)



    def GCcontent_plot (self):
        """ Create the HTML for the Homer GC content plot """

        pconfig = {
            'id': 'Homer Tag Directory GC Content',
            'title': 'Per Sequence GC Content',
            'ylab': 'Normalized Count',
            'xlab': '% GC',
            'ymin': 0,
            'xmax': 1,
            'xmin': 0,
            'yDecimals': True,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}'
        }

        datasets = [
            self.tagdir_data['GCcontent']
        ]


        return linegraph.plot(datasets, pconfig)


    def tag_info_chart (self):

        """ Make the taginfo.txt plot """

        ## TODO: human chrs on hg19. How will this work with GRCh genome or other, non human, genomes? nice if they are ordered by size
        chrs = ["chr" + str(i) for i in range(1,23)]
        chrs.extend([ "chrX", "chrY", "chrM"])

        pconfig = {
            'id': 'tagInfo',
            'title': 'Tag Info Distribution',
            'ylab': 'Tags',
            'cpswitch_counts_label': 'Number of Tags'
        }

        datasets = [
            self.tagdir_data['taginfo_total'],
        ]

        return bargraph.plot(datasets, chrs, pconfig)


    def FreqDist_chart (self):

        """ Make the petag.FreqDistribution_1000 plot """

        pconfig = {
            'id': 'FreqDistribution',
            'title': 'Frequency Distribution',
            'ylab': 'Log10(Fraction of Reads)',
            'xlab': 'Log10(Distance between regions)',
            'data_labels': ['Reads', 'Percent'],
            'smooth_points': 2000,    
            'smooth_points_sumcounts': False,
            ## TODO: xLog does seem to work? 
            'xLog' : True,
            'yLog' : True
       }
        datasets = [
            self.tagdir_data['FreqDistribution']
        ]

        return linegraph.plot(datasets, pconfig)

