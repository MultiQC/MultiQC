#!/usr/bin/env python

""" MultiQC module to parse output from gffcompare """

from __future__ import print_function
import logging
import os
from collections import OrderedDict

from multiqc.plots import scatter, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='GffCompare', anchor='gffcompare',
                href='https://ccb.jhu.edu/software/stringtie/gffcompare.shtml',
        info="is a tool to compare, merge and annotate one or more GFF files with a reference annotation in GFF format.")

        # Parse stats file
        # Everything is hardcoded with linenumbers, needs to be adjusted if output should change. 
        # Other versions of GffCompare with different ouptut formats might not work.
        self.gffcompare_data = {}
        for f in self.find_log_files('gffcompare'):
            #print(f)
            sample = f['s_name']
            self.gffcompare_data[sample]= {}
            lines = f['f'].splitlines()

            ## Transcript and loci numbers:
            # Query
            c = [int(s) for s in lines[5].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts'] = {}
            self.gffcompare_data[sample]['counts']['query'] = {}
            self.gffcompare_data[sample]['counts']['query']['mRNA'] =c[0]
            self.gffcompare_data[sample]['counts']['query']['loci'] = c[1]
            self.gffcompare_data[sample]['counts']['query']['multi_exon_transcripts'] = c[2]
            c = [int(s) for s in lines[6].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['query']['multi_transcript_loci'] = c[0]
            
            # Reference
            c = [int(s) for s in lines[7].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['reference'] = {}
            self.gffcompare_data[sample]['counts']['reference']['mRNA'] = c[0]
            self.gffcompare_data[sample]['counts']['reference']['loci'] = c[1]
            self.gffcompare_data[sample]['counts']['reference']['multi_exon_transcripts'] = c[2]
            
            c = [int(s) for s in lines[8].replace("(", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['super_loci'] = c[0]

            ## Accuracy metrics (Sensitivity/Precision)
            self.gffcompare_data[sample]['accuracy'] = {}
            self.gffcompare_data[sample]['missed'] = {}
            self.gffcompare_data[sample]['novel'] = {}
            for line in lines[10:16]:
                split = line.replace("|", "").replace("Intron chain", "Intron_chain").split()
                self.gffcompare_data[sample]['accuracy'][split[0]] = {}
                self.gffcompare_data[sample]['accuracy'][split[0]]['sensitivity'] = float(split[2])
                self.gffcompare_data[sample]['accuracy'][split[0]]['precision'] = float(split[3])

            ## Additional count data
            self.gffcompare_data[sample]['counts']['matching_intron_chains'] = [int(s) for s in lines[17].replace("(", " ").split() if s.isdigit()][0]
            self.gffcompare_data[sample]['counts']['matching_transcripts'] = [int(s) for s in lines[18].replace("(", " ").split() if s.isdigit()][0]
            self.gffcompare_data[sample]['counts']['matching_loci'] = [int(s) for s in lines[19].replace("(", " ").split() if s.isdigit()][0]

            ## Additional count data
            self.gffcompare_data[sample]['missed']['Exons'] = [int(s) for s in lines[21].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['novel']['Exons'] = [int(s) for s in lines[22].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['missed']['Introns'] = [int(s) for s in lines[23].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['novel']['Introns'] = [int(s) for s in lines[24].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['missed']['Loci'] = [int(s) for s in lines[25].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['novel']['Loci'] = [int(s) for s in lines[26].replace("/", " ").split() if s.isdigit()]

        print(self.gffcompare_data)

        # Filter to strip out ignored sample names
        self.gffcompare_data = self.ignore_samples(self.gffcompare_data)

        # Raise user warning if no data found
        if len(self.gffcompare_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.gffcompare_data)))

        # Add nothing to general statistics table  general statistics table
       
        # Write data file
        self.write_data_file(self.gffcompare_data, 'multiqc_gffcompare')

        # Report sectionsq
        self.add_section (
            name = "Gffcompare comparison accuracy",
            description = (
                """
                Accuracy values for comparison of 
                """
            ),
            helptext = (
                """
                There are three possible cases:
                * **Primers found**: Full length cDNA reads with correct primers at both ends.
                * **Rescued reads**: Split fusion reads.
                * **Unusable**: Reads without correct primer combinations.
                """
            ),
            anchor = 'gffcompare_accuracy',
            plot = self.plot_accuracy()
        )

        self.add_section (
            name = "Gffcompare novel loci",
            description = (
                """
                Comparison of samples with new loci
                """
            ),
            helptext = (
                """
                There are three possible cases:
                * **Primers found**: Full length cDNA reads with correct primers at both ends.
                * **Rescued reads**: Split fusion reads.
                * **Unusable**: Reads without correct primer combinations.
                """
            ),
            anchor = 'gffcompare_novel',
            plot = self.plot_novel()
        )

        self.add_section (
            name = "Gffcompare missed loci",
            description = (
                """
                Accuracy values for comparison of 
                """
            ),
            helptext = (
                """
                There are three possible cases:
                * **Primers found**: Full length cDNA reads with correct primers at both ends.
                * **Rescued reads**: Split fusion reads.
                * **Unusable**: Reads without correct primer combinations.
                """
            ),
            anchor = 'gffcompare_missed',
            plot = self.plot_missed()
        )


    # Plotting functions
    def plot_accuracy(self):
        """ Generate GffCompare accuracy plot"""
        
        datasets = ['Base', 'Exon', 'Intron', 'Intron_chain', 'Transcript', 'Locus']
        
        pconfig = {
            'id': 'gffcompare_accuracy_plot',
            'title': 'GffCompare: Accuracy values',
            'ylab': 'Precision',
            'xlab': 'Sensitivity',
            'ymin': 0,
            'ymax': 1,
            'xmin': 0,
            'xmax': 1,
            'data_labels' : [{'name' : x} for x in datasets]
        }

        data_classification = [{
                sample: {
                    'x' : self.gffcompare_data[sample]['accuracy'][dataset]['sensitivity']/100,
                    'y' : self.gffcompare_data[sample]['accuracy'][dataset]['precision']/100,
                    'name' : dataset 
                }
            for sample in self.gffcompare_data.keys()
            }
        for dataset in datasets
        ]
                        
        print(data_classification)
        return scatter.plot(data_classification, pconfig)

    # Plot number of novel and missing transcripts
    def plot_novel(self):
        """ Generate GffCompare novel elements plot"""

        cats = ['novel', 'found']
        datasets = ['Exons', 'Introns', 'Loci']
        pconfig = {
            'ymin': 0,    
            'data_labels': datasets
        }

        data_novel = [{
                sample:{
                    'novel' : self.gffcompare_data[sample]['novel'][dataset][0],
                    'found' : (self.gffcompare_data[sample]['novel'][dataset][1] -
                        self.gffcompare_data[sample]['novel'][dataset][0])
                }
            for sample in self.gffcompare_data.keys()
            }
        for dataset in datasets
        ]
        
        print(data_novel)
        return bargraph.plot(data_novel, cats, pconfig)
        
    
    def plot_missed(self):
        """ Generate GffCompare missed elements plot"""

        cats = ['missed', 'found']
        datasets = ['Exons', 'Introns', 'Loci']
        pconfig = {
            'ymin': 0,    
            'data_labels': datasets
        }

        data_missed = [{
                sample:{
                    'missed' : self.gffcompare_data[sample]['missed'][dataset][0],
                    'found' : (self.gffcompare_data[sample]['missed'][dataset][1] -
                        self.gffcompare_data[sample]['missed'][dataset][0])
                }
            for sample in self.gffcompare_data.keys()
            }
        for dataset in datasets
        ]
        
        print(data_missed)
        return bargraph.plot(data_missed, cats, pconfig)
