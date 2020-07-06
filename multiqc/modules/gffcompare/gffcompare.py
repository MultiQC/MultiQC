#!/usr/bin/env python

""" MultiQC module to parse output from gffcompare """

from __future__ import print_function
import logging
import os
from collections import OrderedDict

from multiqc.plots import linegraph, bargraph
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
        # Everything is hardcoded, needs to be adjusted if output should change. 
        # I hope there are no other versions of gffcompare that produce different output formats
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
            self.gffcompare_data[sample]['counts']['missed_exons'] = [int(s) for s in lines[21].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['novel_exons'] = [int(s) for s in lines[22].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['missed_introns'] = [int(s) for s in lines[23].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['novel_introns'] = [int(s) for s in lines[24].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['missed_loci'] = [int(s) for s in lines[25].replace("/", " ").split() if s.isdigit()]
            self.gffcompare_data[sample]['counts']['novel_loci'] = [int(s) for s in lines[26].replace("/", " ").split() if s.isdigit()]

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

        # Report sections
        self.add_section (
            name = "Gffcompare comparison accuracy",
            description = (
                """
                This plot shows the cDNA read categories identified by Pychopper </br>
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
            name = "GffCompare  Strand Orientation",
            description = (
                """
                This plot shows the strand orientation of full length cDNA reads
                """
            ),
            helptext = (
                """
                Nanopore cDNA reads are always read forward. To estimate their original strand, 
                Pychopper searches for the location of the start and end primers and assigns the reads accordingly.
                """
            ),
            anchor = 'pychopper_orientation',
            plot = self.plot_orientation()
        )

    # Plotting functions
    def plot_classification(self):
        """ Generate the cDNA read classification plot """

        pconfig = {
            'id': 'pychopper_classification_plot',
            'title': 'Pychopper: Read classification',
            'ylab': '',
            'xDecimals': False,
            'ymin': 0
        }

        data_classification = {}
        for sample in self.pychopper_data.keys():
            data_classification[sample] = {}
            data_classification[sample]=self.pychopper_data[sample]['Classification']

        cats = ['Primers_found', 'Rescue', 'Unusable']
        return bargraph.plot(data_classification, cats, pconfig)

    def plot_orientation(self):
        """ Generate the read strand orientation plot """

        pconfig = {
            'id': 'pychopper_orientation_plot',
            'title': 'Pychopper: Strand Orientation',
            'ylab': '',
            'cpswitch_c_active': False,  
            'xDecimals': False,
            'ymin': 0
        }

        data_orientation = {}
        for sample in self.pychopper_data.keys():
            data_orientation[sample] = {}
            data_orientation[sample]=self.pychopper_data[sample]['Strand']

        cats = ['+', '-']
        return bargraph.plot(data_orientation, cats, pconfig)