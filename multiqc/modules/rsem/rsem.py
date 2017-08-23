#!/usr/bin/env python

""" MultiQC module to parse output from RSEM/rsem-calculate-expression """

from __future__ import print_function
import logging
from collections import OrderedDict

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    RSEM module class, parses .cnt file .
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Rsem', anchor='rsem',
        href='https://deweylab.github.io/RSEM/',
        info="RSEM (RNA-Seq by Expectation-Maximization) is a software package for"\
             "estimating gene and isoform expression levels from RNA-Seq data.")
    

        self.rsem_mapped_data = dict()
        # Find and load any count file
        for f in self.find_log_files('rsem'):
            self.rsem_mapped_data[f['s_name']] = self.parse_rsem_report(f)
            self.add_data_source(f)
        # Filter to strip out ignored sample names
        self.rsem_mapped_data = self.ignore_samples(self.rsem_mapped_data)

        if len(self.rsem_mapped_data) == 0  :
            log.info("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.rsem_mapped_data)))

        # Write parsed report data to a file
        self.write_data_file(self.rsem_mapped_data, 'multiqc_rsem')

        # Basic Stats Table
        self.rsem_stats_table()

        # Assignment bar plot
        self.add_section(
            name = 'RSEM mapped reads',
            anchor = 'Transcript_associated',
            helptext = 'Repartition of reads: uniquely mapped, multimapped, filtered, unalignable and uncertain.',
            plot = self.rsem_chart() )


    def parse_rsem_report (self, f):
        """ Parse the rsem cnt stat file.
        Description of cnt file found : https://github.com/deweylab/RSEM/blob/master/cnt_file_description.txt
        """
        data = dict()
        for l in f['f'].splitlines():
            s = l.split()
            if len(s) > 3:
                #Line: N0 N1 N2 N_tot
                # N0, number of unalignable reads;
                # N1, number of alignable reads; (nUnique + nMulti = N1;)
                # N2, number of filtered reads due to too many alignments; N_tot = N0 + N1 + N2
                data['Unalignable'] = int(s[0])
                data['Filtered'] = int(s[2])
                data['Alignable'] = int(s[1])
                data['alignable_percent'] = (data['Alignable'] / (data['Unalignable']+data['Filtered']+data['Alignable'])) *100
            elif len(s) == 3:
                #Line: nUnique    nMulti  nUncertain
                # nUnique, number of reads aligned uniquely to a gene; nMulti, number of reads aligned to multiple genes; nUnique + nMulti = N1;
                # nUncertain, number of reads aligned to multiple locations in the given reference sequences, which include isoform-level multi-mapping reads
                # Numbers are not concordant with description, for the moment I'd rather not using them.
                data['Unique'] = s[0]
                data['Multi'] = s[1]
                data['Uncertain'] = s[2]
            else:
                break
        return data
    
    def rsem_stats_table(self):
        """ Take the parsed stats from the rsem report and add them to the
        basic stats table at the top of the report """
        headers = OrderedDict()

        headers['alignable_percent'] = {
            'title': '% alignable'.format(config.read_count_prefix),
            'description': '% alignable mapped reads'.format(config.read_count_desc),
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        self.general_stats_addcols(self.rsem_mapped_data, headers)


    def rsem_chart (self):
        """ Make the rsem assignment rates plot """
        
        # Config for the plot
        config = {
            'id': 'rsem_assignment_plot',
            'title': 'RSEM assignments',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        return bargraph.plot(self.rsem_mapped_data, ['Unalignable','Alignable','Filtered'],config)
