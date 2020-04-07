#!/usr/bin/env python

""" MultiQC module to parse output from mirtop"""

from __future__ import print_function
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule
from collections import OrderedDict
import logging
import json

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='mirtop', anchor='mirtop',
        href='https://github.com/miRTop/mirtop/',
        info="is a command line tool to annotate miRNAs and isomiRs and compute general statistics using the mirGFF3 format."
        )
        
        # Find and load any mirtop reports
        self.mirtop_data = dict()
        for f in self.find_log_files('mirtop'):
            self.parse_mirtop_report(f)
        # Filter out ignored samples (given with --ignore-samples option)
        self.mirtop_data = self.ignore_samples(self.mirtop_data)

        # Raise error if dict is empty
        if len(self.mirtop_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.mirtop_data)))

        # Write parsed report data to a file
        self.write_data_file(self.mirtop_data, 'multiqc_mirtop')

        # Create summary table
        self.mirtop_stats_table()

        # Create detailed plots
        self.mirtop_barplot_section(aggregate_snps=True)
        
    def parse_mirtop_report (self, f):
        """ Parse the mirtop log file. """
        
        log.info("Processing file " + f['fn']  )
        content = json.loads(f['f'])
        for sample_name in content['metrics'].keys():
            log.info("Importing sample " + sample_name)
            ## Check for sample name duplicates
            if sample_name in self.mirtop_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))
            parsed_data = dict()
            parsed_data = content['metrics'][sample_name]
            parsed_data['read_count'] = parsed_data['isomiR_sum'] + parsed_data['ref_miRNA_sum']
            parsed_data['isomiR_perc'] = (parsed_data['isomiR_sum'] / parsed_data['read_count']) * 100
            self.mirtop_data[sample_name] = parsed_data

    def mirtop_stats_table(self):
        """ Take the parsed stats from the mirtop report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['ref_miRNA_sum'] = {
            'title': 'Reference reads',
            'description': 'Read counts summed over all reference miRNAs',
            'shared_key': 'reads',
            'scale': 'PuBu'
        }
        headers['isomiR_sum'] = {
            'title': 'IsomiR reads',
            'description': 'Read counts summed over all isomiRs',
            'shared_key': 'reads',
            'scale': 'Oranges'
        }
        headers['read_count'] = {
            'title': 'Total reads',
            'description': 'Total read counts (both isomiRs and reference miRNA)',
            'shared_key': 'reads',
            'scale': 'BuGn'
        } 
        headers['isomiR_perc'] = {
            'title': 'IsomiR %',
            'description': '% of total read counts corresponding to isomiRs',
            'min':0,
            'max':100,
            'suffix':'%',
            'scale': 'YlOrRd'
        }
        self.general_stats_addcols(self.mirtop_data, headers)

    def mirtop_barplot_section(self, aggregate_snps):
        """ Generate barplots for a given stat type"""

        # Define info for each plot
        plots_info = OrderedDict()
        general_helptext = " The different isomiR types are: **iso_3p**, a sequence with a 3' end difference because of trimming or templated tailing; **iso_5p**, a sequence with a 5' end difference because of trimming or templated tailing; **iso_add3p**, a sequence with non templated tailing in the 3' end; **iso_add5p**, a sequence with non templated tailing in the 5' end; **iso_snv**, a sequence with a single nucleotide variant. The **ref_miRNA** label corresponds to the reference miRNA (canonical sequence)."
        plots_info["sum"] = {
            "name" : "Read counts for each isomiR type",
            "description" : "Total counts of reads aligned for each isomiR type, over all detected miRNAs. The total counts of reads detected as reference miRNA sequences is also shown. Since a read can belong to 2 (or more) different isomiRs types (e.g iso_3p and iso_5p), the cumulated read counts shown in this plot for a sample can be higher than its total read count shown in the general statistics.",
            "helptext" : "For each sample, the mean counts of each type of isomiRs over all detected miRNAs is displayed in a different color." + general_helptext}
        plots_info["count"] = {
            "name" : "Number of unique sequences for each isomiR type",
            "description" : "For each isomiR type, number of distinct sequences detected, over all miRNAs. The number of reference miRNA sequences detected is also shown.",
            "helptext" : "For each sample, the number of miRNAs with each type of isomiRs, is displayed in a different color." + general_helptext}
        plots_info["mean"] = {
            "name" : "Means for each isomiR type",
            "description" : "Mean counts, for each isomiR type, over all detected miRNAs. The mean counts of reads detected as reference miRNA sequences is also shown.",
            "helptext" : "For each sample, the mean counts of each type of isomiRs over all detected miRNAs is displayed in a different color." + general_helptext}

        # Each isomiR cat (combining mirtop current and previous versions)
        cats = ["ref_miRNA", "iso_3p", "iso_5p", "iso_add3p", "iso_add5p", "iso_add", "iso_snv", "iso_snv_seed", "iso_snv_central_offset", "iso_snv_central", "iso_snv_central_supp", "iso_snp", "iso_snp_seed", "iso_snp_central_offset", "iso_snp_central", "iso_snp_central_supp"]

        # Y axis label definition
        ylab_dict = {"sum": "Read counts (stacked)", "count": "Unique sequences (stacked)", "mean": "Means (stacked)"}
        
        # Aggregate infos for iso_snp isomiRs (for clarity). "Mean" section will be recomputed
        def aggregate_snps_in_samples():
            snv_aggr = {} ## sub dict with all infos except for snps
            for sample in self.mirtop_data:
                snv_aggr[sample]={key:self.mirtop_data[sample][key] for key in self.mirtop_data[sample] if "iso_snp" not in key}
                snv_aggr[sample]['iso_snv_sum'] = sum([self.mirtop_data[sample][key] for key in self.mirtop_data[sample] if "iso_snp" in key and "sum" in key])
                snv_aggr[sample]['iso_snv_count'] = sum([self.mirtop_data[sample][key] for key in self.mirtop_data[sample] if "iso_snp" in key and "count" in key])
                if snv_aggr[sample]['iso_snv_count'] > 0:
                    snv_aggr[sample]['iso_snv_mean'] = snv_aggr[sample]['iso_snv_sum'] / snv_aggr[sample]['iso_snv_count']
                else:
                    snv_aggr[sample]['iso_snv_mean'] = 0
            return snv_aggr
                
        if aggregate_snps:
            data = aggregate_snps_in_samples()
        else:
            data = self.mirtop_data 
            
        # Gather data and add plot section
        for plot_key in plots_info:
            plot_data = dict()
            log.info("Plotting " + plot_key + " section." )
            for sample in data:
                # Select keys with plot_key for each sample (except keys with "isomiRs_*" and "read_*")
                plot_data[sample]={key:data[sample][key] for key in data[sample] if plot_key in key and "isomiR" not in key and "read" not in key}
            cats_section = {key + "_" + plot_key:{'name':key} for key in cats}
            config = {"id": plot_key, "title": "mirtop: " + plots_info[plot_key]["name"], "ylab": ylab_dict[plot_key]}
            self.add_section (
                name = plots_info[plot_key]["name"],
                description = plots_info[plot_key]["description"],
                helptext = plots_info[plot_key]["helptext"],
                plot = bargraph.plot(plot_data, cats_section, config)
            )

