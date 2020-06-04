#!/usr/bin/env python

""" MultiQC module to parse output from Kaiju """

from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Kaiju',
            anchor='kaiju',
            href="http://kaiju.binf.ku.dk/",
            info="a fast and sensitive taxonomic classification for metagenomics."
        )

        # Set up data structures
        self.kaiju_data = { }
        self.taxonomic_ranks = ["phylum","class","order","family","genus","species"]

        # Find and parse kaiju2table report
        for f in self.find_log_files('kaiju', filehandles=True):
            taxo_rank, parsed_data = self.parse_kaiju2table_report(f)
            if taxo_rank is not None and parsed_data is not None:
                # One report can have several samples, raise error if sample already viewed with this taxonomic rank
                s_names = parsed_data.keys()
                for s_name in s_names:
                    if taxo_rank in self.kaiju_data and s_name in self.kaiju_data[taxo_rank].keys() :
                        log.debug("Duplicate sample found in logs at {} rank! Overwriting sample: {}".format(taxo_rank, s_name))
                self.add_data_source(f)
                if taxo_rank in self.kaiju_data.keys() :
                    self.kaiju_data[taxo_rank].update(parsed_data)
                else :
                    self.kaiju_data[taxo_rank] = parsed_data

        # Filters to strip out ignored sample names
        for taxo_rank in self.kaiju_data :
            self.kaiju_data[taxo_rank] = self.ignore_samples(self.kaiju_data[taxo_rank])
            self.write_data_file(self.kaiju_data[taxo_rank], 'multiqc_kaiju_'+taxo_rank)

        # Number of samples found
        try:
            num_samples = max([len(taxo_rank) for taxo_rank in self.kaiju_data.values()])
        except ValueError:
            # No log files so didn't get any taxo_ranks
            raise UserWarning

        # no file found
        if num_samples == 0:
            raise UserWarning

        log.info("Found {} reports".format(num_samples))

        self.kaiju_stats_table()
        self.kaiju_taxo_rank_plot()

    def parse_kaiju2table_report(self, report):
        """ Search a kaiju with a set of regexes """
        parsed_data = {}
        s_names = []
        for l in report['f']:
            if l.startswith("file\t") :
                continue
            (s_file, pct, reads, taxon_id, taxon_names) = l.rstrip().split("\t")
            s_name = self.clean_s_name(s_file, report['root'])
            if s_name not in s_names :
                s_names.append(s_name)
                parsed_data[s_name]={"assigned":0, "unclassified":0, "cannot_be_assigned":0}
            if taxon_names.startswith("unclassified") :
                parsed_data[s_name]["unclassified"] = int(reads)
            elif taxon_names.startswith("cannot be assigned") :
                parsed_data[s_name]["cannot_be_assigned"] = int(reads)
                taxo_rank = taxon_names.split()[-1]
            else:
                parsed_data[s_name]["assigned"] += int(reads)

        for s_name in parsed_data :
            total_r = parsed_data[s_name]["assigned"] + parsed_data[s_name]["cannot_be_assigned"] + parsed_data[s_name]["unclassified"]
            parsed_data[s_name]['percentage_assigned'] = (float(parsed_data[s_name]["assigned"]) / float(total_r)) * 100.0

        return taxo_rank, parsed_data

    def kaiju_stats_table(self):
        """ Take the parsed stats from the Kaiju reports and add them to the
        basic stats table at the top of the report """
        headers = {}
        general_data ={}
        taxo_ranks = self.kaiju_data.keys()

        # print only phylum rank in general table.
        if len(taxo_ranks) >= 1 and "phylum" in taxo_ranks:
            general_data = self.kaiju_data["phylum"]
            general_taxo_rank = "Phylum"
        else:
            general_data = self.kaiju_data[self.kaiju_data.keys()[0]]
            general_taxo_rank = self.kaiju_data.keys()[0].capitalize()

        headers['percentage_assigned'] = {
            'title': '% Reads assigned {}'.format(general_taxo_rank),
            'description': 'Percentage of reads assigned at {} rank'.format(general_taxo_rank),
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        headers['assigned'] = {
            'title': '{} Reads assigned {} '.format(config.read_count_prefix, general_taxo_rank),
            'description': 'Number of reads assigned ({})  at {} rank'.format(config.read_count_desc, general_taxo_rank),
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'Blues'
        }
        self.general_stats_addcols(general_data, headers)

    def kaiju_taxo_rank_plot (self):
        """ Make the taxonomic rank barplot """

        pconfig = {
            'id': 'kaiju_taxo_rank_plot',
            'title': 'Kaiju: Reads assigned',
            'ylab': 'Number of reads',
            'cpswitch_counts_label': 'Number of reads',
            'data_labels':[]
        }
        datasets = []
        cats = []
        for rank in self.taxonomic_ranks :
            if rank in self.kaiju_data :
                pconfig['data_labels'].append( {'name': rank.capitalize(), 'ylab': 'Number of Reads at {} rank'.format(rank.capitalize()) })
                cats.append(OrderedDict ())
                cats[-1]['assigned'] = {'name': 'Assigned', 'color': '#1f78b4'}
                cats[-1]['cannot_be_assigned'] = {'name': 'Cannot be assigned at this rank', 'color': '#a6cee3'}
                cats[-1]['unclassified'] = {'name': 'Unclassified', 'color': '#ff7f00'}
                datasets.append(self.kaiju_data[rank])

        self.add_section (
            name = 'Read classification at different taxonomic levels',
            anchor = 'kaiju-classification',
            description = 'This plot shows the proportion of reads assigned, unclassified or not classified at different taxonomic rank.',
            plot = bargraph.plot(datasets, cats, pconfig)
        )
