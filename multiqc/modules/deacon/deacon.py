import logging
import json 

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph 

log = logging.getLogger(__name__) 
#logger

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name = "Deacon",
            anchor = "deacon", 
            href = "https://github.com/bede/deacon",
            info = "Search and depletion of FASTA/FASTQ files and streams using accelerated minimizer matching."

    )
        
        self.deacon_data = dict() 
        #self. makes deacon.data an instance variable, not a local for __init__

        #iterate through all detected Deacon log files
        for f in self.find_log_files("deacon", filehandles = True): #filehandles = opens the files
            sample_name = self.clean_s_name(f["s_name"], f) #f["s_name"] -> samplename; self.clean_s_name -> normalizes the name to avoid duplicates
            try:
                data = json.load(f["f"]) #load JSON data from filehandles
            except Exception as e: #catch parsing-errors in e
                log.error(f"failed to parse data from {f['f']} : {e}") #{e} -> error message
                continue

            #extract metrics for each sample
            self.deacon_data[sample_name] = {
                "k" : data.get("k"), #k-mere length
                "w" : data.get("w"), #wondow size
                "abs_threshold" : data.get("abs_threshold"), #number of k-mers that count as hit
                "rel_threshold" : data.get("rel_threshold"),
                #"prefix_length" : data.get("prefix_length"), #prefix length used for hashing/filtering
                #"deplete" : data.get("deplete"), #boolean, depleted reads (True/False)
                #"rename" : data.get("rename"), #boolean, renamed reads after export (True/False)
                "seqs_in" : data.get("seqs_in"), #number of input reads
                "seqs_out" : data.get("seqs_out"), #number of output reads
                "seqs_out_proportion" : data.get("seqs_out_proportion"), 
                "seqs_removed" : data.get("seqs_removed"), #number of removed reads
                "seqs_removed_proportion" : data.get("seqs_removed_proportion"),
                "bp_in" : data.get("bp_in"), #number of input base-pairs
                "bp_out" : data.get("bp_out"), #number of output base-pairs
                "bp_out_proportion" : data.get("bp_out_proportion"),
                "bp_removed" : data.get("bp_removed"), #number of removed base-pairs
                "bp_removed_proportion" : data.get("bp_removed_proportion"),
                "time" : data.get("time"), #runtime
                "seqs_per_second" : data.get("seqs_per_second"),
                "bp_per_second" : data.get("bp_per_second") 
            }
        
        if len(self.deacon_data) == 0:
            raise UserWarning("no deacon reports found")
        
        #create table in the report
        self.add_data_table(
            self.deacon_data,
            headers = {
                "k" : {"title" : "k-mer length"},
                "w" : {"title" : "window size"},
                "abs_threshold" : {"title" : "absolute threshold"},
                "rel_threshold" : {"title" : "relative threshold"},
                #"prefix_length" : {"title" : "prefix length"},
                #"deplete" : {"title" : "title"},
                #"rename" : {"title" : "rename"},
                #-> optionale headers?
                "seqs_in" : {"title" : "reads in"},
                "seqs_out" : {"title" : "reads out"},
                "seqs_removed" : {"title" : "reads removed"},
                "seqs_removed_proportion" : {"title" : "reads removed (%)", "format" : "{:%}"}, #Angabe in Prozent
                "bp_in" : {"title" : "bp in"},
                "bp_out" : {"title" : "bp out"},
                "bp_removed" : {"title" : "bp removed"},
                "bp_removed_proportion" : {"title" : "bp removed (%)", "format" : "{:%}"},
                "time": {"title" : "time (in s)"},
                "seqs_per_second" : {"title" : "reads/s"},
                "bp_per_second" : {"title" : "bp/s"}
            },
            title = "deacon statistics",
            description = "statistics parsed from deacon reports"
        )

        #add plot in report for: reads removed
        plot_data = {}
        for sample, stats in self.deacon_data.items(): #iterate through all samples and their statistics
            if stats.get("seqs_removed_proportion") is not None: #check, if report contains "seqs_removed_proportion"
                plot_data[sample] = stats["seqs_removed_proportion"] * 100 #convert to %

        if len(plot_data) > 0: #create if there is data
            pconfig = { #plot-configuration, dictionary; no default configuration
                "id" : "deacon_removed_reads",
                "title" : "% Reads removed (Deacon)",
                "ylab" : "% removed"
            }

            self.add_section( #new section in MultiQC report
                name = "Reads removed",
                anchor = "deacon_removed_reads",
                description = "percentage of removed reads per sample",
                plot = bargraph.plot(plot_data, pconfig) #generates a barplot
            )

