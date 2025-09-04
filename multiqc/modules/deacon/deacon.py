import logging
import json 
import re

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph 
from multiqc.plots import table

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
            if not f["fn"].endswith("deacon.json"): #every json-report ends with deacon.json
                continue
            try:
                data = json.load(f["f"]) #load JSON data from filehandles
            except Exception as e: #catch parsing-errors in e
                log.error(f"failed to parse data from {f['f']} : {e}") #{e} -> error message
                continue
            
            #check, if data is a dict
            if not isinstance(data, dict):
                log.warning("fSkipping {f['fn']} : parsed JSON is not an oject")
                continue 
            
            #check, if version is version
            version = data.get("version")
            if version is None or not version.startswith("deacon"):
                log.warning(f"Skipping {f['fn']} : not a deacon report (version = {version})")
                continue
            
            sample_name = self.clean_s_name(f["s_name"], f) #f["s_name"] -> samplename; self.clean_s_name -> normalizes the name to avoid duplicates


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
            log.warning("no deacon reports found")
        
        #create table in the report
        self.add_section(
            name = "Deacon statistics",
            anchor = "deacon_stats",
            description = "Statistics parsed from deacon reports",
            plot = table.plot(self.deacon_data, {
                "k" : {"title" : "k-mer length"},
                "w" : {"title" : "window size"},
                "abs_threshold" : {"title" : "absolute threshold"},
                "rel_threshold" : {"title" : "relative threshold"},
                #"prefix_length" : {"title" : "prefix length"},
                #"deplete" : {"title" : "depleted reads"},
                #"rename" : {"title" : "renamed reads"},
                "seqs_in" : {"title" : "reads in"},
                "seqs_out" : {"title" : "reads out"},
                "seqs_removed" : {"title" : "reads removed"},
                "seqs_removed_proportion" : {"title" : "reads removed (%)", "format" : "{:%}"},
                "bp_in" : {"title" : "bp in"},
                "bp_out" : {"title" : "bp out"},
                "bp_removed" : {"title" : "bp removed"},
                "bp_removed_proportion" : {"title" : "bp removed (%)", "format" : "{:%}"},
                "time": {"title" : "time (in s)"},
                "seqs_per_second" : {"title" : "reads/s"},
                "bp_per_second" : {"title" : "bp/s"}
            })
        )

        #add plot in report for: reads removed
        plot_data = {}
        for sample, stats in self.deacon_data.items(): #iterate through all samples and their statistics
            if stats.get("seqs_removed_proportion") is not None: #check, if report contains "seqs_removed_proportion"
                plot_data[sample] = stats["seqs_removed_proportion"]

        if len(plot_data) > 0: #create if there is data
            pconfig_plot = { #plot-configuration, dictionary; no default configuration
                "id" : "deacon_removed_reads",
                "title" : "% Reads removed (Deacon)",
                "ylab" : "% removed"
            }

            self.add_section( #new section in MultiQC report
                name = "Reads removed",
                anchor = "deacon_removed_reads",
                description = "percentage of removed reads per sample",
                plot = bargraph.plot(plot_data, pconfig_plot) #generates a barplot
            )


        #add plot in report for: Reads removed and remaining in absolute and percentage number
        plot2_data = {}
        for sample, stats in self.deacon_data.items(): #iterate through all samples and their statistics
            if stats.get("seqs_removed") is not None and stats.get("seqs_out") is not None: #check, if report contains "seqs_removed" and "seqs_out"
                plot2_data[sample] = {
                    "Reads removed" : stats["seqs_removed"],
                    "Reads remaining" : stats["seqs_out"]
                }
            
        if len(plot2_data) > 0: #create if there is data
            pconfig_plot2 = {
                "id" : "Reads_removed_remaining",
                "title" : "Reads removed and Reads remaining",
                "ylab" : "Number of Reads",
                "cpswitch" : True,  #switch between absolute and percentage number
                "stacked" : True    #stacked bars
            }

            self.add_section( #new section in MultiQC report
                name = "Reads removed and Reads remaining",
                anchor = "Reads_removed_remaining",
                description = "Comparison between removed and remaining Reads, switch between absolute and percentage number.",
                plot = bargraph.plot(plot2_data, pconfig_plot2) #generate a barplot
            )