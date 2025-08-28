import logging
import json 

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph #für Balkendiagramm

log = logging.getLogger(__name__) 
#Logger

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name = "Deacon",
            anchor = "deacon", 
            href = "https://github.com/bede/deacon",
            info = "Search and depletion of FASTA/FASTQ files and streams using accelerated minimizer matching."

    )
        
        self.deacon_data = dict() 
        #self. verhindert lokale Variable -> für betrachtetes Objekt gültig, nicht für def __init__
        for f in self.find_log_files("deacon", filehandles = True): #filhandles = True sagt, dass Multiqc Dateien öffnen soll
            sample_name = self.clean_s_name(f["s_name"], f) #f["s_name"] -> Samplename; self.clean_s_name -> normiert den Namen(Dateipfade entfernen etc.), verhindert, dass zwei Dateien gleichen Namen haben
            try:
                data = json.load(f["f"]) #json-Inhalt der Datei wird gelesen und soll geparsed werden; f["f"] -> offenes file
            except Exception as e: #Fehler werden abgefangen in e
                log.error("failed to parse data from {f['f']} : {e}") #{e} -> liefert Fehlermeldung
                continue

            #Metriken extrahieren 
            self.deacon_data[sample_name] = {
                "k" : data.get("k"), #k-mer Länge?
                "w" : data.get("w"), #Window Größe?
                "abs_threshold" : data.get("abs_threshold"),
                "rel_threshold" : data.get("rel_threshold"),
                #"prefix_length" : data.get("prefix_length"),
                #"deplete" : data.get("deplete"),
                #"rename" : data.get("rename"),
                "seqs_in" : data.get("seqs_in"),
                "seqs_out" : data.get("seqs_out"),
                "seqs_out_proportion" : data.get("seqs_out_proportion"),
                "seqs_removed" : data.get("seqs_removed"),
                "seqs_removed_proportion" : data.get("seqs_removed_proportion"),
                "bp_in" : data.get("bp_in"),
                "bp_out" : data.get("bp_out"),
                "bp_out_proportion" : data.get("bp_out_proportion"),
                "time" : data.get("time"),
                "seqs_per_second" : data.get("seqs_per_second"),
                "bp_per_second" : data.get("bp_per_second")
            }
        
        if len(self.deacon_data) == 0:
            raise UserWarning("no deacon reports found")
        
        #Tabelle als Bericht mit Metriken
        self.add_data_table(
            self.deacon_data,
            headers = {
                "k" : {"title" : "k-mer length"},
                "w" : {"title" : "window size"},
                "abs_threshold" : {"title" : "absolute threshold"},
                "rel_threshold" : {"title" : "relative threshold"},
                "seqs_in" : {"title" : "reads in"},
                "seqs_out" : {"title" : "reads out"},
                "seqs_removed" : {"title" : "reads removed"},
                "seqs_removed_proportion" : {"title" : "reads removed (%)", "formate" : "{:%}"}, #Angabe in Prozent
                "bp_in" : {"title" : "bp in"},
                "bp_out" : {"title" : "bp out"},
                "bp_removed" : {"title" : "bp removed"},
                "bp_removed_proportion" : {"title" : "bp removed (%)", "formate" : "{:%}"},
                "time": {"title" : "time (in s)"},
                "seqs_per_second" : {"title" : "reads/s"},
                "bp_per_second" : {"title" : "bp/s"}
            },
            title = "deacon statistics",
            description = "statistics parsed from deacon reports"
        )

        #Plots hinzufügen -> Reads entfernt
        plot_data = {}
        for sample, stats in self.deacon_data.items(): #für alle samples und deren Statistiken
            if stats.get("seqs_removed_proprtion") is not None: #überprüfen, ob der Report "seqs_removed_proportion" vorhanden ist
                plot_data[sample] = stats["seqs_removed_proportion"] * 100 #Umwandlung in %

        if len(plot_data) > 0: #wenn samples mit Daten enthalten sind
            pconfig = { #Plot-Konfiguration, Dictionary; Report-Konfiguration nicht nach default (ohne pconfig wäre default)
                "id" : "deacon_removed_reads",
                "title" : "% Reads removed (Deacon)",
                "ylab" : "% removed"
            }

            self.add_section( #neue report-Sektion für MultiQC
                name = "Reads removed",
                anchor = "deacon_removed_reads",
                description = "percentage of removed reads per sample",
                plot = bargraph.plot(plot_data, pconfig) #erzeugt das Balkendiagramm
            )