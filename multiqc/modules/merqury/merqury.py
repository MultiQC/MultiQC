#!/usr/bin/env python

""" MultiQC module to parse output from MERQURY """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
from multiqc.plots import linegraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Merqury module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="MERQURY",
            anchor="merqury",
            href="https://github.com/marbl/merqury",
            info="Evaluate genome assemblies with k-mers and more.",
            doi="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02134-9",
        )


        # Find and load any MERQURY reports
        self.completeness_data = dict()
        for f in self.find_log_files("merqury/completeness"):
            self.parse_completeness_log(f)
        
        self.qv_data = dict()
        for f in self.find_log_files("merqury/qv"):
            self.parse_qv_log(f)
        
        self.spectra_data = dict()
        for f in self.find_log_files("merqury/spectra"):
            self.parse_spectra_log(f)

        log.info("Found {} completeness reports".format(len(self.completeness_data)))
        log.info("Found {} qv reports".format(len(self.qv_data)))
        log.info("Found {} spectra histograms".format(len(self.spectra_data)))

        # Write parsed report data to a file
        self.write_data_file(self.completeness_data, "multiqc_merqury_completeness")
        
        # Write parsed report data to a file
        self.write_data_file(self.qv_data, "multiqc_merqury_qv")
        
        # Write parsed report data to a file
        self.write_data_file(self.spectra_data, "multiqc_merqury_spectra")
        
        #self.merqury_general_stats_table()
        self.add_section(name="k-mer completeness (recovery rate)", anchor="merqury-completeness", plot=self.completeness_table())
        
        self.add_section(name="QV estimation", anchor="merqury-qv", plot=self.qv_table())
        
        for s_name in self.spectra_data.keys():
            self.add_section(name="Spectra plots for "+s_name, anchor="merqury-spectra", plot=self.spectra_plot(self.spectra_data[s_name],s_name)) 
    
    def parse_completeness_log(self, f):
        self.add_data_source(f,f["s_name"])
        s_name=f["s_name"].replace("_output_merqury.completeness.tab","")
        for l in f["f"].splitlines():
            s = l.strip().split("\t")
            suffix = s[0]
            s_val1= s[2]
            s_val2=s[3]
            s_val3=s[4]
            self.completeness_data[s_name+"_"+suffix] = dict()
            self.completeness_data[s_name+"_"+suffix]["val1"]=s_val1
            self.completeness_data[s_name+"_"+suffix]["val2"]=s_val2
            self.completeness_data[s_name+"_"+suffix]["percent"]=s_val3
  
    def parse_qv_log(self, f):
        self.add_data_source(f,f["s_name"])
        s_name=f["s_name"].replace("_output_merqury.qv.tab","")
        for l in f["f"].splitlines():
            s = l.strip().split("\t")
            suffix = s[0]
            s_val1= s[1]
            s_val2=s[2]
            s_val3=s[3]
            s_val4=s[4]
            self.qv_data[s_name+"_"+suffix] = dict()
            self.qv_data[s_name+"_"+suffix]["val1"]=s_val1
            self.qv_data[s_name+"_"+suffix]["val2"]=s_val2
            self.qv_data[s_name+"_"+suffix]["qv"]=s_val3
            self.qv_data[s_name+"_"+suffix]["error"]=s_val4
   
    def parse_spectra_log(self,f):
        self.add_data_source(f,f["s_name"])
        s_name=f["s_name"].split(".")[0]
        p_name=".".join(f["s_name"].split(".")[1:])
        nameddata = dict()
        nameddata[s_name]=dict()
        data = dict()
        for l in f["f"].splitlines()[1:]:
            s = l.strip().split("\t")
            serie=s[0]
            if serie not in data.keys(): data[serie]={}
            data[serie][int(s[1])] = int(s[2])
        if s_name not in self.spectra_data.keys(): self.spectra_data[s_name]={}
        self.spectra_data[s_name][p_name]=data

    def spectra_plot(self,d,f):
        config={'data_labels':[]}
        out=[]
        xmax=0
        ymax=0
        for plot_key in d.keys():
            data=d[plot_key]
            colors=["#000000","#cc0000","#246bce","#93c47d","#8e7cc3","#e69138"]
            colorscale= {}
            for i in range(len(data.keys())):
                serie=list(data.keys())[i]
                colorscale[serie]=colors[i]
                if i>0 and i<3:
                    ymax=max([ymax,max(data[serie].values())])
                    peak=max(data[serie].values())
                    xmax=max([xmax,[k for k, v in data[serie].items() if v == peak][0]])
            config["data_labels"].append({'name':plot_key,"ylab": "kmer count", "xlab": "kmer multiplicity","colors":colorscale,"legend":True,"showInLegend":False,"use_legend":False})
            config["xmax"]=xmax*2.1
            config["ymax"]=ymax*1.1
            config['height']=350
            out.append(data)
        return linegraph.plot(out, config)
        
    def completeness_table(self):
        """Take the parsed stats from the QUAST report and add some to the
        General Statistics table at the top of the report"""

        headers = OrderedDict()
        headers["val1"] = {
            "title": "k-mers in assembly",
            "description": "k-mers in assembly",
            "min": 0,
            #"suffix": self.contig_length_suffix,
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        headers["val2"] = {
            "title": "k-mers in read set",
            "description": "k-mers in read set",
            "min": 0,
            #"suffix": self.contig_length_suffix,
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        headers["percent"] = {
            "title": "completeness (%)",
            "description": "kmer completeness",
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        config = {
            "id": "completeness_table",
            "namespace": "MERQURY",
            "min": 0,
        }
        #self.general_stats_addcols(self.merqury_data, headers)
        return table.plot(self.completeness_data, headers, config)
        

    def qv_table(self):
        """Take the parsed stats from the QUAST report and add some to the
        General Statistics table at the top of the report"""

        headers = OrderedDict()
        headers["val1"] = {
            "title": "uniq k-mers",
            "description": "k-mers uniquely found only in the assembly",
            "min": 0,
            #"suffix": self.contig_length_suffix,
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        headers["val2"] = {
            "title": "common k-mers",
            "description": "k-mers found in both assembly and the read set",
            "min": 0,
            #"suffix": self.contig_length_suffix,
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        headers["qv"] = {
            "title": "QV",
            "description": "estimated quality value",
            "min": 0,
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        headers["error"] = {
            "title": "error rate",
            "description": "estimated error rate",
            "min": 0,
            "scale": "RdYlGn",
            "format": "{:,.8f}",
        }
        config = {
            "id": "qv_table",
            "namespace": "MERQURY",
            "min": 0,
        }
        #self.general_stats_addcols(self.merqury_data, headers)
        return table.plot(self.qv_data, headers, config)
        

