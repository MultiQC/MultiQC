#!/usr/bin/env python

""" MultiQC module to parse output from checkatlas """

from __future__ import print_function

import os.path
from collections import OrderedDict
import logging
import re

from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """checkatlas module, parses stderr logs."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="CheckAtlas",
            anchor="checkatlas",
            href="https://github.com/becavin-lab/checkatlas",
            info="An easy to use program to check the quality of your single-cell atlases.",
            doi="",
        )

        self.data_summary = dict()
        for f in self.find_log_files("checkatlas/summary"):
            input_fname = f['s_name'].replace('_checkatlas','')
            s_name = self.clean_s_name(input_fname, f)
            self.data_summary[s_name] = self.parse_table_logs(f['f'])
            self.add_data_source(f, s_name)

        self.data_umap = dict()
        for f in self.find_log_files("checkatlas/umap"):
            input_fname = f['s_name'].replace('_checkatlas_umap','')
            s_name = self.clean_s_name(input_fname, f)
            print(s_name)
            print(f)
            self.data_umap[s_name] = f['f']
            self.add_data_source(f, s_name)

        self.data_adata = dict()
        for f in self.find_log_files("checkatlas/adata"):
            input_fname = f['s_name'].replace('_checkatlas_adata','')
            s_name = self.clean_s_name(input_fname, f)
            self.data_adata[s_name] = self.parse_table_logs(f['f'])
            self.add_data_source(f, s_name)


        if len(self.data_summary) == 0:
            raise UserWarning
        if len(self.data_adata) == 0:
            raise UserWarning
        if len(self.data_umap) == 0:
            raise UserWarning

        if len(self.data_summary) > 0:
            log.info("Found {} reports".format(len(self.data_summary)))

        self.write_data_file(self.data_summary, 'multiqc_checkatlas-summary')
        self.write_data_file(self.data_adata, 'multiqc_checkatlas_adata')
        # self.general_stats_addcols(self.checkatlas_data)

        self.add_sections()


    def parse_table_logs(self, f):
        """
        Parse logs which are .tsv files
        Ex: _checkatlas.tsv, adata_checkatlas.tsv
        :param f:
        :return:
        """
        data = {}
        lines = f.splitlines()
        headers = lines[0].split(',')
        for i in range(1, len(lines)):
            line = lines[i].split(',')
            for j in range(0, len(line)):
                data[headers[j]] = line[j]
        return data

    def add_sections(self):
        self.add_summary_section()
        self.add_adata_section()
        self.add_umap_section()
        self.add_clustermetrics_section()
        self.add_annotationmetrics_section()
        self.add_specificitymetrics_section()
        self.add_dimredmetrics_section()


    def add_summary_section(self):
        self.add_section(
            name='Atlas overview',
            anchor='checkatlas-summary',
            description='Overview of your single-cell atlases',
            helptext="""
                
                """,
            content=table.plot(self.data_summary)
        )

    def add_adata_section(self):
        self.add_section(
            name='AnnData explorer',
            anchor='checkatlas-anndata',
            description='Exploration of your AnnData',
            helptext="""
            
            """,
            content=table.plot(self.data_adata)
        )

    def add_umap_section(self):
        html_content = """
        <div class="tab">\n"""
        for key, value in self.data_umap.items():
            html_content += """  <button class="tablinks" onclick="openUmap(event, '"""+key+"""')">"""+key+ \
            """</button>\n"""
        """
        </div>
        """
        for key, value in self.data_umap.items():
            path_fig = '../' + value.name.replace('./','')
            print(path_fig)
            html_content += """<div id=\""""+key+"""\" class="tabcontent">
            <img src=\""""+path_fig+"""\" >\n</div>\n"""

        html_content += """
        <script>
        function openUmap(evt, pngName) {
          var i, tabcontent, tablinks;
          tabcontent = document.getElementsByClassName("tabcontent");
          for (i = 0; i < tabcontent.length; i++) {
            tabcontent[i].style.display = "none";
          }
          tablinks = document.getElementsByClassName("tablinks");
          for (i = 0; i < tablinks.length; i++) {
            tablinks[i].className = tablinks[i].className.replace(" active", "");
          }
          document.getElementById(pngName).style.display = "block";
          evt.currentTarget.className += " active";
        }
        </script>
        """

        self.add_section(
            name='UMAP visualisation',
            anchor='checkatlas-umap',
            description='UMAP of your atlases.',
            helptext="""
            
            """,
            content=html_content
        )

    def add_clustermetrics_section(self):
        self.add_section(
            name='Classification metrics',
            anchor='checkatlas-clustmetrics',
            description='Quality control metrics calculated on your atlases.',
            helptext="""
            
            """,
            content="",
        )

    def add_annotationmetrics_section(self):
        self.add_section(
            name='Annotaton metrics',
            anchor='checkatlas-annotmetrics',
            description='Quality control metrics calculated on your atlases.',
            helptext="""
            
            """,
            content='',
        )

    def add_specificitymetrics_section(self):
        self.add_section(
            name='Specificity metrics',
            anchor='checkatlas-specmetrics',
            description='Quality control metrics calculated on your atlases.',
            helptext="""
            
            """,
            content='',
        )

    def add_dimredmetrics_section(self):
        self.add_section(
            name='Dimensionality reduction metrics',
            anchor='checkatlas-dimredmetrics',
            description='Quality control metrics calculated on your atlases.',
            helptext="""
            
            """,
            content='',
        )