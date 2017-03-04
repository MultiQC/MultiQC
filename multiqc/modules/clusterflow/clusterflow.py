#!/usr/bin/env python

""" MultiQC module to parse output from Cluster Flow """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Cluster Flow module class, parses run logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Cluster Flow', anchor='clusterflow',
        href='http://clusterflow.io',
        info="is a simple and flexible bioinformatics pipeline tool.")

        # Find and load any Cluster Flow reports
        self.clusterflow_commands = dict()

        for f in self.find_log_files(config.sp['clusterflow'], filehandles=True):
            self.parse_clusterflow_logs(f)
            self.add_data_source(f)

        if len(self.clusterflow_commands) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} pipelines".format(len(self.clusterflow_commands)))

        self.sections = list()

        # Common Commands
        self.sections.append({
            'name': 'Commands',
            'anchor': 'clusterflow-commands',
            'content': self.clusterflow_commands_table()
        })


    def parse_clusterflow_logs(self, f):
        """ Parse Clusterflow logs """
        module = None
        job_id = None
        pipeline_id = None
        for l in f['f']:

            # Get pipeline ID
            module_r = re.match(r'Module:\s+(.+)$', l)
            if module_r:
                module = module_r.group(1)
            job_id_r = re.match(r'Job ID:\s+(.+)$', l)
            if job_id_r:
                job_id = job_id_r.group(1)
                if module is not None:
                    pipeline_r = re.match(r"(cf_.+)_"+re.escape(module)+r"_\d+$", job_id)
                    if pipeline_r:
                        pipeline_id = pipeline_r.group(1)
            
            # Get commands that have been run
            if l.startswith('###CFCMD'):
                if pipeline_id is None:
                    pipeline_id = 'unknown'
                if pipeline_id not in self.clusterflow_commands.keys():
                    self.clusterflow_commands[pipeline_id] = list()
                self.clusterflow_commands[pipeline_id].append(l[8:])


    def clusterflow_commands_table (self):
        """ Generate the trimming length plot """
        
        # I wrote this when I was tired. Sorry if it's incomprehensible.

        html = '''<p>Every Cluster Flow run will have many different commands.
            MultiQC splits these by whitespace, collects by the tool name
            and shows the first command found. Any terms not found in <em>all</em> subsequent
            calls are replaced with <code>[variable]</code>
            <em>(typically input and ouput filenames)</em>.</p>'''

        # Loop through pipelines
        tool_cmds = OrderedDict()
        for pipeline_id, commands in self.clusterflow_commands.items():

            self.var_html = '<span style="background-color:#dedede; color:#999;">[variable]</span>'
            tool_cmd_parts = OrderedDict()
            for cmd in commands:
                s = cmd.split()
                if s[0] not in tool_cmd_parts.keys():
                    tool_cmd_parts[s[0]] = list()
                tool_cmd_parts[s[0]].append(s)

            
            for tool, cmds in tool_cmd_parts.items():
                cons_cmd = self._replace_variable_chunks(cmds)
                # Try again with first two blocks if all variable
                variable_count = cons_cmd.count(self.var_html)
                if variable_count == len(cmds[0]) - 1 and len(cmds[0]) > 2:
                    for subcmd in set([s[1] for s in cmds]):
                        sub_cons_cmd = self._replace_variable_chunks([cmd for cmd in cmds if cmd[1] == subcmd])
                        tool = "{} {}".format(tool, subcmd)
                        if tool not in tool_cmds:
                            tool_cmds[tool] = dict()
                        tool_cmds[tool][pipeline_id] = '<code style="white-space:nowrap;">{}</code>'.format(" ".join(sub_cons_cmd) )
                else:
                    if tool not in tool_cmds:
                        tool_cmds[tool] = dict()
                    tool_cmds[tool][pipeline_id] = '<code style="white-space:nowrap;">{}</code>'.format(" ".join(cons_cmd) )
        
        table_config = {
            'id': 'clusterflow-commands',
            'table_title': 'Cluster Flow Commands',
            'col1_header': 'Base Tool',
            'sortRows': False,
            'no_beeswarm': True
        }
        return html + table.plot(tool_cmds, None, table_config)


    def _replace_variable_chunks(self, cmds):
        """ List through a list of command chunks. Return a single list
        with any variable bits blanked out. """
        
        cons_cmd = None
        while cons_cmd is None:
            for cmd in cmds:
                if cons_cmd is None:
                    cons_cmd = cmd[:]
                else:
                    for idx, s in enumerate(cons_cmd):
                        if s not in cmd:
                            cons_cmd[idx] = self.var_html
        return cons_cmd
