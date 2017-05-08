#!/usr/bin/env python

""" MultiQC module to parse output from Cluster Flow """

from __future__ import print_function
from collections import OrderedDict
import datetime
import logging
import re
import os
import time

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
        self.clusterflow_runfiles = dict()

        for f in self.find_log_files('clusterflow/logs', filehandles=True):
            self.parse_clusterflow_logs(f)
            self.add_data_source(f, 'log')
        for f in self.find_log_files('clusterflow/runfiles', filehandles=True):
            parsed_data = self.parse_clusterflow_runfiles(f)
            if parsed_data is not None:
                self.clusterflow_runfiles[f['s_name']] = parsed_data
                self.add_data_source(f, 'runfile')

        # Filters to strip out ignored sample names
        self.clusterflow_commands = self.ignore_samples(self.clusterflow_commands)
        self.clusterflow_runfiles = self.ignore_samples(self.clusterflow_runfiles)

        if len(self.clusterflow_commands) == 0 and len(self.clusterflow_runfiles) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Count pipelines
        num_log_pipelines = len(self.clusterflow_commands)
        num_runfile_pipelines = len(set(d.get('pipeline_id', 'unknown') for d in self.clusterflow_runfiles.values()))
        log.info("Found {} pipelines".format(max(num_log_pipelines, num_runfile_pipelines)))
        log.debug("Found {} log pipelines".format(num_log_pipelines))
        log.debug("Found {} runfile pipelines".format(num_runfile_pipelines))

        # Pipeline Info
        if len(self.clusterflow_runfiles) > 0:
            self.clusterflow_pipelines_section()

        # Commands
        if len(self.clusterflow_commands) > 0:
            self.clusterflow_commands_table()


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
        """ Make a table of the Cluster Flow commands """

        # I wrote this when I was tired. Sorry if it's incomprehensible.

        desc = '''Every Cluster Flow run will have many different commands.
            MultiQC splits these by whitespace, collects by the tool name
            and shows the first command found. Any terms not found in <em>all</em> subsequent
            calls are replaced with <code>[variable]</code>
            <em>(typically input and ouput filenames)</em>. Each column is for one Cluster Flow run.'''

        # Loop through pipelines
        tool_cmds = OrderedDict()
        headers = dict()
        for pipeline_id, commands in self.clusterflow_commands.items():
            headers[pipeline_id] = {'scale': False}
            self.var_html = '<span style="background-color:#dedede; color:#999;">[variable]</span>'
            tool_cmd_parts = OrderedDict()
            for cmd in commands:
                s = cmd.split()
                tool = self._guess_cmd_name(s)
                if tool not in tool_cmd_parts.keys():
                    tool_cmd_parts[tool] = list()
                tool_cmd_parts[tool].append(s)


            for tool, cmds in tool_cmd_parts.items():
                cons_cmd = self._replace_variable_chunks(cmds)
                # Try again with first two blocks if all variable
                variable_count = cons_cmd.count(self.var_html)
                if variable_count == len(cmds[0]) - 1 and len(cmds[0]) > 2:
                    for subcmd in set([x[1] for x in cmds]):
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
            'namespace': 'Cluster Flow',
            'id': 'clusterflow-commands',
            'table_title': 'Cluster Flow Commands',
            'col1_header': 'Tool',
            'sortRows': False,
            'no_beeswarm': True
        }
        self.add_section (
            name = 'Commands',
            anchor = 'clusterflow-commands',
            description = desc,
            plot = table.plot(tool_cmds, headers, table_config)
        )


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

    def _guess_cmd_name(self, cmd):
        """ Manually guess some known command names, where we can
        do a better job than the automatic parsing. """

        # zcat to bowtie
        if cmd[0] == 'zcat' and 'bowtie' in cmd:
            return 'bowtie'
        # samtools
        if cmd[0] == 'samtools':
            return ' '.join(cmd[0:2])
        # java (eg. picard)
        if cmd[0] == 'java':
            jars = [s for s in cmd if '.jar' in s]
            return os.path.basename(jars[0].replace('.jar', ''))
        return cmd[0]


    def parse_clusterflow_runfiles(self, f):
        """ Parse run files generated by Cluster Flow.
        Currently gets pipeline IDs and associated steps."""

        data = dict()
        in_comment = False
        seen_pipeline = False
        cf_file = False
        for l in f['f']:
            l = l.rstrip()
            # Check that this is from Cluster Flow
            if 'Cluster Flow' in l:
                cf_file = True
            # Header
            if l.startswith('Pipeline: '):
                data['pipeline_name'] = l[10:]
            if l.startswith('Pipeline ID: '):
                data['pipeline_id'] = l[13:]
            if l.startswith('Created at '):
                data['pipeline_start'] = l[11:]
            # Config settings
            if l.startswith('@'):
                s = l.split(None, 1)
                key = s[0].replace('@', '').strip()
                try:
                    data[key] = "\t".join(s[1:])
                except IndexError:
                    data[key] = True
            # Comments
            if l.startswith('/*'):
                in_comment = True
            if l.startswith('*/'):
                in_comment = False
            if in_comment:
                if 'comment' not in data:
                    data['comment'] = ''
                data['comment'] += l+"\n"
            # Pipeline steps
            if l.strip().startswith('#'):
                if 'pipeline_steps' not in data:
                    data['pipeline_steps'] = []
                data['pipeline_steps'].append(l)
                seen_pipeline = True
            # Step output files
            elif seen_pipeline:
                s = l.split("\t")
                if len(s) > 1:
                    if 'files' not in data:
                        data['files'] = OrderedDict()
                    if s[0] not in data['files']:
                        data['files'][s[0]] = []
                    data['files'][s[0]].append(s[1:])
        # Parse the start date
        dt = None
        if 'pipeline_id' in data:
            s = data['pipeline_id'].split('_')
            dt = datetime.datetime.fromtimestamp(int(s[-1]))
        elif 'pipeline_start' in data:
            dt_r = re.match(r'(\d{2}):(\d{2}), (\d{2})-(\d{2})-(\d{4})', data['pipeline_start'])
            if dt_r:
                dt = datetime.datetime(
                    int(dt_r.group(5)), # year
                    int(dt_r.group(4)), # month
                    int(dt_r.group(3)), # day
                    int(dt_r.group(1)), # hour
                    int(dt_r.group(2))  # minute
                )

        # Not a Cluster Flow file (eg. Nextflow .run file)
        if not cf_file:
            return None

        if dt is not None:
            data['pipeline_start_dateparts'] = {
                'year':        dt.year,
                'month':       dt.month,
                'day':         dt.day,
                'hour':        dt.hour,
                'minute':      dt.minute,
                'second':      dt.second,
                'microsecond': dt.microsecond,
                'timestamp':   time.mktime(dt.timetuple())
            }
        # Cluster Flow v0.4 and before did not print the pipeline ID in run files
        # Try to guess - will be wrong as no microsecond info, but hopefully unique
        # and reproducible for other run files
        if 'pipeline_id' not in data:
            if 'pipeline_name' in data and 'pipeline_start_dateparts' in data:
                log.debug('Trying to guess pipeline ID for file "{}"'.format(f['fn']))
                data['pipeline_id'] = 'cf_{}_{}'.format(data['pipeline_name'], data['pipeline_start_dateparts']['timestamp'])
        return data

    def clusterflow_pipelines_section(self):
        """ Generate HTML for section about pipelines, generated from
        information parsed from run files. """
        data = dict()
        pids_guessed = ''
        for f,d in self.clusterflow_runfiles.items():
            pid = d.get('pipeline_id', 'unknown')
            if d.get('pipeline_id_guess', False) is True:
                pid += '*'
                pids_guessed = ' Project IDs with an asterisk may be inaccurate.'
            # Count the number of files going into the first module
            num_starting_files = 0
            for step_name, files in d.get('files',{}).items():
                if step_name.startswith('start'):
                    num_starting_files += len(files)
            # Reformat the date so that column sorting works nicely
            if 'pipeline_start_dateparts' in d:
                dt = d['pipeline_start_dateparts']
                d['pipeline_start'] = '{}-{:02d}-{:02d} {:02d}:{:02d}'.format(dt['year'], dt['month'], dt['day'], dt['hour'], dt['minute'])
            if pid not in data:
                data[pid] = d
                data[pid]['num_starting_files'] = int(num_starting_files)
            else:
                data[pid]['num_starting_files'] += int(num_starting_files)

        headers = OrderedDict()
        headers['pipeline_name'] = {'title': 'Pipeline Name'}
        headers['pipeline_start'] = {'title': 'Date Started', 'description': 'Date and time that pipeline was started (YYYY-MM-DD HH:SS)'}
        headers['genome'] = {'title': 'Genome ID', 'description': 'ID of reference genome used'}
        headers['num_starting_files'] = {'title': '# Starting Files', 'format': '{:,.0f}', 'description': 'Number of input files at start of pipeline run.'}
        table_config = {
            'namespace': 'Cluster Flow',
            'id': 'clusterflow-pipelines',
            'table_title': 'Cluster Flow Pipelines',
            'col1_header': 'Pipeline ID',
            'no_beeswarm': True,
            'save_file': True
        }
        self.add_section (
            name = 'Pipelines',
            anchor = 'clusterflow-pipelines',
            description = 'Information about pipelines is parsed from <code>*.run</code> files. {}'.format(pids_guessed),
            plot = table.plot(data, headers, table_config),
            content = self.clusterflow_pipelines_printout()
        )


    def clusterflow_pipelines_printout(self):
        """ Print the steps used in each Cluster Flow pipeline """
        data = dict()
        html = ''
        for f,d in self.clusterflow_runfiles.items():
            pid = d.get('pipeline_id', 'unknown')
            data[pid] = [
                d.get('pipeline_name'),
                "\n".join(d.get('pipeline_steps', []))
            ]
        for pid, d in data.items():
            html += '''
                <div class="panel panel-default">
                    <div class="panel-heading"><h3 class="panel-title">Pipeline Steps: {} (<code>{}</code>)</h3></div>
                    <pre class="panel-body" style="border:0; background-color:transparent; padding:0 15px; margin:0; color:#666; font-size:90%;">{}</pre>
                </div>
                '''.format(pid, d[0], d[1])
        return html

