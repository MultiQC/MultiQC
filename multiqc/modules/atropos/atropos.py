"""MultiQC module that reads JSON output from Atropos.

Note: this code is mostly borrowed from the Cutadapt [1] and FastQC [2] MultiQC
modules.
1. https://github.com/ewels/MultiQC/blob/master/multiqc/modules/fastqc/fastqc.py
2. https://github.com/ewels/MultiQC/blob/master/multiqc/modules/cutadapt/cutadapt.py
"""
from __future__ import print_function, division, absolute_import
from collections import OrderedDict, Iterable, defaultdict
import logging
import io
import math
from numbers import Number
import operator
import os
import json
from numpy import mean, multiply, cumsum
from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

## Hard-coded URLs
# TODO: eventually these need pointers to specific sections of the docs
# TODO: handle multiplexed output
# TODO: compute read1 and read2 status separately
# TODO: look at whether there are any plots to duplicate from here: https://github.com/pnnl/fqc

ATROPOS_GITHUB_URL = "https://github.com/jdidion/atropos"
ATROPOS_DOC_URL = "http://atropos.readthedocs.org/en/latest/guide.html"

class Status(object):
    def __init__(self, name, val, color):
        self.name = name
        self.val = val
        self.color = color
    
    def __lt__(self, other):
        return self.val < other.val
    
    def __eq__(self, other):
        return self.val == other.val
    
    def __repr__(self):
        return self.name

FAIL = Status('fail', 0, '#d9534f')
WARN = Status('warn', 1, '#f0ad4e')
PASS = Status('pass', 2, '#5cb85c')

DEFAULT_COLOR = '#999'

ADD_LISTENERS = '<script type="text/javascript">$(function () {{ add_atropos_listeners("{phase}", atropos_passfails_{phase})  }});</script>'
PASSFAILS = '<script type="text/javascript">atropos_passfails_{phase} = {statuses};</script>'

# Initialise the logger
log = logging.getLogger(__name__)

## Module

class MultiqcModule(BaseMultiqcModule):
    """Base Atropos module class. Loads JSON ouptut. Submodule classes are
    responsible for extracting data and generating summary table and plots.
    
    Args:
        name: Module name.
        anchor: HTML anchor name.
        from_config: Whether to load data from configured locations. This should
            only be set to False in tests.
    """
    def __init__(self, name='', anchor='', info='', from_config=True):
        # Create the submodule based on the phase
        mod_cust_config = getattr(self, 'mod_cust_config', {})
        phase = mod_cust_config.get('phase', 'trim')
        submodule_class = SUBMODULES[phase]
        self.submodule = submodule_class()

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name=name or self.submodule.name,
            anchor=anchor or self.submodule.anchor,
            href=ATROPOS_GITHUB_URL,
            info=info or self.submodule.info + (
                "is a general-purpose NGS pre-processing tool that "
                 "specializes in adatper- and quality-trimming, written "
                 "by John Didion at NHGRI. It is a fork of Cutadapt, "
                 "written by Marcel Martin at Stockholm University."))

        self.css = {}
        self.js = {}

        # Load data from logs and generate the report
        if from_config:
            self.init_from_config()
            self.submodule.atropos_report(self)

    def __getattr__(self, name):
        """Fetch all missing attributes from the submodule.
        """
        return getattr(self.submodule, name)

    def init_from_config(self):
        """Initialize module from JSON files discovered via MultiQC
        configuration.
        """
        log_files = self.find_log_files('atropos', filehandles=True)
        if not log_files:
            # Atropos summary files are JSON files (with .json extension) and
            # always have '"program": "Atropos"' as the first key-value pair.
            # Note that this way is deprecated; currently it produces a warning
            # but it may produce an error in the future.
            patterns = dict(fn='*.json', contents='"program": "Atropos"')
            log_files = self.find_log_files(patterns, filehandles=True)

        for file_dict in log_files:
            fileobj = file_dict['f']
            data = json.load(fileobj)
            # We always run the trim submodule because it adds to the general
            # stats table.
            if self.submodule.can_add_data(data):
                self.add_atropos_data(data, fileobj, file_dict)

        if self.num_samples == 0:
            log.debug("Could not find any reports in %s", config.analysis_dir)
            raise UserWarning
        else:
            log.info("Found %s reports", self.num_samples)

    def add_atropos_data(self, data, data_source, file_dict):
        """For each sample, Atropos generates a <sample>.json file. That
        file has summary information (numbers of files processed, description of
        the inputs (e.g. FASTA vs FASTQ, single vs paired-end), and a summary
        of the analyses performed.

        Args:
            data: Atropos summary dict.
            data_source: The file from which 'data' was loaded.
        """
        close = False
        if data_source and isinstance(data_source, str):
            data_source = open(data_source)
            close = True
        if data is None:
            if data_source:
                data = json.load(data_source)
                if close:
                    data_source.close()
            else:
                raise ValueError("One of 'data' or 'data_source' must be provided")

        sample_id = data['sample_id']
        if file_dict:
            self.add_data_source(file_dict, sample_id)

        self.submodule.add_data(sample_id, data)

# Submodules

class Submodule(object):
    def __init__(self):
        # List of sample IDs
        self.atropos_sample_ids = []
        # General data (for summary table)
        self.atropos_general_data = OrderedDict()

    @property
    def num_samples(self):
        return len(self.atropos_sample_ids)

    def can_add_data(self, data):
        return True

    def add_data(self, sample_id, data):
        self.atropos_sample_ids.append(sample_id)

    def atropos_report(self, parent):
        """Generate the report.
        """
        if not self.atropos_sample_ids:
            log.debug("No reports to process; atropos_report raising UserWarning")
            raise UserWarning

        # Set assets on parent
        self.atropos_assets(parent)

        # Add to general stats
        if self.atropos_general_data:
            self.atropos_general(parent)

        self.atropos_plots(parent)

    def atropos_assets(self, parent):
        """Add assets to parent."""
        pass

    def atropos_general(self, parent):
        """Add items to the parent summary table.
        """
        raise NotImplementedError()

    def atropos_plots(self, parent):
        """Add plots to the parent report.
        """
        raise NotImplementedError()


class TrimModule(Submodule):
    name = 'Atropos: Trimming'
    anchor = 'atropos'
    phase = 'trim'
    info = ''

    def __init__(self):
        super(TrimModule, self).__init__()
        self.atropos_trim_data = OrderedDict()
        self.trim_sections = [
            TrimmedLength()
        ]

    def add_data(self, sample_id, data):
        super(TrimModule, self).add_data(sample_id, data)
        self.atropos_general_data[sample_id] = data['derived'].copy()
        self.atropos_general_data[sample_id]['mean_sequence_length'] = \
            mean(data['derived']['mean_sequence_lengths'])
        self.atropos_general_data[sample_id].update(dict(
            (key, data[key])
            for key in (
                'total_record_count',
                'total_bp_counts')))
        
        if 'trim' in data:
            self.atropos_trim_data[sample_id] = data['trim']
            self.atropos_general_data[sample_id]['fraction_bp_trimmed'] = \
                1.0 - data['trim']['formatters']['fraction_total_bp_written']
            for mod_dict in data['trim']['modifiers'].values():
                if 'fraction_records_with_adapters' in mod_dict:
                    self.atropos_general_data[sample_id]['fraction_records_with_adapters'] = \
                        mean(mod_dict['fraction_records_with_adapters'])
                    break

    def atropos_assets(self, parent):
        parent.css['atropos-assets/css/multiqc_atropos.css'] = os.path.join(
            os.path.dirname(__file__), 'assets', 'css', 'multiqc_atropos.css')
        parent.js['atropos-assets/js/multiqc_atropos.js'] = os.path.join(
            os.path.dirname(__file__), 'assets', 'js', 'multiqc_atropos.js')

    def atropos_general(self, parent):
        """Add some single-number stats to the basic statistics table at the
        top of the report.
        """
        headers = OrderedDict()
        headers['input_format'] = {
            'title': 'Format',
            'description': 'Input File Format'
        }
        headers['total_record_count'] = {
            'title': 'M Seqs',
            'description': 'Total Sequences (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['sum_total_bp_count'] = {
            'title': 'M bp',
            'description': 'Total Base Pairs (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'base_count'
        }
        headers['mean_sequence_length'] = {
            'title': 'Length',
            'description': 'Average Sequence Length (bp)',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'format': '{:.0f}'
        }
        headers['fraction_bp_trimmed'] = {
            'title': '% Base Pairs Trimmed',
            'description': '% Total Base Pairs trimmed',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev',
            'modify': lambda x: x * 100,
            'format': '{:.1f}%'
        }
        headers['fraction_records_with_adapters'] = {
            'title': '% Reads w/ Adapters',
            'description': '% Total Reads with Adapters',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev',
            'modify': lambda x: x * 100,
            'format': '{:.1f}%'
        }
        print('adding general stats', str(self.atropos_general_data), str(headers))
        parent.general_stats_addcols(self.atropos_general_data, headers)
    
    def atropos_plots(self, parent):
        for section in self.trim_sections:
            context = section(self.atropos_trim_data, self.atropos_general_data)
            if 'plot' in context:
                parent.add_section(**context['plot'])


class QcModule(Submodule):
    def __init__(self):
        super(QcModule, self).__init__()
        self.atropos_qc_data = OrderedDict()
        self.qc_sections = [
            PerBaseQuality(),
            #PerTileQuality()
            PerSequenceQuality(),
            PerBaseContent(),
            PerSequenceGC(),
            PerBaseN(),
            #SequenceLength(),
            #Duplication(),
            #Contaminants(),
            #Kmers()
        ]

    def can_add_data(self, data):
        return self.phase in data

    def add_stats_data(self, sample_id, pairing, input_names, source, data):
        self.atropos_qc_data[sample_id] = OrderedDict()
        source_files = input_names
        if pairing == 3:
            self.atropos_qc_data[sample_id] = [None, None]
            for idx in range(2):
                read = 'read{}'.format(idx+1)
                if read in data:
                    self.atropos_qc_data[sample_id][idx] = data[read]
                    self.atropos_qc_data[sample_id][idx]['source'] = source_files[idx]
        else:
            read = 'read{}'.format(pairing)
            data[read]['source'] = source_files[pairing-1]
            self.atropos_qc_data[sample_id] = [data[read], None]

    def atropos_assets(self, parent):
        # Colours to be used for plotting lines
        parent.status_colours = dict(
            (status.name, status.color)
            for status in (PASS, WARN, FAIL))
        parent.status_colours['default'] = DEFAULT_COLOR

    def atropos_general(self, parent):
        headers = OrderedDict()
        headers['percent_fails_{}'.format(self.phase)] = {
            'title': '% Failed ({}-trimming)'.format(self.phase),
            'description': 'Percentage of Atropos QC modules failed ({}-trimming)'.format(self.phase),
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
            'format': '{:.0f}%',
            'hidden': True
        }
        parent.general_stats_addcols(self.atropos_general_data, headers)
    
    def atropos_plots(self, parent):
        # Add statuses to intro. Note that this is slightly different than
        # FastQC: Atropos reports the relevant statistic, and the threshold
        # for pass/warn/fail is configured in MutliQC (defaulting to
        # the thresholds defined in FastQC).
        statuses = {}
        fails = defaultdict(int)
        for section in self.qc_sections:
            section_name = "{}_{}".format(section.name, self.phase)
            context = section(self.atropos_qc_data, self.atropos_general_data, phase=self.phase)
            if all(key in context for key in ('statuses', 'plot')):
                statuses[section_name] = context['statuses']
                for sample_id, status in statuses[section_name].items():
                    if status == FAIL:
                        fails[sample_id] += 1
                parent.add_section(**context['plot'])
        
        for sample_id, num_fails in fails.items():
            self.atropos_general_data[sample_id]['percent_fails_{}'.format(self.phase)] = \
                num_fails * 100 / len(self.qc_sections)
        
        parent.intro += ADD_LISTENERS.format(phase=self.phase)
        parent.intro += PASSFAILS.format(phase=self.phase, statuses=json.dumps(simplify(statuses)))
    
class PreModule(QcModule):
    name = 'Atropos: Pre-trim QC'
    info = 'reports pre-trimming statistics. Atropos '
    anchor = 'atropos-pre'
    phase = 'pre'

    def add_data(self, sample_id, data):
        super(PreModule, self).add_data(sample_id, data)
        self.atropos_general_data[sample_id] = {}
        pairing = data['input']['input_read']
        input_names = data['input']['input_names']
        for source, phase_data in data['pre'].items():
            self.add_stats_data(sample_id, pairing, input_names, source, phase_data)

class PostModule(QcModule):
    name = 'Atropos: Post-trim QC'
    info = 'reports post-trimming statistics. Atropos '
    anchor = 'atropos-post'
    phase = 'post'

    def can_add_data(self, data):
        return super(PostModule, self).can_add_data(data) and 'NoFilter' in data['post']
    
    def add_data(self, sample_id, data):
        super(PostModule, self).add_data(sample_id, data)
        self.atropos_general_data[sample_id] = {}
        pairing = data['input']['input_read']
        input_names = data['input']['input_names']
        for source, post_data in data['post']['NoFilter'].items():
            self.add_stats_data(sample_id, pairing, input_names, source, post_data)

SUBMODULES = dict(
    trim=TrimModule,
    pre=PreModule,
    post=PostModule
)

# Base Section class

class Section(object):
    """Base class for QC sections.

    The __call__ method returns a context that includes both the statuses
    and the plot. There are several steps that are all methods that subclasses
    can override for customization.

    statuses: get_statuses -> get_status (for each sample) -> get_status_for_pair -> get_status_for
    plot: get_plot -> get_html -> add_plot_data, add_html_variables
    """
    name = \
    display = \
    anchor = \
    plot_type = \
    plot_config = \
    html = None
    #headers = {}
    
    def __call__(self, data, general, **kwargs):
        section_data = self.prepare_data(data, general)
        context = kwargs.copy()
        self.add_to_context(context, section_data)
        return context
    
    def prepare_data(self, data, general):
        """Returns the data for this section from the parent dict.
        """
        raise NotImplementedError()
    
    def add_to_context(self, context, data):
        """Add plot and any intermediate data to context.
        """
        if data is not None:
            context['plot'] = self.get_plot(context, data)
    
    def get_plot(self, context, data):
        """Returns the plot dict.
        """
        return {
            'name': self.display,
            'anchor': self.anchor,
            'content': self.get_html(context, data)
        }
    
    def get_html(self, context, data):
        """Returns the formatted HTML.
        """
        self.add_plot_data(context, data)
        return self.format_html(context)
    
    def add_plot_data(self, context, data):
        """Add plot data to the context. Plot data must be named 
        'plot_data{read}', where read is 1 or 2.
        """
        plot_data1 = OrderedDict()
        plot_data2 = OrderedDict()
        for sample_id, (data1, data2) in data.items():
            plot_data1[sample_id] = self.get_sample_plot_data(context, data1)
            if data2:
                plot_data2[sample_id] = self.get_sample_plot_data(context, data2)
        
        content1 = self.get_read_plot(context, plot_data1, 1)
        if plot_data2:
            content2 = self.get_read_plot(context, plot_data2, 2)
            context['plot'] = self.wrap_plots(content1, content2)
        else:
            context['plot'] = content1
    
    def get_sample_plot_data(self, context, data):
        """Get the plot data for a single sample.
        """
        raise NotImplementedError()
    
    def get_read_plot(self, context, plot_data, read):
        """Returns the plot for a single read.
        """
        return self.plot_type.plot(plot_data, self.get_plot_config(context, read=read))
    
    def get_plot_config(self, context, **kwargs):
        """Returns the plot config dict.
        """
        return self.plot_config
    
    def wrap_plots(self, *content):
        """Wrap a pair of plots in a <div>.
        """
        return '<div class="pair-plot-wrapper">' + ''.join(content) + '</div>'
    
    def format_html(self, context):
        """Format the html string with variables from the context.
        """
        return self.html.format(**context)

# Trim sections

class TrimmedLength(Section):
    """There is one plot (or pair of plots) per adapter. Data is stored as
    adapter_seq: {front/back: [<read1>, <read2>]}, where front/back is the end 
    from which the adapter is trimmed (so there are four plots in the case of 
    adapters that are trimmed at both the front and back (i.e. linked and
    anywhere)), and where each read is [<obs>,<exp>], where obs and exp are 
    each a {sample_id: {x:y}} dict for the observed and expected distributions.
    """
    name = 'trimmed_length'
    display = 'Trimmed length'
    anchor = 'atropos_trimmed_length'
    header_html = """
<p>This plot shows the number of reads with certain lengths of adapter trimmed.
Obs/Exp shows the raw counts divided by the number expected due to sequencing
errors. A defined peak may be related to adapter length. See the
<a href="{}" target="_blank">Atropos documentation</a> for more information on
how these numbers are generated.</p>
""".format(ATROPOS_DOC_URL)
    
    def prepare_data(self, data, general):
        def add_lengths(side, read, sample_id, num_reads, adapter_length, read_dict, dest_dict):
            for key in ('lengths_{}', '{0}_lengths_{0}'):
                key = key.format(side)
                if key in read_dict:
                    if dest_dict[side] is None:
                        dest_dict[side] = [[{},{}], [{},{}]]
                    hist = sort_hist(read_dict[key])
                    dest_dict[side][read][0][sample_id] = hist
                    dest_dict[side][read][1][sample_id] = ordered_dict(
                        (length, num_reads * 0.25 ** min(length, adapter_length))
                        for length, count in hist.items())
                    return
        
        adapter_data = {}
        for sample_id, trim_data in data.items():
            if 'modifiers' not in trim_data:
                continue
            num_reads = general[sample_id]['total_record_count']
            for modifier_dict in trim_data['modifiers'].values():
                if 'adapters' in modifier_dict:
                    for read, read_dict in enumerate(modifier_dict['adapters']):
                        for adapter_id, adapter_dict in read_dict.items():
                            seq = adapter_dict['sequence']
                            if seq not in adapter_data:
                                adapter_data[seq] = dict(front=None, back=None)
                            seqlen = len(seq)
                            add_lengths(
                                'front', read, sample_id, num_reads, seqlen, 
                                adapter_dict, adapter_data[seq])
                            add_lengths(
                                'back', read, sample_id, num_reads, seqlen, 
                                adapter_dict, adapter_data[seq])
        
        return adapter_data
    
    def get_html(self, context, data):
        html = self.header_html
        for adapter_seq, adapter_data in data.items():
            adapter_plots = []
            for side in sorted(adapter_data.keys()):
                reads = adapter_data[side]
                if reads is None:
                    continue
                for read, read_data in enumerate(reads, 1):
                    if not any(read_data):
                        continue
                    adapter_plots.append(linegraph.plot(
                        read_data, self.get_plot_config(
                            context, seq=adapter_seq, side=side, read=read)))
            html += self.wrap_plots(*adapter_plots)
        return html
    
    def get_plot_config(self, context, seq='', side='back', read=1):
        return {
            'id': 'atropos_trimmed_length_plot_{seq}_{side}_{read}'.format(
                seq=seq, side=side, read=read),
            'title': 'Lengths of Trimmed Sequences<br/>Seq={seq}, Side={side}, Read={read}'.format(
                seq=seq, side=side, read=read),
            'ylab': 'Counts',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0, 
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'Counts', 'ylab': 'Count'},
                            {'name': 'Obs/Exp', 'ylab': 'Observed / Expected'}]
        }

# QC sections

class QcSection(Section):
    default_thresholds = None
    compare = operator.lt
    
    def get_plot(self, context, data):
        """Returns the plot dict.
        """
        return {
            'name': "{}".format(self.display),
            'anchor': "{}_{}".format(self.anchor, context['phase']),
            'content': self.get_html(context, data)
        }

    def prepare_data(self, data, general):
        """Returns the data for this section from the parent dict.
        """
        section_data = {}
        for sample_id, (read1_data, read2_data) in data.items():
            sample_data = [None, None]
            if read1_data:
                sample_data[0] = read1_data[self.name]
            if read2_data:
                sample_data[1] = read2_data[self.name]
            section_data[sample_id] = sample_data
        return section_data
    
    def add_to_context(self, context, data):
        context['statuses'] = self.get_statuses(context, data)
        #if section.headers:
        #    self.general_stats_addcols(section_data, section.headers)
        context['plot'] = self.get_plot(context, data)
    
    def get_statuses(self, context, data):
        """Returns {sample_id: status} dict.
        """
        return ordered_dict(
            (sample_id, self.get_status(context, *sample_data))
            for sample_id, sample_data in data.items())
    
    def get_status(self, context, sample_data1, sample_data2=None):
        """Returns the status. If this is paired data, the status should
        be the min of the read1 and read2 statuses.
        """
        stat1, stat2 = [
            self.compute_statistic(context, data)
            for data in (sample_data1, sample_data2)]
        return self.get_status_for_pair(context, stat1, stat2)
    
    def compute_statistic(self, context, data):
        """Compute the threshold statistic.
        """
        raise NotImplementedError()
    
    def get_status_for_pair(self, context, stat1, stat2=None):
        """Returns the status for one or a pair of threshold statistics.
        """
        status1 = self.get_status_for(context, stat1)
        if stat2 is not None:
            status2 = self.get_status_for(context, stat2)
            return min(status1, status2)
        else:
            return status1
    
    def get_status_for(self, context, stat):
        """Returns the status for single threshold statistic.
        """
        # TODO: Add configuration for thresholds.
        for status, threshold in zip((FAIL, WARN), self.default_thresholds):
            if self.compare(stat, threshold):
                return status
        return PASS

## Static HTML/templates

class PerBaseQuality(QcSection):
    name = 'base_qualities'
    display = 'Sequence Quality Histograms'
    anchor = 'atropos_base_qualities'
    html = """
<p>The mean quality value across each base position in the read. See the
<a href="{}" target="_blank">Atropos help</a>.
</p>
{{plot}}""".format(ATROPOS_DOC_URL)
    # TODO: Add boxplots as in FastQC output.
    plot_type = linegraph
    
    def get_plot_config(self, context, read=1):
        return {
            'id': 'atropos_base_qualities_plot_{phase}_{read}'.format(
                phase=context['phase'], read=read),
            'title': 'Mean Quality Scores: Read {}'.format(read),
            'ylab': 'Phred Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            'colors': get_status_cols(context['statuses']),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
    
    def compute_statistic(self, context, data):
        # base_qualities is in the form 
        # { columns: qualities, rows: {pos: counts} }
        # Note: The FastQC documentation is incorrect - they do not compute
        # quartiles/medians on a per-base basis, but rather on a per-bin
        # basis, where bin sizes are determined by the read length.
        # TODO: average over bins
        quals = data['columns']
        min_lower_quartile = None
        min_median = None
        for pos, base_counts in data['rows'].items():
            counts_cumsum = cumsum(base_counts)
            lower_quartile = weighted_lower_quantile(
                quals,
                (float(c) / counts_cumsum[-1] for c in base_counts),
                0.25)
            if min_lower_quartile is None or lower_quartile < min_lower_quartile:
                min_lower_quartile = lower_quartile
            median = weighted_median(quals, counts_cumsum)
            if min_median is None or median < min_median:
                min_median = median
        return (min_lower_quartile, min_median)
    
    def get_status_for(self, context, stat):
        # TODO: decide how to average over bins. For now we only look at 
        # medians.
        min_lower_quartile, min_median = stat
        #if min_lower_quartile < 5 or min_median < 20:
        if min_median < 20:
            return FAIL
        #if min_lower_quartile < 10 or min_median < 25:
        if min_median < 25:
            return WARN
        return PASS
    
    def get_sample_plot_data(self, context, data):
        return hist_to_means(data)

class PerTileQuality(QcSection):
    name = 'tile_sequence_qualities'
    default_thresholds = (-5, -2)
    # TODO

class PerSequenceQuality(QcSection):
    name = 'qualities'
    display = 'Per Sequence Quality Scores'
    anchor = 'atropos_qualities'
    default_thresholds = (20, 27)
    html = """
<p>The number of reads with average quality scores. Shows if a subset of reads
has poor quality. See the <a href="{}" target="_blank">Atropos help</a>.</p>
<p>{{plot}}</p>""".format(ATROPOS_DOC_URL)
    plot_type = linegraph
    
    def compute_statistic(self, context, data):
        return data['summary']['modes'][0]
    
    def get_plot_config(self, context, read=1):
        return {
            'id': 'atropos_qualities_plot_{phase}_{read}'.format(
                phase=context['phase'], read=read),
            'title': 'Per Sequence Quality Scores: Read {}'.format(read),
            'ylab': 'Count',
            'xlab': 'Mean Sequence Quality (Phred Score)',
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': get_status_cols(context['statuses']),
            'tt_label': '<b>Phred {point.x}</b>: {point.y} reads',
            'xPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
    
    def get_sample_plot_data(self, context, data):
        return dict((int(pos), count) for pos, count in data['hist'].items())

class PerBaseContent(QcSection):
    name = 'bases'
    display = 'Per Base Sequence Content'
    anchor = 'atropos_bases'
    compare = operator.gt
    default_thresholds = (0.2, 0.1)
    header_html = """
<p>The proportion of each base position for which each of the four normal DNA
bases has been called. See the <a href="{}" target="_blank">Atropos help</a>.</p>
<p class="text-primary"><span class="glyphicon glyphicon-info-sign"></span>
Click a heatmap row to see a line plot for that dataset.</p>
""".format(ATROPOS_DOC_URL)
    plot_html = """
<div id="atropos_bases_plot_{phase}_{read}" class="atropos_bases_plot">
    <h4 class="atropos_seq_heatmap_header">Per Base Sequence Content: Read {read}</h4>
    <h5><span class="s_name"><em class="text-muted">rollover for sample name</em></span></h5>
    <div class="atropos_seq_heatmap_key">
        Position: <span id="atropos_seq_heatmap_{phase}_{read}_key_pos">-</span>
        <div><span id="atropos_seq_heatmap_{phase}_{read}_key_t"> %T: <span>-</span></span></div>
        <div><span id="atropos_seq_heatmap_{phase}_{read}_key_c"> %C: <span>-</span></span></div>
        <div><span id="atropos_seq_heatmap_{phase}_{read}_key_a"> %A: <span>-</span></span></div>
        <div><span id="atropos_seq_heatmap_{phase}_{read}_key_g"> %G: <span>-</span></span></div>
    </div>
    <div id="atropos_seq_heatmap_{phase}_{read}_div" class="atropos-overlay-plot">
        <div id="atropos_seq_{phase}_{read}" class="hc-plot">
            <canvas id="atropos_seq_heatmap_{phase}_{read}" height="100%" width="800px" style="width:100%;"></canvas>
        </div>
    </div>
    <div class="clearfix"></div>
</div>
<script type="text/javascript">
    atropos_seq_content_heatmap_{phase}_{read} = new SeqContentHeatmap("{phase}", {read}, {data}, atropos_passfails_{phase});
    $(function () {{ atropos_seq_content_heatmap_{phase}_{read}.draw(); }});
</script>"""
    
    def compute_statistic(self, context, data):
        max_diff = 0
        for row in data['rows'].values():
            acgt = row[:4]
            frac = tuple(base / sum(acgt) for base in acgt)
            max_diff = max(max_diff,
                abs(frac[0]-frac[3]),
                abs(frac[1]-frac[2]))
        return max_diff
    
    def get_sample_plot_data(self, context, data):
        return dict(
            (pos, dict(zip(('a','c','g','t'), row[:4])))
            for pos, row in data['rows'].items())
    
    def format_html(self, context):
        return self.header_html + context['plot']
    
    def get_read_plot(self, context, plot_data, read):
        """Create the epic HTML for the Atropos sequence content heatmap.
        """
        return self.plot_html.format(
            phase=context['phase'], read=read, data=json.dumps(plot_data))

class PerSequenceGC(QcSection):
    name = 'gc'
    display = 'Per Sequence GC Content'
    anchor = 'atropos_gc'
    compare = operator.gt
    default_thresholds = (0.3, 0.15)
    html =  """
<p>The average GC content of reads. Normal random library typically have a
roughly normal distribution of GC content. See the <a href="{}" target="_blank">
Atropos help</a>.</p>
<p>The dashed black line shows theoretical GC content: {{theoretical_gc_name}}.</p>
{{plot}}""".format(ATROPOS_DOC_URL)
    
    def add_to_context(self, context, data):
        theoretical_gc = theoretical_gc_name = None
        
        tgc = getattr(config, 'atropos_config', {}).get('atropos_theoretical_gc', None)
        if tgc is not None:
            theoretical_gc_name = os.path.basename(tgc)
            tgc_fn = 'atropos_theoretical_gc_{}.txt'.format(tgc)
            tgc_path = os.path.join(os.path.dirname(__file__), 'atropos_theoretical_gc', tgc_fn)
            if not os.path.isfile(tgc_path):
                tgc_path = tgc
            try:
                with io.open (tgc_path, "r", encoding='utf-8') as f:
                    theoretical_gc_raw = f.read()
                    theoretical_gc = []
                    for l in theoretical_gc_raw.splitlines():
                        if '# FastQC theoretical GC content curve:' in l:
                            theoretical_gc_name = l[39:]
                        elif not l.startswith('#'):
                            s = l.split()
                            try:
                                theoretical_gc.append(
                                    [float(s[0]), float(s[1])])
                            except (TypeError, IndexError):
                                pass
            except IOError:
                log.warning(
                    "Couldn't open FastQC Theoretical GC Content file %s",
                    tgc_path)
        
        means = []
        sds = []
        max_total = 0

        for sample_id, (data1, data2) in data.items():
            for read_data in (data1, data2):
                if theoretical_gc is None:
                    means.append(read_data['summary']['mean'])
                    sds.append(read_data['summary']['stdev'])
                max_total = max(max_total, sum(read_data['hist'].values()))
        
        if theoretical_gc is None:
            normdist = NormalDistribution(mean(means), mean(sds))
            theoretical_gc = [(i, normdist(i)) for i in range(101)]
            theoretical_gc_name = "Empirical"
        
        context['theoretical_gc'] = theoretical_gc
        context['theoretical_gc_name'] = theoretical_gc_name
        context['max_total'] = max_total
        
        super(PerSequenceGC, self).add_to_context(context, data)
    
    def compute_statistic(self, context, data):
        maxcount = max(data['hist'].values())
        dev_pct = 0
        for pct in range(101):
            dev_pct += abs(
                min(
                    context['theoretical_gc'][pct][1] * 
                    context['max_total'], maxcount) -
                data['hist'].get(pct, 0))
        return dev_pct * 100 / context['max_total']

    def get_sample_plot_data(self, context, data):
        def normalize_gc(gc, total):
            return ordered_dict(
                (pct, count * 100 / total)
                for pct, count in gc)
        gchist = dict((int(pos), count) for pos, count in data['hist'].items())
        total = sum(gchist.values())
        gcnorm = normalize_gc(gchist.items(), total)
        return (gchist, gcnorm)
    
    def get_read_plot(self, context, data, read):
        plot_data = dict((sample_id, values[0]) for sample_id, values in data.items())
        plot_data_norm = dict((sample_id, values[1]) for sample_id, values in data.items())
        plot_config = {
            'id': 'atropos_gc_plot_{phase}_{read}'.format(
                phase=context['phase'], read=read),
            'title': 'Per Sequence GC Content: Read {}'.format(read),
            'ylab': 'Count',
            'xlab': '%GC',
            'ymin': 0,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            'colors': get_status_cols(context['statuses']),
            'data_labels': [
                {'name': 'Percentages', 'ylab': 'Percentage'},
                {'name': 'Counts', 'ylab': 'Count'}
            ]
        }
        esconfig = {
            'name': 'Theoretical GC Content',
            'dashStyle': 'Dash',
            'lineWidth': 2,
            'color': '#000000',
            'marker': { 'enabled': False },
            'enableMouseTracking': False,
            'showInLegend': False,
        }
        plot_config['extra_series'] = [ [dict(esconfig)], [dict(esconfig)] ]
        plot_config['extra_series'][0][0]['data'] = [ 
            (d[0], d[1] * 100)
            for d in context['theoretical_gc'] ]
        plot_config['extra_series'][1][0]['data'] = [ 
            (d[0], d[1] * context['max_total'])
            for d in context['theoretical_gc'] ]
        return linegraph.plot([plot_data_norm, plot_data], plot_config)

class PerBaseN(QcSection):
    name = 'per_base_N_content'
    display = 'Per Base N Content'
    anchor = 'atropos_per_base_N_content'
    compare = operator.gt
    threshold_statistic = 'frac_N'
    default_thresholds = (0.2, 0.05)
    html = """
<p>The percentage of base calls at each position for which an N was called.
See the <a href="{}" target="_blank">Atropos help</a>.</p>
{{plot}}""".format(ATROPOS_DOC_URL)
    plot_type = linegraph
    
    def prepare_data(self, data, general):
        """Returns the data for this section from the parent dict.
        """
        def get_n_data(bases):
            nidx = bases['columns'].index('N')
            return dict(
                n_counts = ordered_dict(
                    (int(pos), bases['rows'][pos][nidx])
                    for pos in bases['rows'].keys()),
                total_counts = ordered_dict(
                    (int(pos), sum(bases['rows'][pos]))
                    for pos in bases['rows'].keys()))
        
        section_data = {}
        for sample_id, (read1_data, read2_data) in data.items():
            sample_data = [None, None]
            if read1_data:
                sample_data[0] = get_n_data(read1_data['bases'])
            if read2_data:
                sample_data[1] = get_n_data(read2_data['bases'])
            section_data[sample_id] = sample_data
        return section_data
    
    def compute_statistic(self, context, data):
        return max(
            data['n_counts'][pos] / data['total_counts'][pos]
            for pos in data['total_counts'].keys())

    def get_plot_config(self, context, read=1):
        return {
            'id': 'atropos_per_base_N_content_plot_{phase}_{read}'.format(
                phase=context['phase'], read=read),
            'title': 'Per Base N Content: Read {}'.format(read),
            'ylab': 'Percentage N-Count',
            'xlab': 'Position in Read (bp)',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': get_status_cols(context['statuses']),
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ]
        }
    
    def get_sample_plot_data(self, context, data):
        """Create the HTML for the per base N content plot.
        """
        return ordered_dict(
            (key, data['n_counts'].get(key, 0) * 100 / data['total_counts'].get(key))
            for key in data['total_counts'].keys())

class SequenceLength(QcSection):
    name = 'lengths'
    display = 'Sequence Length Distribution'
    anchor = 'atropos_lengths'
    threshold_statistic = 'length_range'
    html = """
<p>The distribution of fragment sizes (read lengths) found.
See the <a href="{}" target="_blank">Atropos help</a>.</p>
{{plot}}""".format(ATROPOS_DOC_URL)
    all_same_html = """
<p>All samples have sequences of a single length ({length} bp).</p>"""
    all_same_within_samples_html = """
<p>All samples have sequences of a single length ({length} bp).
See the <a href="#general_stats">General Statistics Table</a>.</p>"""
    
    def get_status_for(self, context, stat):
        """This always returns 'pass' for post-trimming.
        """
        if context['phase'] == 'pre':
            min_len, max_len = stat
            if min_len == 0:
                return FAIL
            if min_len != max_len:
                return WARN
        return PASS
    
    def get_html(self, context, data):
        unique_lengths = set()
        multiple_lengths = False
        
        def _handle_read(read_data, unique_lengths):
            unique_lengths.update(read_data['hist'].keys())
            return len(read_data['hist']) > 1
        
        for sample_id, (data1, data2) in data.items():
            multiple_lengths |= _handle_read(data1, unique_lengths)
            multiple_lengths |= _handle_read(data2, unique_lengths)
        
        if not multiple_lengths:
            if len(unique_lengths) == 1:
                html = self.all_same_html
            else:
                html = self.all_same_within_samples_html
            return html.format(",".join(unique_lengths))
        else:
            return super(SequenceLength, self).get_html(context, data)
    
    def get_plot_config(self, context, read=1):
        return {
            'id': 'atropos_lengths_plot_{phase}_{read}'.format(
                phase=context['phase'], read=read),
            'title': 'Sequence Length Distribution: Read {}'.format(read),
            'ylab': 'Read Count',
            'xlab': 'Sequence Length (bp)',
            'ymin': 0,
            'yMinTickInterval': 0.1,
            'xDecimals': False,
            'colors': get_status_cols(context['statuses']),
            'tt_label': '<b>{point.x} bp</b>: {point.y}',
        }
    
    def get_sample_plot_data(self, context, data):
        return data['lengths']['hist']

## Utils

class NormalDistribution(object):
    def __init__(self, mean, stdev):
        self.mean = mean
        self.sd2 = 2 * stdev * stdev
    
    def __call__(self, value):
        return (
            (math.e ** -(((value - self.mean) ** 2) / self.sd2)) /
            math.sqrt(self.sd2 * math.pi))

def ordered_dict(iterable):
    d = OrderedDict()
    for key, value in iterable:
        d[key] = value
    return d

def sort_hist(hist):
    return ordered_dict(sorted(
        ((int(length), count) for length, count in hist.items()),
        key=lambda x: x[0]))

def weighted_lower_quantile(values, freqs, quantile):
    for i, frac in enumerate(freqs):
        if frac > quantile:
            return values[i-1] if i > 0 else values[0]
    return values[-1]

def weighted_mean(vals, counts):
    return sum(multiply(vals, counts)) / sum(counts)

def weighted_median(vals, counts_cumsum):
    total = counts_cumsum[-1]
    if total == 0:
        return None
    mid1 = mid2 = (total // 2) + 1
    if total % 2 == 0:
        mid1 -= 1
    val1 = val2 = None
    for i, val in enumerate(counts_cumsum):
        if val1 is None and mid1 <= val:
            val1 = vals[i]
        if mid2 <= val:
            val2 = vals[i]
            break
    return float(val1 + val2) / 2

def mode(data):
    """Given a sequence of (value, count) tuples, returns the value associated
    with the largest count.
    """
    return tuple(sorted(data, key=lambda x: x[1]))[-1][0]

def hist_to_means(hist):
    """Get weighted means for each column in a position histogram.
    """
    return dict(
        (int(pos), weighted_mean(hist['columns'], counts))
        for pos, counts in hist['rows'].items())

def get_status_cols(statuses):
    """Helper function - returns a list of colours according to the
    status of this module for each sample.
    """
    return ordered_dict(
        (sample_id, status.color)
        for sample_id, status in statuses.items())

def simplify(obj):
    """Simplify an object for json serialization.
    """
    if (
            obj is None or
            isinstance(obj, str) or 
            isinstance(obj, Number) or 
            isinstance(obj, bool)):
        return obj
    
    if isinstance(obj, dict):
        newobj = {}
        for key, value in obj.items():
            newobj[str(key)] = simplify(value)
        return newobj
    
    if isinstance(obj, Iterable):
        return tuple(simplify(item) for item in obj)
    
    return str(obj)

    # These stats are not yet implemented in Atropos.
    #
    # def seq_dup_levels_plot (self):
    #     """ Create the HTML for the Sequence Duplication Levels plot """
    #
    #     data = dict()
    #     for s_name in self.fastqc_data:
    #         try:
    #             d = {d['duplication_level']: d['percentage_of_total'] for d in self.fastqc_data[s_name]['sequence_duplication_levels']}
    #             data[s_name] = OrderedDict()
    #             for k in self.dup_keys:
    #                 try:
    #                     data[s_name][k] = d[k]
    #                 except KeyError:
    #                     pass
    #         except KeyError:
    #             pass
    #     if len(data) == 0:
    #         log.debug('lengths not found in FastQC reports')
    #         return None
    #
    #     plot_config = {
    #         'id': 'fastqc_sequence_duplication_levels_plot',
    #         'title': 'Sequence Duplication Levels',
    #         'categories': True,
    #         'ylab': '% of Library',
    #         'xlab': 'Sequence Duplication Level',
    #         'ymax': 100,
    #         'ymin': 0,
    #         'yMinTickInterval': 0.1,
    #         'colors': self.get_status_cols('sequence_duplication_levels'),
    #         'tt_label': '<b>{point.x}</b>: {point.y:.1f}%',
    #     }
    #
    #     self.sections.append({
    #         'name': 'Sequence Duplication Levels',
    #         'anchor': 'fastqc_sequence_duplication_levels',
    #         'content': '<p>The relative level of duplication found for every sequence. ' +
    #                     'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html" target="_blank">FastQC help</a>.</p>' +
    #                     linegraph.plot(data, plot_config)
    #     })
    #
    # def overrepresented_sequences (self):
    #     """Sum the percentages of overrepresented sequences and display them in a bar plot"""
    #
    #     data = dict()
    #     for s_name in self.fastqc_data:
    #         data[s_name] = dict()
    #         try:
    #             max_pcnt   = max( [ float(d['percentage']) for d in self.fastqc_data[s_name]['overrepresented_sequences']] )
    #             total_pcnt = sum( [ float(d['percentage']) for d in self.fastqc_data[s_name]['overrepresented_sequences']] )
    #             data[s_name]['total_overrepresented'] = total_pcnt
    #             data[s_name]['top_overrepresented'] = max_pcnt
    #             data[s_name]['remaining_overrepresented'] = total_pcnt - max_pcnt
    #         except KeyError:
    #             if self.fastqc_data[s_name]['statuses']['overrepresented_sequences'] == 'pass':
    #                 data[s_name]['total_overrepresented'] = 0
    #                 data[s_name]['top_overrepresented'] = 0
    #                 data[s_name]['remaining_overrepresented'] = 0
    #             else:
    #                 log.debug("Couldn't find data for {}, invalid Key".format(s_name))
    #
    #     cats = OrderedDict()
    #     cats['top_overrepresented'] = { 'name': 'Top over-represented sequence' }
    #     cats['remaining_overrepresented'] = { 'name': 'Sum of remaining over-represented sequences' }
    #
    #     # Config for the plot
    #     plot_config = {
    #         'id': 'fastqc_overrepresented_sequencesi_plot',
    #         'title': 'Overrepresented sequences',
    #         'ymin': 0,
    #         'yCeiling': 100,
    #         'yMinRange': 20,
    #         'tt_decimals': 2,
    #         'tt_suffix': '%',
    #         'tt_percentages': False,
    #         'ylab_format': '{value}%',
    #         'cpswitch': False,
    #         'ylab': 'Percentage of Total Sequences'
    #     }
    #
    #     # Check if any samples have more than 1% overrepresented sequences, else don't make plot.
    #     if max([ x['total_overrepresented'] for x in data.values()]) < 1:
    #         plot_html = '<div class="alert alert-info">{} samples had less than 1% of reads made up of overrepresented sequences</div>'.format(len(data))
    #     else:
    #         plot_html = bargraph.plot(data, cats, plot_config)
    #
    #     self.sections.append({
    #         'name': 'Overrepresented sequences',
    #         'anchor': 'fastqc_overrepresented_sequences',
    #         'content': '<p> The total amount of overrepresented sequences found in each library. ' +
    #                 'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html" target="_blank">FastQC help for further information</a>.</p>'
    #                 + plot_html
    #         })
    #
    # def adapter_content_plot (self):
    #     """ Create the HTML for the FastQC adapter plot """
    #
    #     data = dict()
    #     for s_name in self.fastqc_data:
    #         try:
    #             for d in self.fastqc_data[s_name]['adapter_content']:
    #                 pos = self.avg_bp_from_range(d['position'])
    #                 for r in self.fastqc_data[s_name]['adapter_content']:
    #                     pos = self.avg_bp_from_range(r['position'])
    #                     for a in r.keys():
    #                         k = "{} - {}".format(s_name, a)
    #                         if a != 'position':
    #                             try:
    #                                 data[k][pos] = r[a]
    #                             except KeyError:
    #                                 data[k] = {pos: r[a]}
    #         except KeyError:
    #             pass
    #     if len(data) == 0:
    #         log.debug('adapter_content not found in FastQC reports')
    #         return None
    #
    #     # Lots of these datasets will be all zeros.
    #     # Only take datasets with > 0.1% adapter contamination
    #     data = {k:d for k,d in data.items() if max(data[k].values()) >= 0.1 }
    #
    #     plot_config = {
    #         'id': 'fastqc_adapter_content_plot',
    #         'title': 'Adapter Content',
    #         'ylab': '% of Sequences',
    #         'xlab': 'Position',
    #         'yCeiling': 100,
    #         'yMinRange': 5,
    #         'ymin': 0,
    #         'xDecimals': False,
    #         'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
    #         'hide_empty': True,
    #         'yPlotBands': [
    #             {'from': 20, 'to': 100, 'color': '#e6c3c3'},
    #             {'from': 5, 'to': 20, 'color': '#e6dcc3'},
    #             {'from': 0, 'to': 5, 'color': '#c3e6c3'},
    #         ],
    #     }
    #
    #     if len(data) > 0:
    #         plot_html = linegraph.plot(data, plot_config)
    #     else:
    #         plot_html = '<div class="alert alert-warning">No samples found with any adapter contamination > 0.1%</div>'
    #
    #     # Note - colours are messy as we've added adapter names here. Not
    #     # possible to break down pass / warn / fail for each adapter, which
    #     # could lead to misleading labelling (fails on adapter types with
    #     # little or no contamination)
    #
    #     self.sections.append({
    #         'name': 'Adapter Content',
    #         'anchor': 'fastqc_adapter_content',
    #         'content': '<p>The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position. ' +
    #                     'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html" target="_blank">FastQC help</a>. ' +
    #                     'Only samples with &ge; 0.1% adapter contamination are shown.</p>' + plot_html
    #     })
