import json
import logging
import operator
import os
from collections import OrderedDict, defaultdict
from future.utils import iteritems
from itertools import islice

from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='bcl2fastq',
            anchor='bcl2fastq',
            href="https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html",
            info="can be used to both demultiplex data and convert BCL files"
                 " to FASTQ file formats for downstream analysis."
        )

        # Gather data from all json files
        self.bcl2fastq_data = dict()
        for myfile in self.find_log_files('bcl2fastq'):
            self.parse_file_as_json(myfile)

        # Collect counts by lane and sample (+source_files)
        self.bcl2fastq_bylane = dict()
        self.bcl2fastq_bysample = dict()
        self.bcl2fastq_bysample_lane = dict()
        self.source_files = dict()
        self.split_data_by_lane_and_sample()

        # Filter to strip out ignored sample names
        self.bcl2fastq_bylane = self.ignore_samples(self.bcl2fastq_bylane)
        self.bcl2fastq_bysample = self.ignore_samples(self.bcl2fastq_bysample)
        self.bcl2fastq_bysample_lane = self.ignore_samples(self.bcl2fastq_bysample_lane)

        # Return with Warning if no files are found
        if len(self.bcl2fastq_bylane) == 0 and len(self.bcl2fastq_bysample) == 0:
            raise UserWarning

        # Print source files
        for s in self.source_files.keys():
            self.add_data_source(
                s_name=s,
                source=",".join(list(set(self.source_files[s]))),
                module='bcl2fastq',
                section='bcl2fastq-bysample'
            )

        # Add sample counts to general stats table
        self.add_general_stats()
        self.write_data_file(
            {str(k): self.bcl2fastq_bylane[k] for k in self.bcl2fastq_bylane.keys()},
            'multiqc_bcl2fastq_bylane'
        )
        self.write_data_file(self.bcl2fastq_bysample, 'multiqc_bcl2fastq_bysample')

        # Add section for summary stats per flow cell
        self.add_section (
            name = 'Lane Statistics',
            anchor = 'bcl2fastq-lanestats',
            description = 'Statistics about each lane for each flowcell',
            plot = self.lane_stats_table()
        )

        # Add section for counts by lane
        cats = OrderedDict()
        cats["perfect"] = {'name': 'Perfect Index Reads'}
        cats["imperfect"] = {'name': 'Mismatched Index Reads'}
        cats["undetermined"] = {'name': 'Undetermined Reads'}
        self.add_section (
            name = 'Clusters by lane',
            anchor = 'bcl2fastq-bylane',
            description = 'Number of reads per lane (with number of perfect index reads).',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.""",
            plot = bargraph.plot(
                self.get_bar_data_from_counts(self.bcl2fastq_bylane),
                cats,
                {
                    'id': 'bcl2fastq_lane_counts',
                    'title': 'bcl2fastq: Clusters by lane',
                    'ylab': 'Number of clusters',
                    'hide_zero_cats': False
                }
            )
        )

        # Add section for counts by sample
        # get cats for per-lane tab
        lcats = set()
        for s_name in self.bcl2fastq_bysample_lane:
            lcats.update(self.bcl2fastq_bysample_lane[s_name].keys())
        lcats = sorted(list(lcats))
        self.add_section (
            name='Clusters by sample',
            anchor = 'bcl2fastq-bysample',
            description = 'Number of reads per sample.',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combinned. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot = bargraph.plot(
                [
                    self.get_bar_data_from_counts(self.bcl2fastq_bysample),
                    self.bcl2fastq_bysample_lane
                ],
                [cats, lcats],
                {
                    'id': 'bcl2fastq_sample_counts',
                    'title': 'bcl2fastq: Clusters by sample',
                    'hide_zero_cats': False,
                    'ylab': 'Number of clusters',
                    'data_labels': ['Index mismatches', 'Counts per lane']
                }
            )
        )

        # Add section with undetermine barcodes
        self.add_section(
            name="Undetermined barcodes by lane",
            anchor="undetermine_by_lane",
            description="Count of the top twenty most abundant undetermined"
                        " barcodes by lanes",
            plot=bargraph.plot(
                self.get_bar_data_from_undetermined(self.bcl2fastq_bylane),
                None,
                {
                    'id': 'bcl2fastq_undetermined',
                    'title': 'bcl2fastq: Undetermined barcodes by lane',
                    'ylab': 'Count',
                    'cpswitch': False,
                    'tt_percentages': False,
                    'use_legend': True,
                    'tt_suffix': 'reads'
                }
            )
        )

    def parse_file_as_json(self, myfile):
        try:
            content = json.loads(myfile["f"])
        except ValueError:
            log.warn('Could not parse file as json: {}'.format(myfile["fn"]))
            return
        runId = content["RunId"]
        if runId not in self.bcl2fastq_data:
            self.bcl2fastq_data[runId] = dict()
        run_data = self.bcl2fastq_data[runId]
        for conversionResult in content.get("ConversionResults", []):
            l = conversionResult["LaneNumber"]
            lane = 'L{}'.format(l)
            if lane in run_data:
                log.debug(
                    "Duplicate runId/lane combination found! Overwriting: {}".format(
                        self.prepend_runid(runId, lane)
                    )
                )
            run_data[lane] = {
                "total": 0,
                "total_yield": 0,
                "perfectIndex": 0,
                "samples": dict(),
                "yieldQ30": 0,
                "qscore_sum": 0
            }
            # simplify the population of dictionnaries
            rlane = run_data[lane]
            # Add undetermine barcodes
            try:
                unknown_barcode = content['UnknownBarcodes'][l - 1]['Barcodes']
            except IndexError:
                unknown_barcode = next(
                    (item['Barcodes'] for item in content['UnknownBarcodes']
                     if item['Lane'] == 8),
                    None
                )

            run_data[lane]['unknown_barcodes'] = unknown_barcode
            for demuxResult in conversionResult.get("DemuxResults", []):
                sample = "{}-{}".format(demuxResult["SampleName"],
                                        demuxResult["SampleId"])
                if sample in run_data[lane]["samples"]:
                    log.debug(
                        "Duplicate runId/lane/sample combination found! Overwriting: {}, {}".format(
                            self.prepend_runid(runId, lane), sample
                        )
                    )
                run_data[lane]["samples"][sample] = {
                    "total": 0,
                    "total_yield": 0,
                    "perfectIndex": 0,
                    "filename": os.path.join(myfile['root'], myfile["fn"]),
                    "yieldQ30": 0,
                    "qscore_sum": 0,
                    "R1_yield": 0,
                    "R2_yield": 0,
                    "R1_Q30": 0,
                    "R2_Q30": 0,
                    "R1_trimmed_bases": 0,
                    "R2_trimmed_bases": 0
                }
                # simplify the population of dictionnaries
                lsample = run_data[lane]["samples"][sample]
                rlane["total"] += demuxResult["NumberReads"]
                rlane["total_yield"] += demuxResult["Yield"]
                lsample["total"] += demuxResult["NumberReads"]
                lsample["total_yield"] += demuxResult["Yield"]
                for indexMetric in demuxResult.get("IndexMetrics", []):
                    rlane["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                    lsample["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                for readMetric in demuxResult.get("ReadMetrics", []):
                    r = readMetric["ReadNumber"]
                    rlane["yieldQ30"] += readMetric["YieldQ30"]
                    rlane["qscore_sum"] += readMetric["QualityScoreSum"]
                    lsample["yieldQ30"] += readMetric["YieldQ30"]
                    lsample["qscore_sum"] += readMetric["QualityScoreSum"]
                    lsample["R{}_yield".format(r)] += readMetric["Yield"]
                    lsample["R{}_Q30".format(r)] += readMetric["YieldQ30"]
                    lsample["R{}_trimmed_bases".format(r)] += readMetric["TrimmedBases"]
            undeterminedYieldQ30 = 0
            undeterminedQscoreSum = 0
            undeterminedTrimmedBases = 0
            if "Undetermined" in conversionResult:
                for readMetric in conversionResult["Undetermined"]["ReadMetrics"]:
                    undeterminedYieldQ30 += readMetric["YieldQ30"]
                    undeterminedQscoreSum += readMetric["QualityScoreSum"]
                    undeterminedTrimmedBases += readMetric["TrimmedBases"]
                run_data[lane]["samples"]["undetermined"] = {
                    "total": conversionResult["Undetermined"]["NumberReads"],
                    "total_yield": conversionResult["Undetermined"]["Yield"],
                    "perfectIndex": 0,
                    "yieldQ30": undeterminedYieldQ30,
                    "qscore_sum": undeterminedQscoreSum,
                    "trimmed_bases": undeterminedTrimmedBases
                }

        # Calculate Percents and averages
        for lane_id, lane in iteritems(run_data):
            try:
                lane["percent_Q30"] = (float(lane["yieldQ30"])
                    / float(lane["total_yield"])) * 100.0
            except ZeroDivisionError:
                lane["percent_Q30"] = "NA"
            try:
                lane["percent_perfectIndex"] = (float(lane["perfectIndex"])
                    / float(lane["total"])) * 100.0
            except ZeroDivisionError:
                lane["percent_perfectIndex"] = "NA"
            try:
                lane["mean_qscore"] = float(lane["qscore_sum"]) / float(lane["total_yield"])
            except ZeroDivisionError:
                lane["mean_qscore"] = "NA"
            for sample_id, sample in iteritems(lane["samples"]):
                try:
                    sample["percent_Q30"] = (float(sample["yieldQ30"])
                        / float(sample["total_yield"])) * 100.0
                except ZeroDivisionError:
                    sample["percent_Q30"] = "NA"
                try:
                    sample["percent_perfectIndex"] = (float(sample["perfectIndex"])
                        / float(sample["total"])) * 100.0
                except ZeroDivisionError:
                    sample["percent_perfectIndex"] = "NA"
                try:
                    sample["mean_qscore"] = float(sample["qscore_sum"]) / float(sample["total_yield"])
                except ZeroDivisionError:
                    sample["mean_qscore"] = "NA"

    def split_data_by_lane_and_sample(self):
        for run_id, r in iteritems(self.bcl2fastq_data):
            for lane_id, lane in iteritems(r):
                uniqLaneName = self.prepend_runid(run_id, lane_id)
                self.bcl2fastq_bylane[uniqLaneName] = {
                    "total": lane["total"],
                    "total_yield": lane["total_yield"],
                    "perfectIndex": lane["perfectIndex"],
                    "undetermined": lane["samples"].get("undetermined", {}).get("total", "NA"),
                    "yieldQ30": lane["yieldQ30"],
                    "qscore_sum": lane["qscore_sum"],
                    "percent_Q30": lane["percent_Q30"],
                    "percent_perfectIndex": lane["percent_perfectIndex"],
                    "mean_qscore": lane["mean_qscore"],
                    "unknown_barcodes": self.get_unknown_barcodes(lane['unknown_barcodes']),
                }
                for sample_id, sample in iteritems(lane["samples"]):
                    if sample_id not in self.bcl2fastq_bysample:
                        self.bcl2fastq_bysample[sample_id] = {
                            "total": 0,
                            "total_yield": 0,
                            "R1_yield": 0,
                            "R2_yield": 0,
                            "perfectIndex": 0,
                            "yieldQ30": 0,
                            "R1_Q30": 0,
                            "R2_Q30": 0,
                            "R1_trimmed_bases": 0,
                            "R2_trimmed_bases": 0,
                            "qscore_sum": 0
                        }
                    s = self.bcl2fastq_bysample[sample_id]
                    s["total"] += sample["total"]
                    s["total_yield"] += sample["total_yield"]
                    s["perfectIndex"] += sample["perfectIndex"]
                    s["yieldQ30"] += sample["yieldQ30"]
                    s["qscore_sum"] += sample["qscore_sum"]
                    # Undetermined samples did not have R1 and R2 information
                    try:
                        s["R1_yield"] += sample["R1_yield"]
                        s["R2_yield"] += sample["R2_yield"]
                        s["R1_Q30"] += sample["R1_Q30"]
                        s["R2_Q30"] += sample["R2_Q30"]
                        s["R1_trimmed_bases"] += sample["R1_trimmed_bases"]
                        s["R2_trimmed_bases"] += sample["R2_trimmed_bases"]
                    except KeyError:
                        pass
                    try:
                        s["percent_Q30"] = (float(s["yieldQ30"])
                            / float(s["total_yield"])) * 100.0
                    except ZeroDivisionError:
                        s["percent_Q30"] = "NA"
                    try:
                        s["percent_perfectIndex"] = (float(s["perfectIndex"])
                            / float(s["total"])) * 100.0
                    except ZeroDivisionError:
                        s["percent_perfectIndex"] = "NA"
                    try:
                        s["mean_qscore"] = (float(s["qscore_sum"]) / float(s["total_yield"]))
                    except ZeroDivisionError:
                        s["mean_qscore"] = "NA"
                    if sample_id != "undetermined":
                        if sample_id not in self.source_files:
                            self.source_files[sample_id] = []
                        self.source_files[sample_id].append(sample["filename"])

    def get_unknown_barcodes(self, lane_unknown_barcode):
        """ Python3.6 dictionnaries keep the order but not in other version.
        This function return an `OrderedDict` sorted by barcode count.
        """
        try:
            sorted_barcodes = OrderedDict(
                sorted(
                    iteritems(lane_unknown_barcode),
                    key=operator.itemgetter(1),
                    reverse=True
                )
            )
        except AttributeError:
            sorted_barcodes=None
        return sorted_barcodes

    def add_general_stats(self):
        data = dict()
        for sample_id, sample in iteritems(self.bcl2fastq_bysample):
            # Zero division is possible
            try:
                percent_R1_Q30 = '{0:.1f}'.format(
                    float(100.0 * sample["R1_Q30"] / sample["R1_yield"])
                )
            except ZeroDivisionError:
                percent_R1_Q30 = '0.0'
            try:
                percent_R2_Q30 = '{0:.1f}'.format(
                    float(100.0 * sample["R2_Q30"] / sample["R2_yield"])
                )
            except ZeroDivisionError:
                percent_R2_Q30 = '0.0'
            try:
                perfect_percent = '{0:.1f}'.format(
                    float(100.0 * sample["perfectIndex"] / sample["total"])
                )
            except ZeroDivisionError:
                perfect_percent = '0.0'

            data[sample_id] = {
                "yieldQ30": sample["yieldQ30"],
                "percent_R1_Q30": percent_R1_Q30,
                "percent_R2_Q30": percent_R2_Q30,
                "total": sample["total"],
                "perfectPercent": perfect_percent,
                "R1_trimmed_bases": sample["R1_trimmed_bases"],
                "R2_trimmed_bases": sample["R2_trimmed_bases"]
            }

        headers = OrderedDict()
        headers['total'] = {
            'title': '{} Clusters'.format(config.read_count_prefix),
            'description': 'Total number of reads for this sample as determined by bcl2fastq demultiplexing ({})'.format(config.read_count_desc),
            'scale': 'RdYlGn',
            'shared_key': 'read_count'
        }
        headers['yieldQ30'] = {
            'title': '{} Yield &ge; Q30'.format(config.base_count_prefix),
            'description': 'Number of bases with a Phred score of 30 or higher ({})'.format(config.base_count_desc),
            'scale': 'RdYlGn',
            'shared_key': 'base_count'
        }
        headers['percent_R1_Q30'] = {
            'title': '% R1 Yield &ge; Q30',
            'description': 'Percent of bases in R1 with a Phred score of 30 or higher',
            'scale': 'RdYlGn',
            'max': 100,
            'min': 0,
            'suffix': '%'
        }
        headers['percent_R2_Q30'] = {
            'title': '% R2 Yield &ge; Q30',
            'description': 'Percent of bases in R2 with a Phred score of 30 or higher',
            'scale': 'RdYlGn',
            'max': 100,
            'min': 0,
            'suffix': '%'
        }
        headers['perfectPercent'] = {
            'title': '% Perfect Index',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        headers['R1_trimmed_bases'] = {
            'title': '{} R1 trimmed'.format(config.base_count_prefix),
            'description': 'Number of bases trimmed ({})'.format(config.base_count_desc),
            'scale': 'RdYlGn',
            'modify': lambda x: x * 0.000001
        }
        headers['R2_trimmed_bases'] = {
            'title': '{} R2 trimmed'.format(config.base_count_prefix),
            'description': 'Number of bases trimmed ({})'.format(config.base_count_desc),
            'scale': 'RdYlGn',
            'modify': lambda x: x * 0.000001
        }
        self.general_stats_addcols(data, headers)

    def lane_stats_table(self):
        """ Return a table with overview stats for each bcl2fastq lane for a single flow cell """
        headers = OrderedDict()
        headers['total_yield'] = {
            'title': '{} Total Yield'.format(config.base_count_prefix),
            'description': 'Number of bases ({})'.format(config.base_count_desc),
            'scale': 'RdYlGn',
            'shared_key': 'base_count'
        }
        headers['total'] = {
            'title': '{} Total Clusters'.format(config.read_count_prefix),
            'description': 'Total number of clusters for this lane ({})'.format(config.read_count_desc),
            'scale': 'RdYlGn',
            'shared_key': 'read_count'
        }
        headers['percent_Q30'] = {
            'title': '% bases &ge; Q30',
            'description': 'Percentage of bases with greater than or equal to Q30 quality score',
            'suffix': '%',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn'
        }
        headers['mean_qscore'] = {
            'title': 'Mean Quality',
            'description': 'Average phred qualty score',
            'min': 0,
            'scale': 'RdYlGn'
        }
        headers['percent_perfectIndex'] = {
            'title': '% Perfect Index',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        table_config = {
            'namespace': 'bcl2fastq',
            'id': 'bcl2fastq-lane-stats-table',
            'table_title': 'bcl2fastq Lane Statistics',
            'col1_header': 'Run ID - Lane',
            'no_beeswarm': True
        }
        return table.plot(self.bcl2fastq_bylane, headers, table_config)

    def prepend_runid(self, runId, rest):
        return str(runId)+" - "+str(rest)

    def get_bar_data_from_counts(self, counts):
        bar_data = {}
        for key, value in iteritems(counts):
            bar_data[key] = {
                "perfect": value["perfectIndex"],
                "imperfect": value["total"] - value["perfectIndex"],
            }
            if "undetermined" in value:
                bar_data[key]["undetermined"] = value["undetermined"]
        return bar_data

    def get_bar_data_from_undetermined(self, flowcells):
        """ Get data to plot for undetermined barcodes.
        """
        bar_data = defaultdict(dict)
        # get undetermined barcodes for each lanes
        for lane_id, lane in iteritems(flowcells):
            try:
                for barcode, count in islice(iteritems(lane['unknown_barcodes']), 20):
                    bar_data[barcode][lane_id] = count
            except AttributeError:
                pass

        # sort results
        bar_data = OrderedDict(sorted(
            iteritems(bar_data),
            key=lambda x: sum(x[1].values()),
            reverse=True
        ))
        return OrderedDict(
            (key, value) for key, value in islice(iteritems(bar_data), 20)
        )
