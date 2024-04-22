import json
import logging
import os
from collections import defaultdict
from itertools import islice

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="bcl2fastq",
            anchor="bcl2fastq",
            href="https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html",
            info="can be used to both demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.",
            # Can't find a DOI // doi=
        )

        # Gather data from all json files
        self.bcl2fastq_data = dict()
        for f in self.find_log_files("bcl2fastq"):
            self.parse_file_as_json(f)

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
            raise ModuleNoSamplesFound

        # Print source files
        for s in self.source_files.keys():
            self.add_data_source(
                s_name=s,
                source=",".join(list(set(self.source_files[s]))),
                module="bcl2fastq",
                section="bcl2fastq-bysample",
            )

        # Add sample counts to general stats table
        self.add_general_stats()
        self.write_data_file(
            {str(k): self.bcl2fastq_bylane[k] for k in self.bcl2fastq_bylane.keys()}, "multiqc_bcl2fastq_bylane"
        )
        self.write_data_file(self.bcl2fastq_bysample, "multiqc_bcl2fastq_bysample")

        # Add section for summary stats per flow cell
        self.add_section(
            name="Lane Statistics",
            anchor="bcl2fastq-lanestats",
            description="Statistics about each lane for each flowcell",
            plot=self.lane_stats_table(),
        )

        # Add section for counts by lane
        cats = {
            "perfect": {"name": "Perfect Index Reads"},
            "imperfect": {"name": "Mismatched Index Reads"},
            "undetermined": {"name": "Undetermined Reads"},
        }
        self.add_section(
            name="Clusters by lane",
            anchor="bcl2fastq-bylane",
            description="Number of reads per lane (with number of perfect index reads).",
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.""",
            plot=bargraph.plot(
                self.get_bar_data_from_counts(self.bcl2fastq_bylane),
                cats,
                {
                    "id": "bcl2fastq_lane_counts",
                    "title": "bcl2fastq: Clusters by lane",
                    "ylab": "Number of clusters",
                    "hide_zero_cats": False,
                },
            ),
        )

        # Add section for counts by sample
        # get cats for per-lane tab
        lcats = set()
        for s_name in self.bcl2fastq_bysample_lane:
            lcats.update(self.bcl2fastq_bysample_lane[s_name].keys())
        lcats = sorted(list(lcats))
        self.add_section(
            name="Clusters by sample",
            anchor="bcl2fastq-bysample",
            description="Number of reads per sample.",
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combined. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot=bargraph.plot(
                [self.get_bar_data_from_counts(self.bcl2fastq_bysample), self.bcl2fastq_bysample_lane],
                [cats, lcats],
                {
                    "id": "bcl2fastq_sample_counts",
                    "title": "bcl2fastq: Clusters by sample",
                    "hide_zero_cats": False,
                    "ylab": "Number of clusters",
                    "data_labels": ["Index mismatches", "Counts per lane"],
                },
            ),
        )

        # Add section with undetermined barcodes
        self.add_section(
            name="Undetermined barcodes by lane",
            anchor="undetermined_by_lane",
            description="Count of the top twenty most abundant undetermined barcodes by lanes",
            plot=bargraph.plot(
                self.get_bar_data_from_undetermined(self.bcl2fastq_bylane),
                None,
                {
                    "id": "bcl2fastq_undetermined",
                    "title": "bcl2fastq: Undetermined barcodes by lane",
                    "ylab": "Reads",
                    "use_legend": True,
                    "sort_samples": False,  # keep top barcode on top
                },
            ),
        )

    def parse_file_as_json(self, f):
        try:
            content = json.loads(f["f"])
        except ValueError:
            log.warning(f"Could not parse file as json: {f['fn']}")
            return

        # Clean / prepend directories to sample names
        runId = self.clean_s_name(content["RunId"], f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None, runId)

        if runId not in self.bcl2fastq_data:
            self.bcl2fastq_data[runId] = dict()
        run_data = self.bcl2fastq_data[runId]
        for conversionResult in content.get("ConversionResults", []):
            lane_num = conversionResult["LaneNumber"]
            lane = f"L{conversionResult['LaneNumber']}"
            if lane in run_data:
                log.debug(f"Duplicate runId/lane combination found! Overwriting: {self.prepend_runid(runId, lane)}")
            run_data[lane] = {
                "total": 0,
                "total_yield": 0,
                "perfectIndex": 0,
                "samples": dict(),
                "yieldQ30": 0,
                "qscore_sum": 0,
            }
            # simplify the population of dictionaries
            rlane = run_data[lane]

            # Add undetermined barcodes
            unknown_barcode = dict()
            for lane_data in content.get("UnknownBarcodes", list()):
                if lane_data["Lane"] == lane_num:
                    unknown_barcode = lane_data["Barcodes"]
                    break
            run_data[lane]["unknown_barcodes"] = unknown_barcode

            for demuxResult in conversionResult.get("DemuxResults", []):
                if demuxResult["SampleName"] == demuxResult["SampleId"]:
                    sample = demuxResult["SampleName"]
                else:
                    sample = f"{demuxResult['SampleId']}-{demuxResult['SampleName']}"
                sample = self.clean_s_name(sample, f)
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
                    "filename": os.path.join(f["root"], f["fn"]),
                    "yieldQ30": 0,
                    "qscore_sum": 0,
                }
                # simplify the population of dictionaries
                lsample = run_data[lane]["samples"][sample]
                for r in range(1, 5):
                    lsample[f"R{r}_yield"] = 0
                    lsample[f"R{r}_Q30"] = 0
                    lsample[f"R{r}_trimmed_bases"] = 0
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
                    lsample[f"R{r}_yield"] += readMetric["Yield"]
                    lsample[f"R{r}_Q30"] += readMetric["YieldQ30"]
                    lsample[f"R{r}_trimmed_bases"] += readMetric["TrimmedBases"]
                # Remove unpopulated read keys
                for r in range(1, 5):
                    if not lsample[f"R{r}_yield"] and not lsample[f"R{r}_Q30"] and not lsample[f"R{r}_trimmed_bases"]:
                        lsample.pop(f"R{r}_yield")
                        lsample.pop(f"R{r}_Q30")
                        lsample.pop(f"R{r}_trimmed_bases")
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
                    "trimmed_bases": undeterminedTrimmedBases,
                }

        # Calculate Percents and averages
        for lane_id, lane in run_data.items():
            try:
                lane["percent_Q30"] = (float(lane["yieldQ30"]) / float(lane["total_yield"])) * 100.0
            except ZeroDivisionError:
                lane["percent_Q30"] = "NA"
            try:
                lane["percent_perfectIndex"] = (float(lane["perfectIndex"]) / float(lane["total"])) * 100.0
            except ZeroDivisionError:
                lane["percent_perfectIndex"] = "NA"
            try:
                lane["mean_qscore"] = float(lane["qscore_sum"]) / float(lane["total_yield"])
            except ZeroDivisionError:
                lane["mean_qscore"] = "NA"

            for sample_id, sample in lane["samples"].items():
                try:
                    sample["percent_Q30"] = (float(sample["yieldQ30"]) / float(sample["total_yield"])) * 100.0
                except ZeroDivisionError:
                    sample["percent_Q30"] = "NA"
                try:
                    sample["percent_perfectIndex"] = (float(sample["perfectIndex"]) / float(sample["total"])) * 100.0
                except ZeroDivisionError:
                    sample["percent_perfectIndex"] = "NA"
                try:
                    sample["mean_qscore"] = float(sample["qscore_sum"]) / float(sample["total_yield"])
                except ZeroDivisionError:
                    sample["mean_qscore"] = "NA"

    def split_data_by_lane_and_sample(self):
        for run_id, r in self.bcl2fastq_data.items():
            for lane_id, lane in r.items():
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
                    "unknown_barcodes": lane["unknown_barcodes"],
                }
                for sample_id, sample in lane["samples"].items():
                    if sample_id not in self.bcl2fastq_bysample:
                        self.bcl2fastq_bysample[sample_id] = {
                            "total": 0,
                            "total_yield": 0,
                            "perfectIndex": 0,
                            "yieldQ30": 0,
                            "qscore_sum": 0,
                        }
                        for r in range(1, 5):
                            self.bcl2fastq_bysample[sample_id][f"R{r}_yield"] = 0
                            self.bcl2fastq_bysample[sample_id][f"R{r}_Q30"] = 0
                            self.bcl2fastq_bysample[sample_id][f"R{r}_trimmed_bases"] = 0
                    s = self.bcl2fastq_bysample[sample_id]
                    s["total"] += sample["total"]
                    s["total_yield"] += sample["total_yield"]
                    s["perfectIndex"] += sample["perfectIndex"]
                    s["yieldQ30"] += sample["yieldQ30"]
                    s["qscore_sum"] += sample["qscore_sum"]
                    # Undetermined samples did not have R1 and R2 information
                    for r in range(1, 5):
                        try:
                            s[f"R{r}_yield"] += sample[f"R{r}_yield"]
                            s[f"R{r}_Q30"] += sample[f"R{r}_Q30"]
                            s[f"R{r}_trimmed_bases"] += sample[f"R{r}_trimmed_bases"]
                        except KeyError:
                            pass
                    try:
                        s["percent_Q30"] = (float(s["yieldQ30"]) / float(s["total_yield"])) * 100.0
                    except ZeroDivisionError:
                        s["percent_Q30"] = "NA"
                    try:
                        s["percent_perfectIndex"] = (float(s["perfectIndex"]) / float(s["total"])) * 100.0
                    except ZeroDivisionError:
                        s["percent_perfectIndex"] = "NA"
                    try:
                        s["mean_qscore"] = float(s["qscore_sum"]) / float(s["total_yield"])
                    except ZeroDivisionError:
                        s["mean_qscore"] = "NA"
                    if sample_id != "undetermined":
                        if sample_id not in self.source_files:
                            self.source_files[sample_id] = []
                        self.source_files[sample_id].append(sample["filename"])
                # Remove unpopulated read keys
                for sample_id, sample in lane["samples"].items():
                    for r in range(1, 5):
                        try:
                            if (
                                not self.bcl2fastq_bysample[sample_id][f"R{r}_yield"]
                                and not self.bcl2fastq_bysample[sample_id][f"R{r}_Q30"]
                                and not self.bcl2fastq_bysample[sample_id][f"R{r}_trimmed_bases"]
                            ):
                                self.bcl2fastq_bysample[sample_id].pop(f"R{r}_yield")
                                self.bcl2fastq_bysample[sample_id].pop(f"R{r}_Q30")
                                self.bcl2fastq_bysample[sample_id].pop(f"R{r}_trimmed_bases")
                        except KeyError:
                            pass

    def add_general_stats(self):
        data = dict()
        for sample_id, sample in self.bcl2fastq_bysample.items():
            percent_R_Q30 = dict()
            for r in range(1, 5):
                # Zero division is possible
                try:
                    percent_R_Q30[r] = "{0:.1f}".format(
                        float(100.0 * sample["R{}_Q30".format(r)] / sample["R{}_yield".format(r)])
                    )
                except ZeroDivisionError:
                    percent_R_Q30[r] = "0.0"
                except KeyError:
                    pass
            try:
                perfect_percent = f"{float(100.0 * sample['perfectIndex'] / sample['total']):.1f}"
            except ZeroDivisionError:
                perfect_percent = "0.0"

            data[sample_id] = {
                "yieldQ30": sample["yieldQ30"],
                "total": sample["total"],
                "perfectPercent": perfect_percent,
            }
            for r in range(1, 5):
                try:
                    data[sample_id][f"percent_R{r}_Q30"] = percent_R_Q30[r]
                    data[sample_id][f"R{r}_trimmed_bases"] = sample[f"R{r}_trimmed_bases"]
                except KeyError:
                    pass

        headers = {
            "total": {
                "title": f"{config.read_count_prefix} Clusters",
                "description": "Total number of reads for this sample as determined by bcl2fastq demultiplexing ({})".format(
                    config.read_count_desc
                ),
                "scale": "Blues",
                "shared_key": "read_count",
            },
            "yieldQ30": {
                "title": f"Yield ({config.base_count_prefix}) ≥ Q30",
                "description": f"Number of bases with a Phred score of 30 or higher ({config.base_count_desc})",
                "scale": "Greens",
                "shared_key": "base_count",
            },
        }
        # If no data for a column, header will be automatically ignored
        for r in range(1, 5):
            headers[f"percent_R{r}_Q30"] = {
                "title": f"% R{r} Yield ≥ Q30",
                "description": f"Percent of bases in R{r} with a Phred score of 30 or higher",
                "scale": "RdYlGn",
                "max": 100,
                "min": 0,
                "suffix": "%",
            }
        headers["perfectPercent"] = {
            "title": "% Perfect Index",
            "description": "Percent of reads with perfect index (0 mismatches)",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        # If no data for a column, header will be automatically ignored
        for r in range(1, 5):
            hideCol = True
            for s in data:
                try:
                    if data[s][f"R{r}_trimmed_bases"] > 0:
                        hideCol = False
                except KeyError:
                    pass
            try:
                headers[f"R{r}_trimmed_bases"] = {
                    "title": f"{config.base_count_prefix} R{r} trimmed",
                    "description": f"Number of bases trimmed ({config.base_count_desc})",
                    "scale": "RdYlBu",
                    "modify": lambda x: x * 0.000001,
                    "hidden": hideCol,
                }
            except KeyError:
                pass
        self.general_stats_addcols(data, headers)

    def lane_stats_table(self):
        """Return a table with overview stats for each bcl2fastq lane for a single flow cell"""
        headers = {
            "total_yield": {
                "title": f"{config.base_count_prefix} Total Yield",
                "description": f"Number of bases ({config.base_count_desc})",
                "scale": "Greens",
                "shared_key": "base_count",
            },
            "total": {
                "title": f"{config.read_count_prefix} Total Clusters",
                "description": f"Total number of clusters for this lane ({config.read_count_desc})",
                "scale": "Blues",
                "shared_key": "read_count",
            },
            "percent_Q30": {
                "title": "% bases ≥ Q30",
                "description": "Percentage of bases with greater than or equal to Q30 quality score",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
            },
            "mean_qscore": {
                "title": "Mean Quality",
                "description": "Average phred qualty score",
                "min": 0,
                "scale": "Spectral",
            },
            "percent_perfectIndex": {
                "title": "% Perfect Index",
                "description": "Percent of reads with perfect index (0 mismatches)",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "suffix": "%",
            },
        }
        table_config = {
            "namespace": "bcl2fastq",
            "id": "bcl2fastq-lane-stats-table",
            "table_title": "bcl2fastq Lane Statistics",
            "col1_header": "Run ID - Lane",
        }
        return table.plot(self.bcl2fastq_bylane, headers, table_config)

    @staticmethod
    def prepend_runid(runId, rest):
        return str(runId) + " - " + str(rest)

    @staticmethod
    def get_bar_data_from_counts(data_by_flowcell):
        bar_data = {}
        for flowcell_id, data in data_by_flowcell.items():
            bar_data[flowcell_id] = {
                "perfect": data["perfectIndex"],
                "imperfect": data["total"] - data["perfectIndex"],
            }
            if "undetermined" in data:
                bar_data[flowcell_id]["undetermined"] = data["undetermined"]
        return bar_data

    @staticmethod
    def get_bar_data_from_undetermined(data_by_flowcell):
        """Get data to plot for undetermined barcodes."""
        bar_data = defaultdict(dict)
        # get undetermined barcodes for each lanes
        for flowcell_id, data in data_by_flowcell.items():
            try:
                for barcode, count in islice(data["unknown_barcodes"].items(), 20):
                    bar_data[barcode][flowcell_id] = count
            except AttributeError:
                pass

        sorted_items = sorted(bar_data.items(), key=lambda kv: sum(kv[1].values()), reverse=True)
        if len(sorted_items) > 20:
            sorted_items = sorted_items[:20]
        return dict(sorted_items)
