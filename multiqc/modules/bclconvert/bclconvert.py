import csv
import logging
import os
import xml.etree.ElementTree as ET
from collections import defaultdict
import functools
from itertools import islice
from typing import Tuple, Dict, Optional

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This BclConvert module is based on the bcl2fastq multiqc module. It can parse multiple
    bclconvert run outputs as long as they are from the same sequencing run. When doing this,
    the undetermined reads will be 'corrected' and re-calculated (as an unknown read from
    one bclcovnert run might not be truly unknown, but simply from another run).

    #### Calculate estimated depth

    You can specify a genome size in config

    It's often useful to talk about sequencing yield in terms of estimated depth of coverage.
    In order to make MultiQC show the estimated depth for each sample, specify the reference genome/target size in your [MultiQC configuration](http://multiqc.info/docs/#configuring-multiqc):

    ```yaml
    bclconvert:
      genome_size: 3049315783
    ```

    The coverage depth will be estimated as the yield Q30 dvivided by the genome size.

    MultiQC comes with effective genome size presets for Human and Mouse, so you can
    provide the genome build name instead, like this: `genome_size: hg38_genome`. The
    following values are supported: `hg19_genome`, `hg38_genome`, `mm10_genome`.

    #### Add barplots containing undetermined barcodes

    By default, the bar plot of undetermined barcodes is only shown when reporting from a single demultiplexing run.

    If you would like to show it with multiple runs (eg. bclconvert runs are split by lane),
    you can specify following parameter in your MultiQC config:

    ```yaml
    bclconvert:
      create_undetermined_barcode_barplots: True
    ```

    The default of this configuration value is `False`
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="BCL Convert",
            anchor="bclconvert",
            href="https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html",
            info="Demultiplexes data and converts BCL files to FASTQ file formats for downstream analysis.",
            # Can't find a DOI // doi=
        )

        # Set up and collate bclconvert run and demux files
        demuxes, qmetrics, multiple_sequencing_runs, last_run_id = self._collate_log_files()

        # variables to store reads for undetermined read recalculation
        self.per_lane_undetermined_reads: Optional[Dict] = dict()
        self.total_reads_in_lane_per_file: Dict = dict()

        bclconvert_data: Dict = dict()
        for demux in demuxes.values():
            self.parse_demux_data(demux, bclconvert_data, len(demuxes))
        for qmetric in qmetrics.values():
            self.parse_qmetrics_data(bclconvert_data, qmetric)

        if len(demuxes) == 0:
            raise ModuleNoSamplesFound
        elif len(demuxes) > 1 and not multiple_sequencing_runs:
            log.warning("Found multiple runs from the same sequencer output")
            self.intro += """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Detected multiple bclconvert runs from the same sequencer output.
                    Runs were merged and undetermined stats were recalculated.
                </div>
            """
            self._recalculate_undetermined(bclconvert_data, last_run_id)
        elif multiple_sequencing_runs:
            # If we have data from multiple sequencing runs, the recalculation in
            # _recalculate_undetermined(last_run_id) wont work. In this case we
            # suppress/hide the info.
            log.warning("Found multiple sequencer runs")
            self.intro += """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Detected multiple sequencer runs.
                    Sample stats were merged.
                </div>
            """
            self.per_lane_undetermined_reads = None

        create_undetermined_barplots = (
            getattr(config, "bclconvert", {}).get("create_undetermined_barcode_barplots", False) or len(demuxes) == 1
        )
        if create_undetermined_barplots:
            self._parse_top_unknown_barcodes(bclconvert_data, last_run_id)

        # Collect counts by lane and sample
        (
            bclconvert_by_lane,
            bclconvert_by_sample,
            counts_by_sample_by_lane,
        ) = self._split_data_by_lane_and_sample(bclconvert_data)

        # Filter to strip out ignored sample names
        bclconvert_by_lane = self.ignore_samples(bclconvert_by_lane)
        bclconvert_by_sample = self.ignore_samples(bclconvert_by_sample)
        counts_by_sample_by_lane = self.ignore_samples(counts_by_sample_by_lane)

        # Return with Warning if no files are found
        if len(bclconvert_by_lane) == 0 and len(bclconvert_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"{len(bclconvert_by_lane)} lanes and {len(bclconvert_by_sample)} samples found")

        # Calculate mean quality scores
        for sample_id, sample in bclconvert_by_sample.items():
            if sample["yield"] > 0:
                sample["mean_quality"] = sample["_quality_score_sum"] / sample["yield"]
            del sample["_quality_score_sum"]

        for lane_id, lane in bclconvert_by_lane.items():
            if lane["yield"] > 0:
                lane["mean_quality"] = lane["_quality_score_sum"] / lane["yield"]
            del lane["_quality_score_sum"]

        self.write_data_file(bclconvert_by_lane, "multiqc_bclconvert_bylane")
        self.write_data_file(bclconvert_by_sample, "multiqc_bclconvert_bysample")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add sections for summary stats per flow cell
        self.add_section(
            name="Sample Statistics",
            anchor="bclconvert-samplestats",
            description="Statistics about each sample for each flowcell",
            plot=self.sample_stats_table(bclconvert_data, bclconvert_by_sample),
        )

        self.add_section(
            name="Lane Statistics",
            anchor="bclconvert-lanestats",
            description="Statistics about each lane for each flowcell",
            plot=self.lane_stats_table(bclconvert_by_lane),
        )

        # Add section for counts by lane
        cats = {
            "perfect": {"name": "Perfect Index Reads"},
            "imperfect": {"name": "Mismatched Index Reads"},
            "undetermined": {"name": "Undetermined Reads"},
        }
        extra = ""
        if len(demuxes) > 1 and not multiple_sequencing_runs:
            extra = """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Found multiple runs from the same sequencer output.
                    Runs were merged and <em>Undetermined Reads</em> were recalculated.
                </div>
            """
        elif multiple_sequencing_runs:
            extra = """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Found multiple sequencer runs.
                    <em>Undetermined Reads</em> cannot be recalculated for multiple sequencing runs and are not shown.
                </div>
            """
        self.add_section(
            name="Clusters by lane",
            anchor="bclconvert-bylane",
            description="Number of reads per lane (with number of perfect index reads)." + extra,
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.""",
            plot=bargraph.plot(
                self.get_bar_data_from_counts(bclconvert_data, bclconvert_by_lane, last_run_id),
                cats,
                {
                    "id": "bclconvert_lane_counts",
                    "title": "bclconvert: Clusters by lane",
                    "ylab": "Number of clusters",
                    "hide_empty": False,
                },
            ),
        )

        self.add_section(
            name="Clusters by sample",
            anchor="bclconvert-bysample",
            description="Number of reads per sample.",
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combined. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot=bargraph.plot(
                [
                    self.get_bar_data_from_counts(bclconvert_data, bclconvert_by_sample, last_run_id),
                    counts_by_sample_by_lane,
                ],
                [cats, sorted(bclconvert_by_lane.keys())],
                {
                    "id": "bclconvert_sample_counts",
                    "title": "bclconvert: Clusters by sample",
                    "hide_empty": False,
                    "ylab": "Number of clusters",
                    "data_labels": ["Index mismatches", "Counts per lane"],
                },
            ),
        )

        # Add section with undetermined barcodes
        if create_undetermined_barplots:
            undetermined_data = self.get_bar_data_from_undetermined(bclconvert_by_lane)
            if undetermined_data:
                self.add_section(
                    name="Undetermined barcodes by lane",
                    anchor="undetermine_by_lane",
                    description="Undetermined barcodes by lanes",
                    plot=bargraph.plot(
                        undetermined_data,
                        None,
                        {
                            "id": "bclconvert_undetermined",
                            "title": "bclconvert: Undetermined barcodes by lane",
                            "ylab": "Count",
                            "use_legend": True,
                            "tt_suffix": "reads",
                        },
                    ),
                )
            else:
                self.add_section(
                    name="Undetermined barcodes by lane",
                    anchor="undetermine_by_lane",
                    content="<div class='alert alert-info'>No undetermined barcodes found</div>",
                )

    @staticmethod
    @functools.lru_cache
    def _get_genome_size() -> Optional[int]:
        gs = getattr(config, "bclconvert", {}).get("genome_size")
        if gs:
            try:
                gs = int(float(gs))
            except ValueError:
                presets = {"hg19_genome": 2897310462, "hg38_genome": 3049315783, "mm10_genome": 2652783500}
                if gs in presets:
                    gs = presets[gs]
                else:
                    log.warning(
                        f"Genome '{gs}' is unknown to MultiQC, please specify size explicitly or choose one of the "
                        f"following: {', '.join(presets.keys())}."
                    )
                    gs = None
        log.debug(f"Determined genome size as: {gs}")
        return gs

    @staticmethod
    def _reads_dictionary():
        return {
            "clusters": 0,
            "yield": 0,
            "_calculated_yield": 0,  # re-calculated yield from demux stats where it is not provided explicitly and is calculated from # Reads and Read Length
            "perfect_index_reads": 0,
            "one_mismatch_index_reads": 0,
            "basesQ30": 0,
            "depth": 0,
            "_quality_score_sum": 0,  # used to re-calculate mean_quality
            "_calculated_quality_score_sum": 0,  # used to re-calculate mean_quality from demux stats where Yield is not provided explicitly and is calculated from # Reads and Read Length
        }

    @staticmethod
    def _is_single_end_reads(root):
        ncount = 0
        # we have single-end reads if we have exactly one IsIndexedRead=N element in our RunInfo.xml
        for element in root.findall("./Run/Reads/Read"):
            if element.get("IsIndexedRead") == "N":
                ncount += 1
        if ncount == 1:
            return True
        return False

    @staticmethod
    def _get_r2_length(root):
        for element in root.findall("./Run/Reads/Read"):
            if element.get("Number") != "1" and element.get("IsIndexedRead") == "N":
                return element.get("NumCycles")

        log.error("Could not figure out read 2 length from RunInfo.xml")
        raise ModuleNoSamplesFound

    def _parse_single_runinfo_file(self, runinfo_file):
        """
        Get run id and cluster length from RunInfo.xml
        """
        # Find all reads with IsIndexedRead = N
        root = ET.fromstring(runinfo_file["f"])
        run_id = root.find("Run").get("Id")
        reads = root.findall("./Run/Reads/Read")

        non_index_reads = [element for element in reads if element.get("IsIndexedRead") == "N"]
        if len(non_index_reads) == 0:
            log.error("No non-index reads found in RunInfo.xml")
            raise ModuleNoSamplesFound

        read_length_r1 = int(non_index_reads[0].get("NumCycles"))
        cluster_length = read_length_r1
        if len(non_index_reads) > 1:
            read_length_r2 = int(non_index_reads[1].get("NumCycles"))
            cluster_length += read_length_r2

        self.add_data_source(
            runinfo_file,
            module="bclconvert",
            section="bclconvert-runinfo-xml",
        )
        return {"run_id": run_id, "cluster_length": cluster_length}

    def _collate_log_files(self):
        # This function returns a list of self.find_log_files('bclconvert/demux') dicts,
        # with the run_id added on, sorted by root directory
        #
        # To get all that we must parse runinfo files and match them with demux files,
        # because demux files don't contain run-ids we need to match demux and runinfo
        # logs from the same directory, but find_log_files() does not guarantee order;
        # however it provides root dir, so we use that.
        demuxes_by_root = {f["root"]: f for f in self.find_log_files("bclconvert/demux")}
        runinfos_by_root = {f["root"]: f for f in self.find_log_files("bclconvert/runinfo")}
        qmetrics_by_root = {f["root"]: f for f in self.find_log_files("bclconvert/quality_metrics")}

        for root in runinfos_by_root.copy().keys():
            if root not in demuxes_by_root:
                log.error(f"Found RunInfo.xml file in {root} but no Demux Stats file, skipping")
                del runinfos_by_root[root]
                continue
        for root in demuxes_by_root.copy().keys():
            if root not in runinfos_by_root:
                log.error(f"Found Demux Stats file in {root} but no RunInfo.xml file, skipping")
                del demuxes_by_root[root]
                continue

        multiple_sequencing_runs = None
        last_run_id = None
        for root, demux in demuxes_by_root.items():
            runinfo = self._parse_single_runinfo_file(runinfos_by_root[root])
            demux["run_id"] = runinfo["run_id"]
            demux["cluster_length"] = runinfo["cluster_length"]
            qmetrics = qmetrics_by_root.get(root)
            if qmetrics:
                qmetrics["run_id"] = runinfo["run_id"]

            if last_run_id and runinfo["run_id"] != last_run_id:
                # this will mean we supress unknown reads, since we can't do a recalculation
                multiple_sequencing_runs = True
            last_run_id = runinfo["run_id"]

        return demuxes_by_root, qmetrics_by_root, multiple_sequencing_runs, last_run_id

    def _recalculate_undetermined(self, bclconvert_data, last_run_id):
        # We have to calculate "corrected" unknown read counts when parsing more than
        # one bclconvert run. To do this: add up all the reads in a lane that were
        # assigned to samples, then take the total reads in a lane (which is taken
        # _from the sum of all reads in a single file_), subtract the former from the
        # latter, and use that as "undetermined samples in lane."
        total_reads_per_lane = dict()
        for filename, lanedata in self.total_reads_in_lane_per_file.items():
            for lane_id, reads in lanedata.items():
                if lane_id not in total_reads_per_lane:
                    total_reads_per_lane[lane_id] = int()
                if total_reads_per_lane[lane_id] and reads != total_reads_per_lane[lane_id]:
                    log.error(
                        "Warning: different amounts of reads per lane across input files! "
                        "Cannot expect calculations to be accurate!"
                    )
                total_reads_per_lane[lane_id] = reads

        run_data = bclconvert_data[last_run_id]  # in this situation we have only one run id
        for lane_id, lane in run_data.items():
            determined_reads = 0
            for sample_id, sample in lane["samples"].items():
                determined_reads += sample["clusters"]
            self.per_lane_undetermined_reads[lane_id] = total_reads_per_lane[lane_id] - determined_reads

    def parse_demux_data(self, demux_file, bclconvert_data, num_demux_files):
        """
        Parse a bclconvert output stats csv, populate variables appropriately
        """
        filename = str(os.path.join(demux_file["root"], demux_file["fn"]))
        total_reads_in_lane = dict()

        run_data = bclconvert_data.get(demux_file["run_id"], dict())
        bclconvert_data[demux_file["run_id"]] = run_data
        with open(filename) as fh:
            reader: csv.DictReader = csv.DictReader(fh, delimiter=",")
            for row in reader:
                lane_id = f"L{row['Lane']}"
                lane = run_data.get(lane_id)
                if lane is None:
                    lane = dict(samples={}, cluster_length=demux_file["cluster_length"], **self._reads_dictionary())
                    run_data[lane_id] = lane
                    self.per_lane_undetermined_reads[lane_id] = 0

                sname = row["SampleID"]
                self.add_data_source(
                    demux_file,
                    s_name=sname,
                    module="bclconvert",
                    section="bclconvert-runinfo-demux-csv",
                )
                if sname != "Undetermined":
                    # Don't include undetermined reads at all in any of the calculations...
                    if sname not in lane["samples"]:
                        lane["samples"][sname] = self._reads_dictionary()
                        lane["samples"][sname]["filename"] = filename
                        lane["samples"][sname]["samples"] = {}

                    sample = lane["samples"][sname]  # this sample in this lane

                    # total lane stats
                    lane["clusters"] += int(row["# Reads"])
                    lane["_calculated_yield"] += int(row["# Reads"]) * demux_file["cluster_length"]
                    lane["perfect_index_reads"] += int(row["# Perfect Index Reads"])
                    lane["one_mismatch_index_reads"] += int(row["# One Mismatch Index Reads"])
                    lane["basesQ30"] += int(row.get("# of >= Q30 Bases (PF)", "0"))  # Column only present pre v3.9.3

                    # stats for this sample in this lane
                    sample["clusters"] += int(row["# Reads"])
                    sample["_calculated_yield"] += int(row["# Reads"]) * demux_file["cluster_length"]
                    sample["perfect_index_reads"] += int(row["# Perfect Index Reads"])
                    sample["one_mismatch_index_reads"] += int(row["# One Mismatch Index Reads"])
                    sample["index"] = str(row["Index"])
                    try:
                        # Not all demux files have Sample_Project column
                        sample["sample_project"] = str(row["Sample_Project"])
                    except KeyError:
                        pass

                    # columns only present pre v3.9.3, after they moved to quality_metrics
                    sample["basesQ30"] += int(row.get("# of >= Q30 Bases (PF)", 0))
                    # Collecting to re-calculate mean_quality:
                    calculated_quality_score_sum = (
                        float(row.get("Mean Quality Score (PF)", 0)) * sample["_calculated_yield"]
                    )
                    sample["_calculated_quality_score_sum"] += calculated_quality_score_sum
                    lane["_calculated_quality_score_sum"] += calculated_quality_score_sum

                if lane_id not in total_reads_in_lane:
                    total_reads_in_lane[lane_id] = 0

                # Add up number of reads, regardless of undetermined or not
                total_reads_in_lane[lane_id] += int(row["# Reads"])

                if num_demux_files == 1 and sname == "Undetermined":
                    self.per_lane_undetermined_reads[lane_id] += int(row["# Reads"])

        self.total_reads_in_lane_per_file[filename] = total_reads_in_lane

    def parse_qmetrics_data(self, bclconvert_data, qmetrics_file):
        """
        Parse a bclconvert output stats CSV, populate variables appropriately
        """
        filename = str(os.path.join(qmetrics_file["root"], qmetrics_file["fn"]))
        self.total_reads_in_lane_per_file[filename] = dict()

        reader: csv.DictReader = csv.DictReader(open(filename), delimiter=",")
        for row in reader:
            run_data = bclconvert_data[qmetrics_file["run_id"]]
            lane_id = f"L{row['Lane']}"
            if lane_id not in run_data:
                log.warning(f"Found unrecognised lane {lane_id} in Quality Metrics file, skipping")
                continue
            lane = run_data[lane_id]
            sample = row["SampleID"]
            self.add_data_source(
                qmetrics_file,
                s_name=sample,
                module="bclconvert",
                section="bclconvert-runinfo-quality-metrics-csv",
            )
            if sample != "Undetermined":  # don't include undetermined reads at all in any of the calculations...
                if sample not in run_data[lane_id]["samples"]:
                    log.warning(f"Found unrecognised sample {sample} in Quality Metrics file, skipping")
                    continue
                lane_sample = run_data[lane_id]["samples"][sample]  # this sample in this lane

                # Parse the stats that moved to this file in v3.9.3
                lane["yield"] += int(row["Yield"])
                lane["basesQ30"] += int(row["YieldQ30"])
                lane_sample["yield"] += int(row["Yield"])
                lane_sample["basesQ30"] += int(row["YieldQ30"])
                # Collecting to re-calculate mean_quality:
                lane["_quality_score_sum"] += float(row["QualityScoreSum"])
                lane_sample["_quality_score_sum"] += float(row["QualityScoreSum"])

    def _parse_top_unknown_barcodes(self, bclconvert_data, last_run_id):
        run_data = bclconvert_data[last_run_id]

        for unknown_barcode_file in self.find_log_files("bclconvert/unknown_barcodes", filehandles=True):
            barcode_reader = csv.DictReader(unknown_barcode_file["f"], delimiter=",")
            for unknown_barcode_row in barcode_reader:
                thislane = "L" + str(unknown_barcode_row["Lane"])
                if "top_unknown_barcodes" not in run_data[thislane]:
                    run_data[thislane]["top_unknown_barcodes"] = {}
                thisbarcode = str(unknown_barcode_row["index"]) + "-" + str(unknown_barcode_row["index2"])
                run_data[thislane]["top_unknown_barcodes"][thisbarcode] = int(unknown_barcode_row["# Reads"])

    @staticmethod
    def _total_reads_for_run(bclconvert_data, run_id):
        totalreads = 0
        for lane_id, lane in bclconvert_data[run_id].items():
            totalreads += lane["clusters"]
        return totalreads

    @staticmethod
    def _total_reads_all_runs(bclconvert_data):
        totalreads = 0
        for key, run_data in bclconvert_data.items():
            for lane_id, lane in run_data.items():
                totalreads += lane["clusters"]
        return totalreads

    def _set_lane_percentage_stats(self, data, cluster_length):
        try:
            data["percent_Q30"] = (float(data["basesQ30"]) / float(data["clusters"] * cluster_length)) * 100.0
        except ZeroDivisionError:
            data["percent_Q30"] = "NA"
        try:
            data["percent_perfectIndex"] = (float(data["perfect_index_reads"]) / float(data["clusters"])) * 100.0
        except ZeroDivisionError:
            data["percent_perfectIndex"] = "NA"
        try:
            data["percent_oneMismatch"] = float(data["one_mismatch_index_reads"]) / float(data["clusters"]) * 100.0
        except ZeroDivisionError:
            data["percent_oneMismatch"] = "NA"
        if self._get_genome_size():
            data["depth"] = float(data["basesQ30"]) / self._get_genome_size()
        else:
            data["depth"] = "NA"

    def _split_data_by_lane_and_sample(self, bclconvert_data) -> Tuple[Dict, Dict, Dict]:
        """
        Populate a collection of "stats across all lanes" and "stats across all samples"
        """
        bclconvert_by_lane = dict()
        bclconvert_by_sample = dict()
        count_by_sample_by_lane = defaultdict(lambda: defaultdict(int))  # counts only - used for a stacked bargraph
        for run_id, r in bclconvert_data.items():
            # set stats for each lane (across all samples) in bclconvert_bylane dictionary
            for lane_id, lane in r.items():
                self._set_lane_percentage_stats(lane, lane["cluster_length"])

                lane_key_name = self.prepend_runid(run_id, lane_id)
                bclconvert_by_lane[lane_key_name] = {
                    "depth": lane["depth"],
                    "clusters": lane["clusters"],
                    "yield": lane["yield"] or lane["_calculated_yield"],
                    "perfect_index_reads": lane["perfect_index_reads"],
                    "one_mismatch_index_reads": lane["one_mismatch_index_reads"],
                    "basesQ30": lane["basesQ30"],
                    "percent_Q30": lane["percent_Q30"],
                    "percent_perfectIndex": lane["percent_perfectIndex"],
                    "percent_oneMismatch": lane["percent_oneMismatch"],
                    "top_unknown_barcodes": lane["top_unknown_barcodes"] if "top_unknown_barcodes" in lane else {},
                    "_quality_score_sum": lane["_quality_score_sum"] or lane["_calculated_quality_score_sum"],
                }

                # now set stats for each sample (across all lanes) in bclconvert_bysample dictionary
                for sample_id, sample in lane["samples"].items():
                    if sample_id not in bclconvert_by_sample:
                        bclconvert_by_sample[sample_id] = self._reads_dictionary()

                    s = bclconvert_by_sample[sample_id]
                    s["clusters"] += int(sample["clusters"])
                    s["yield"] += sample["yield"] or sample["_calculated_yield"]
                    s["perfect_index_reads"] += sample["perfect_index_reads"]
                    s["one_mismatch_index_reads"] += sample["one_mismatch_index_reads"]
                    s["basesQ30"] += sample["basesQ30"]
                    s["cluster_length"] = lane["cluster_length"]
                    s["_quality_score_sum"] += sample["_quality_score_sum"] or sample["_calculated_quality_score_sum"]
                    s["index"] = sample["index"]
                    try:
                        # Not all demux files have Sample_Project column
                        s["sample_project"] = sample["sample_project"]
                    except KeyError:
                        pass

                    if not self._get_genome_size():
                        s["depth"] = "NA"
                    else:
                        s["depth"] += float(sample["basesQ30"]) / self._get_genome_size()

                    count_by_sample_by_lane[sample_id][lane_key_name] += sample["clusters"]

        return bclconvert_by_lane, bclconvert_by_sample, count_by_sample_by_lane

    def sample_stats_table(self, bclconvert_data, bclconvert_by_sample):
        sample_stats_data = dict()
        total_reads = self._total_reads_all_runs(bclconvert_data)
        depth_available = False

        for sample_id, sample in bclconvert_by_sample.items():
            # Percent stats for bclconvert-bysample i.e. stats for sample across all lanes
            try:
                perfect_percent = f"{float(100.0 * sample['perfect_index_reads'] / sample['clusters']):.1f}"
            except ZeroDivisionError:
                perfect_percent = "0.0"
            try:
                one_mismatch_percent = "{0:.1f}".format(
                    float(100.0 * sample["one_mismatch_index_reads"] / sample["clusters"])
                )
            except ZeroDivisionError:
                one_mismatch_percent = "0.0"

            try:
                yield_q30_percent = f"{float(100.0 * (sample['basesQ30'] / sample['yield'])):.1f}"
            except ZeroDivisionError:
                yield_q30_percent = "0.0"  #

            try:
                percent_yield = (float(sample["yield"]) / float(total_reads * sample["cluster_length"])) * 100.0
            except ZeroDivisionError:
                percent_yield = "NA"

            try:
                percent_reads = (float(sample["clusters"]) / float(total_reads)) * 100.0
            except ZeroDivisionError:
                percent_reads = "NA"

            sample_stats_data[sample_id] = {
                "depth": sample["depth"],
                "basesQ30": sample["basesQ30"],
                "clusters": sample["clusters"],
                "percent_reads": percent_reads,
                "yield": sample["yield"],
                "percent_yield": percent_yield,
                "yield_q30_percent": yield_q30_percent,
                # "perfect_index": samle['perfect_index_reads'], # don't need these
                # "one_mismatch_index_reads": sample['one_mismatch_index_reads'],
                "perfect_percent": perfect_percent,
                "one_mismatch_percent": one_mismatch_percent,
                "mean_quality": sample.get("mean_quality"),
                "index": sample["index"],
            }
            if sample["depth"] != "NA":
                depth_available = True
            # Not all demux files have Sample_Project column
            if "sample_project" in sample:
                sample_stats_data[sample_id]["sample_project"] = sample["sample_project"]

        headers = {}
        if depth_available:
            headers["depth"] = {
                "title": "Coverage",
                "description": (
                    "Estimated sequencing depth based on the number of bases with quality score greater or equal to Q30, "
                    + (
                        f", assuming the genome size is {self._get_genome_size()} as provided in config"
                        if self._get_genome_size() is not None
                        else ""
                    )
                ),
                "min": 0,
                "suffix": "X",
                "scale": "BuPu",
            }

        headers["clusters"] = {
            "title": f"{config.read_count_prefix} Clusters",
            "description": "Total number of clusters (read pairs) for this sample as determined by bclconvert demultiplexing ({})".format(
                config.read_count_desc
            ),
            "scale": "Blues",
            "shared_key": "read_count",
        }
        headers["yield"] = {
            "title": f"Yield ({config.base_count_prefix})",
            "description": "Total number of bases for this sample as determined by bclconvert demultiplexing ({})".format(
                config.base_count_desc
            ),
            "scale": "Greens",
            "shared_key": "base_count",
        }
        headers["percent_reads"] = {
            "title": "% Clusters",
            "description": "Percentage of clusters (read pairs) for this sample in this run, as determined by bclconvert demultiplexing",
            "scale": "Blues",
            "max": 100,
            "min": 0,
            "suffix": "%",
        }
        headers["percent_yield"] = {
            "title": "% Yield",
            "description": "Percentage of sequenced bases for this sample in this run",
            "scale": "Greens",
            "max": 100,
            "min": 0,
            "suffix": "%",
        }
        headers["basesQ30"] = {
            "title": f"Bases ({config.base_count_prefix}) ≥ Q30 (PF)",
            "description": "Number of bases with a Phred score of 30 or higher, passing filter ({})".format(
                config.base_count_desc
            ),
            "scale": "Blues",
            "shared_key": "base_count",
        }
        headers["yield_q30_percent"] = {
            "title": "% Bases ≥ Q30 (PF)",
            "description": "Percent of bases with a Phred score of 30 or higher, passing filter ({})".format(
                config.base_count_desc
            ),
            "scale": "Greens",
            "max": 100,
            "min": 0,
            "suffix": "%",
        }
        headers["perfect_percent"] = {
            "title": "% Perfect Index",
            "description": "Percent of reads with perfect index (0 mismatches)",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["one_mismatch_percent"] = {
            "title": "% One Mismatch Index",
            "description": "Percent of reads with one mismatch index",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["mean_quality"] = {
            "title": "Mean Quality Score",
            "description": "Mean quality score of bases",
            "min": 0,
            "max": 40,
            "scale": "RdYlGn",
        }
        headers["index"] = {
            "title": "Index",
            "description": "Sample index",
            "scale": False,
            "hidden": True,
        }
        headers["sample_project"] = {
            "title": "Project",
            "description": "Sample project",
            "scale": False,
            "hidden": True,
        }

        # Table config
        table_config = {
            "namespace": "bclconvert",
            "id": "bclconvert-sample-stats-table",
            "title": "bclconvert Sample Statistics",
        }

        return table.plot(sample_stats_data, headers, table_config)

    def lane_stats_table(self, bclconvert_by_lane):
        depth_available = False
        for lane_id, lane in bclconvert_by_lane.items():
            try:
                yield_q30_percent = f"{float(100.0 * (lane['basesQ30'] / lane['yield'])):.1f}"
            except ZeroDivisionError:
                yield_q30_percent = "0.0"
            bclconvert_by_lane[lane_id]["yield_q30_percent"] = yield_q30_percent
            if lane["depth"] != "NA":
                depth_available = True

        headers = {}
        if depth_available:
            headers["depth-lane"] = {
                "title": "Coverage",
                "description": (
                    "Estimated sequencing depth based on the number of bases with quality score greater or equal to Q30"
                    + (
                        f", assuming the genome size is {self._get_genome_size()} as provided in config"
                        if self._get_genome_size() is not None
                        else ""
                    )
                ),
                "suffix": "X",
                "scale": "BuPu",
            }

        headers["reads-lane"] = {
            "title": f"{config.read_count_prefix} Clusters",
            "description": "Total number of clusters (read pairs) for this sample as determined by bclconvert demultiplexing ({})".format(
                config.read_count_desc
            ),
            "scale": "Blues",
            "shared_key": "read_count",
        }
        headers["yield-lane"] = {
            "title": f"Yield ({config.base_count_prefix})",
            "description": "Total number of bases for this sample as determined by bclconvert demultiplexing ({})".format(
                config.base_count_desc
            ),
            "scale": "Greens",
            "shared_key": "base_count",
        }
        headers["basesQ30-lane"] = {
            "title": f"Bases ({config.base_count_prefix}) ≥ Q30 (PF)",
            "description": "Number of bases with a Phred score of 30 or higher, passing filter ({})".format(
                config.base_count_desc
            ),
            "scale": "Blues",
            "shared_key": "base_count",
        }
        headers["yield_q30_percent-lane"] = {
            "title": "% Bases ≥ Q30 (PF)",
            "description": "Percent of bases with a Phred score of 30 or higher, passing filter",
            "max": 100,
            "min": 0,
            "scale": "Greens",
        }
        headers["perfect_index_reads-lane"] = {
            "title": f"{config.read_count_prefix} Perfect Index",
            "description": f"Reads with perfect index - 0 mismatches ({config.read_count_desc})",
            "scale": "Blues",
            "shared_key": "read_count",
        }

        headers["one_mismatch_index_reads-lane"] = {
            "title": f"{config.read_count_prefix} One Mismatch",
            "description": f"Reads with one mismatch index ({config.read_count_desc})",
            "scale": "Spectral",
            "shared_key": "read_count",
        }
        headers["percent_perfectIndex-lane"] = {
            "title": "% Perfect Index",
            "description": "Percent of reads with perfect index - 0 mismatches",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["percent_oneMismatch-lane"] = {
            "title": "% One Mismatch",
            "description": "Percent of reads with one mismatch",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["mean_quality-lane"] = {
            "title": "Mean Quality Score",
            "description": "Mean quality score of bases",
            "min": 0,
            "max": 40,
            "scale": "PiYG",
        }

        # Table config
        table_config = {
            "namespace": "bclconvert-lane",
            "id": "bclconvert-lane-stats-table",
            "title": "bclconvert Lane Statistics",
            "col1_header": "Run ID - Lane",
        }

        # new dict with matching keys for plotting (this avoids duplicate html id linting errors)
        bclconvert_bylane_foroutput = dict()
        for laneid, lanestats in bclconvert_by_lane.items():
            if laneid not in bclconvert_bylane_foroutput:
                bclconvert_bylane_foroutput[laneid] = dict()
            for key, value in lanestats.items():
                bclconvert_bylane_foroutput[laneid][key + "-lane"] = value

        return table.plot(bclconvert_bylane_foroutput, headers, table_config)

    @staticmethod
    def prepend_runid(runid, rest):
        return str(runid) + " - " + str(rest)

    def get_bar_data_from_counts(self, bclconvert_data, counts, last_run_id) -> Dict[str, Dict[str, int]]:
        # For per-lane stats we fetch undetermined reads, too.
        bar_data = {}
        for key, value in counts.items():
            bar_data[key] = {
                "perfect": value["perfect_index_reads"],
                "imperfect": value["clusters"] - value["perfect_index_reads"],
            }
            try:
                if key.startswith(
                    self.prepend_runid(last_run_id, "")  # this wont run in multiple sequencing run situations
                ):  # per-lane stats start with a prepended run id, this is a per-lane entry
                    this_lane_id = key.replace(self.prepend_runid(last_run_id, ""), "")
                    rundata = bclconvert_data[last_run_id]
                    if this_lane_id in rundata:  # this is definitely a lane
                        bar_data[key]["undetermined"] = self.per_lane_undetermined_reads[this_lane_id]
            except TypeError:
                # do nothing, there is no Undetermined - this will happen in case of multiple run ids
                pass

        return bar_data

    @staticmethod
    def get_bar_data_from_undetermined(flowcells):
        """
        Get data to plot for undetermined barcodes.
        """

        bar_data = defaultdict(dict)
        # get undetermined barcodes for each lanes
        for lane_id, lane in flowcells.items():
            try:
                for barcode, count in islice(lane["top_unknown_barcodes"].items(), 20):
                    bar_data[barcode][lane_id] = count
            except AttributeError:
                pass
            except KeyError:
                pass

        return {key: value for key, value in islice(bar_data.items(), 20)}
