import csv
import functools
import logging
from collections import defaultdict
from itertools import islice
from pathlib import Path
from typing import List, Mapping, Set, Tuple, Dict, Optional, TypedDict, Union
from xml.etree import ElementTree

from pydantic import BaseModel

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.plots.plotly.violin import ViolinPlot
from multiqc.plots.table_object import ColumnDict, InputRow, ValueT
from multiqc.types import ColumnKey, LoadedFileDict, SampleGroup, SampleName

log = logging.getLogger(__name__)


class BaseMetrics(BaseModel):
    run_id: str
    cluster_length: int
    clusters: int = 0
    perfect_index_reads: int = 0
    one_mismatch_index_reads: int = 0
    percent_clusters: Optional[float] = None
    percent_perfect_index_reads: Optional[float] = None
    percent_one_mismatch_index_reads: Optional[float] = None
    yield_: Optional[int] = None
    percent_yield: Optional[float] = None
    yield_q30: Optional[int] = None
    percent_yield_q30: Optional[float] = None
    mean_quality: Optional[float] = None
    top_unknown_barcodes: Dict[str, int] = {}
    depth: Optional[float] = None
    # used to re-calculate mean_quality
    quality_score_sum: Optional[float] = None
    # re-calculated yield from demux stats where it is not provided explicitly and is calculated from # Reads
    # and Read Length:
    calculated_yield: int = 0
    # used to re-calculate mean_quality from demux stats where Yield is not provided explicitly and is
    # calculated from # Reads and Read Length
    calculated_qscore_sum: Optional[float] = None


class ChunkMetrics(BaseMetrics):
    """Data for one chunk (single run, single lane, single sample)"""

    index: Optional[str] = None
    sample_project: Optional[str] = None


class SampleSummary(BaseMetrics):
    """Data for a sample across all runs and lanes"""

    index: Optional[str] = None
    sample_project: Optional[str] = None
    lanes: Dict[str, ChunkMetrics] = {}


class LaneSummary(BaseMetrics):
    """All data that went through a lane on a run"""

    samples: Dict[str, ChunkMetrics] = {}


class RunSummary(BaseMetrics):
    """Summary for a run (all lanes)"""

    lanes: Dict[str, LaneSummary] = {}


class RunInfo(BaseModel):
    s_name: str
    path: Path
    cluster_length: int
    run_id: str


class BarDataDict(TypedDict, total=False):
    perfect: int
    imperfect: int
    undetermined: int


class MultiqcModule(BaseMultiqcModule):
    """
    This BclConvert module is based on the bcl2fastq multiqc module. It can parse multiple
    bclconvert run outputs as long as they are from the same sequencing run. When doing this,
    the undetermined reads will be 'corrected' and re-calculated (as an unknown read from
    one run might not be truly unknown, but simply from another run).

    #### Calculate estimated depth

    You can specify a genome size in config

    It's often useful to talk about sequencing yield in terms of estimated depth of coverage.
    In order to make MultiQC show the estimated depth for each sample, specify the reference genome/target size in
    your [MultiQC configuration](https://docs.seqera.io/multiqc/getting_started/config):

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
        super(MultiqcModule, self).__init__(
            name="BCL Convert",
            anchor="bclconvert",
            href="https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html",
            info="Demultiplexes data and converts BCL files to FASTQ file formats for downstream analysis.",
            # Can't find a DOI // doi=
        )

        # Set up and collate bclconvert run and demux files
        demuxes_files, qmetric_files = self._collate_log_files()
        if len(demuxes_files) == 0:
            raise ModuleNoSamplesFound

        # variables to store reads for undetermined read recalculation
        self.undetermined_reads_per_lane: Dict[str, int] = defaultdict(int)
        self.total_reads_in_lane_per_demuxfile: Dict[Path, Dict[str, int]] = dict()

        data_by_sample: Dict[str, SampleSummary] = dict()
        data_by_run: Dict[str, RunSummary] = dict()
        for demux in demuxes_files.values():
            self.parse_demux_data(
                demux,
                data_by_sample=data_by_sample,
                data_by_run=data_by_run,
                num_demux_files=len(demuxes_files),
            )
        # Return with Warning if no files are found
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found runs: {len(data_by_sample)}, samples: {len(data_by_sample)}")

        for qmetric in qmetric_files.values():
            self.parse_qmetrics_data(data_by_run, data_by_sample, qmetric)

        # Now that we collected all metrics from both demux and qmetrics files, we can calculate derivative metrics
        total_reads = self._total_reads_all_runs(data_by_run)
        for sample in data_by_sample.values():
            self._finalize_metrics(sample, total_reads)
        for run in data_by_run.values():
            self._finalize_metrics(run, total_reads)
            for lane in run.lanes.values():
                self._finalize_metrics(lane, total_reads)
                for chunk in lane.samples.values():
                    self._finalize_metrics(chunk, total_reads)

        if len(demuxes_files) > 1 and len(data_by_run) == 1:
            log.warning("Found multiple files for one sequencer run, recalculating undetermined reads")
            self._recalculate_undetermined(data_by_run)
            self.intro += """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Detected multiple bclconvert files from the same sequencer run.
                    Files were merged and undetermined stats were recalculated.
                </div>
            """
        elif len(data_by_run) > 1:
            # If we have data from multiple sequencing runs, the recalculation in
            # _recalculate_undetermined(last_run_id) won't work. In this case we
            # suppress/hide the info.
            log.warning("Found multiple sequencer runs, do not report undetermined stats")
            self.intro += """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Detected multiple sequencer runs.
                    Sample stats were merged.
                </div>
            """
            self.undetermined_reads_per_lane = {}

        if len(data_by_run) == 1:
            self._parse_top_unknown_barcodes(run=list(data_by_run.values())[0])

        self._write_data_files(data_by_sample, data_by_run)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add sections for summary stats per flow cell
        self.add_section(
            name="Sample Statistics",
            anchor="bclconvert-samplestats",
            description="Statistics about each sample for each flowcell",
            plot=self.sample_stats_table(data_by_sample),
        )

        self.add_section(
            name="Lane Statistics",
            anchor="bclconvert-lanestats",
            description="Statistics about each lane for each flowcell",
            plot=self.lane_stats_table(data_by_run),
        )

        extra = ""
        if len(demuxes_files) > 1 and not len(data_by_run) > 1:
            extra = """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Found multiple runs from the same sequencer output.
                    Runs were merged and <em>Undetermined Reads</em> were recalculated.
                </div>
            """
        elif len(data_by_run) > 1:
            extra = """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Found multiple sequencer runs.
                    <em>Undetermined Reads</em> cannot be recalculated for multiple sequencing runs and are not shown.
                </div>
            """

        self._clusters_by_lane_barplot(data_by_run, extra)

        self._clusters_by_sample_barplot(data_by_sample)

        TOP_N_UNDETERMINED_BARCODES = 40

        # Add section with undetermined barcodes
        if len(data_by_run) == 1:
            undetermined_data = self.get_bar_data_from_undetermined(data_by_run, top_n=TOP_N_UNDETERMINED_BARCODES)
            if undetermined_data:
                self.add_section(
                    name="Undetermined barcodes",
                    anchor="undetermine_by_lane",
                    description=f"Top {TOP_N_UNDETERMINED_BARCODES} undetermined barcodes",
                    plot=bargraph.plot(
                        undetermined_data,
                        None,
                        {
                            "id": "bclconvert_undetermined",
                            "title": f"bclconvert: top {TOP_N_UNDETERMINED_BARCODES} undetermined barcodes",
                            "ylab": "Count",
                            "use_legend": True,
                            "tt_suffix": "reads",
                            "sort_samples": False,
                        },
                    ),
                )
            else:
                self.add_section(
                    name="Undetermined barcodes by lane",
                    anchor="undetermine_by_lane",
                    content="<div class='alert alert-info'>No undetermined barcodes found</div>",
                )

    def _clusters_by_lane_barplot(self, data_by_run: Dict[str, RunSummary], extra: str):
        perferct_imperfect_split: Dict[str, Dict[str, int]] = {}
        for run_id, run in data_by_run.items():
            for lane_id, lane in run.lanes.items():
                runlane_id = f"{run_id} - {lane_id}" if len(data_by_run) > 1 else lane_id
                perferct_imperfect_split[runlane_id] = {
                    "perfect": lane.perfect_index_reads,
                    "imperfect": lane.clusters - lane.perfect_index_reads,
                }
                if len(data_by_run) == 1 and lane_id in self.undetermined_reads_per_lane:
                    assert lane_id == runlane_id
                    perferct_imperfect_split[runlane_id]["undetermined"] = self.undetermined_reads_per_lane[lane_id]

        self.add_section(
            name="Clusters by lane",
            anchor="bclconvert-bylane",
            description="Number of reads per lane (with number of perfect index reads)." + extra,
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.""",
            plot=bargraph.plot(
                data=perferct_imperfect_split,
                cats={
                    "perfect": {"name": "Perfect Index Reads"},
                    "imperfect": {"name": "Mismatched Index Reads"},
                    "undetermined": {"name": "Undetermined Reads"},
                },
                pconfig={
                    "id": "bclconvert_lane_counts",
                    "title": "bclconvert: Clusters by lane",
                    "ylab": "Number of clusters",
                    "hide_empty": False,
                },
            ),
        )

    def _clusters_by_sample_barplot(self, data_by_sample: Dict[str, SampleSummary]):
        perferct_imperfect_split: Dict[str, Dict[str, int]] = {}
        for s_name, sample in data_by_sample.items():
            perferct_imperfect_split[s_name] = {
                "perfect": sample.perfect_index_reads,
                "imperfect": sample.clusters - sample.perfect_index_reads,
            }

        lane_split: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        lane_ids: List[str] = []
        for sample_id, sample in data_by_sample.items():
            for lane_id, lane in sample.lanes.items():
                lane_split[sample_id][lane_id] += lane.clusters
                if lane_id not in lane_ids:
                    lane_ids.append(lane_id)

        self.add_section(
            name="Clusters by sample",
            anchor="bclconvert-bysample",
            description="Number of reads per sample.",
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combined. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot=bargraph.plot(
                data=[
                    perferct_imperfect_split,
                    lane_split,
                ],
                cats=[
                    {
                        "perfect": {"name": "Perfect Index Reads"},
                        "imperfect": {"name": "Mismatched Index Reads"},
                    },
                    sorted(lane_ids),
                ],
                pconfig={
                    "id": "bclconvert_sample_counts",
                    "title": "bclconvert: Clusters by sample",
                    "hide_empty": False,
                    "ylab": "Number of clusters",
                    "data_labels": ["Index mismatches", "Counts per lane"],
                },
            ),
        )

    def _write_data_files(self, data_by_sample: Dict[str, SampleSummary], data_by_run: Dict[str, RunSummary]):
        data_by_sample_flat = {sname: data.model_dump() for sname, data in data_by_sample.items()}

        self.write_data_file(data_by_sample_flat, "multiqc_bclconvert_bysample")

        data_by_lane = {}
        for run_id, run in data_by_run.items():
            for lane_id, lane in run.lanes.items():
                data_by_lane[f"{run_id} - {lane_id}"] = lane.model_dump()

        self.write_data_file(data_by_lane, "multiqc_bclconvert_bylane")

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
    def _is_single_end_reads(root: ElementTree.Element) -> bool:
        ncount = 0
        # we have single-end reads if we have exactly one IsIndexedRead=N element in our RunInfo.xml
        for element in root.findall("./Run/Reads/Read"):
            if element.get("IsIndexedRead") == "N":
                ncount += 1
        if ncount == 1:
            return True
        return False

    @staticmethod
    def _get_r2_length(root: ElementTree.Element) -> Optional[str]:
        for element in root.findall("./Run/Reads/Read"):
            if element.get("Number") != "1" and element.get("IsIndexedRead") == "N":
                return element.get("NumCycles")

        log.error("Could not figure out read 2 length from RunInfo.xml")
        raise ModuleNoSamplesFound

    def _parse_single_runinfo_file(self, runinfo_file: LoadedFileDict[str]) -> RunInfo:
        """
        Get run id and cluster length from RunInfo.xml
        """
        # Find all reads with IsIndexedRead = N
        root: ElementTree.Element = ElementTree.fromstring(runinfo_file["f"])
        run = root.find("Run")
        if run is None:
            log.error("No Run element found in RunInfo.xml")
            raise ModuleNoSamplesFound
        run_id = run.get("Id")
        if run_id is None:
            log.error("No Run ID found in RunInfo.xml")
            raise ModuleNoSamplesFound

        reads = root.findall("./Run/Reads/Read")
        non_index_reads = [element for element in reads if element.get("IsIndexedRead") == "N"]
        if len(non_index_reads) == 0:
            log.error("No non-index reads found in RunInfo.xml")
            raise ModuleNoSamplesFound

        read_length_r1 = non_index_reads[0].get("NumCycles")
        if read_length_r1 is None:
            log.error("No NumCycles field for read 1 found in RunInfo.xml")
            raise ModuleNoSamplesFound

        cluster_length = int(read_length_r1)
        if len(non_index_reads) > 1:
            read_length_r2 = non_index_reads[1].get("NumCycles")
            if read_length_r2 is None:
                log.error("No NumCycles field for read 2 found in RunInfo.xml")
                raise ModuleNoSamplesFound
            cluster_length += int(read_length_r2)

        self.add_data_source(
            runinfo_file,
            module="bclconvert",
            section="bclconvert-runinfo-xml",
        )
        return RunInfo(
            s_name=runinfo_file["s_name"],
            path=Path(runinfo_file["root"]) / runinfo_file["fn"],
            cluster_length=cluster_length,
            run_id=run_id,
        )

    def _collate_log_files(self) -> Tuple[Dict[str, RunInfo], Dict[str, RunInfo]]:
        # This function returns a list of self.find_log_files('bclconvert/demux') dicts,
        # with the run_id added on, sorted by root directory
        #
        # To get all that we must parse runinfo files and match them with demux files,
        # because demux files don't contain run-ids we need to match demux and runinfo
        # logs from the same directory, but find_log_files() does not guarantee order;
        # however it provides root dir, so we use that.
        _demuxes_by_root: Dict[str, LoadedFileDict[str]] = {
            f["root"]: f for f in self.find_log_files("bclconvert/demux")
        }
        _runinfos_by_root: Dict[str, LoadedFileDict[str]] = {
            f["root"]: f for f in self.find_log_files("bclconvert/runinfo")
        }
        _qmetrics_by_root: Dict[str, LoadedFileDict[str]] = {
            f["root"]: f for f in self.find_log_files("bclconvert/quality_metrics")
        }

        for root in _runinfos_by_root.copy().keys():
            if root not in _demuxes_by_root:
                log.error(f"Found RunInfo.xml file in {root} but no Demux Stats file, skipping")
                del _runinfos_by_root[root]
                continue
        for root in _demuxes_by_root.copy().keys():
            if root not in _runinfos_by_root:
                log.error(f"Found Demux Stats file in {root} but no RunInfo.xml file, skipping")
                del _demuxes_by_root[root]
                continue

        demuxes_by_root: Dict[str, RunInfo] = {}
        runinfos_by_root: Dict[str, RunInfo] = {}
        qmetrics_by_root: Dict[str, RunInfo] = {}

        for root, _demux in _demuxes_by_root.items():
            runinfo: RunInfo = self._parse_single_runinfo_file(_runinfos_by_root[root])
            runinfos_by_root[root] = runinfo
            demuxes_by_root[root] = RunInfo(
                s_name=_demux["s_name"],
                path=Path(_demux["root"]) / _demux["fn"],
                cluster_length=runinfo.cluster_length,
                run_id=runinfo.run_id,
            )
            if _qmetrics_f := _qmetrics_by_root.get(root):
                qmetrics_by_root[root] = RunInfo(
                    s_name=_qmetrics_f["s_name"],
                    path=Path(_qmetrics_f["root"]) / _qmetrics_f["fn"],
                    cluster_length=runinfo.cluster_length,
                    run_id=runinfo.run_id,
                )

        return demuxes_by_root, qmetrics_by_root

    def _finalize_metrics(self, metrics: BaseMetrics, total_reads: int):
        """
        Calculate derivative metrics: percentages, depth, mean quality.
        Set yield_ and quality_score_sum from calculated versions if unset.
        """

        if metrics.yield_ is None:
            metrics.yield_ = metrics.calculated_yield
        if metrics.quality_score_sum is None:
            metrics.quality_score_sum = metrics.calculated_qscore_sum

        if metrics.yield_q30 is not None and (metrics.clusters * metrics.cluster_length):
            metrics.percent_yield_q30 = (float(metrics.yield_q30) / (metrics.clusters * metrics.cluster_length)) * 100.0
        if metrics.clusters:
            metrics.percent_perfect_index_reads = (float(metrics.perfect_index_reads) / metrics.clusters) * 100.0
            metrics.percent_one_mismatch_index_reads = (
                float(metrics.one_mismatch_index_reads) / metrics.clusters * 100.0
            )
        if metrics.yield_q30 is not None and (gs := self._get_genome_size()) is not None:
            metrics.depth = float(metrics.yield_q30) / gs

        if metrics.clusters:
            metrics.percent_perfect_index_reads = float(metrics.perfect_index_reads) / metrics.clusters * 100.0
            metrics.percent_one_mismatch_index_reads = (
                float(metrics.one_mismatch_index_reads) / metrics.clusters * 100.0
            )

        if metrics.yield_:
            if metrics.yield_q30 is not None:
                metrics.percent_yield_q30 = float(metrics.yield_q30) / metrics.yield_ * 100.0

        if metrics.yield_ and total_reads * metrics.cluster_length:
            metrics.percent_yield = metrics.yield_ / (total_reads * metrics.cluster_length) * 100.0

        if total_reads:
            metrics.percent_clusters = float(metrics.clusters) / total_reads * 100.0

        if metrics.yield_ and metrics.quality_score_sum is not None:
            metrics.mean_quality = metrics.quality_score_sum / metrics.yield_
        metrics.quality_score_sum = None

    def _calculate_mean_quality(
        self,
        data_by_runlane: Dict[str, LaneSummary],
        data_by_sample: Dict[str, SampleSummary],
        data_by_run: Dict[str, RunSummary],
    ):
        for _, runlane in data_by_runlane.items():
            if runlane.yield_ and runlane.quality_score_sum is not None:
                runlane.mean_quality = runlane.quality_score_sum / runlane.yield_
            runlane.quality_score_sum = None
            for _, sample in runlane.samples.items():
                if sample.yield_ and sample.quality_score_sum is not None:
                    sample.mean_quality = sample.quality_score_sum / sample.yield_
                sample.quality_score_sum = None
        for _, sample in data_by_sample.items():
            if sample.yield_ and sample.quality_score_sum is not None:
                sample.mean_quality = sample.quality_score_sum / sample.yield_
            sample.quality_score_sum = None
        for _, run in data_by_run.items():
            if run.yield_ and run.quality_score_sum is not None:
                run.mean_quality = run.quality_score_sum / run.yield_
            run.quality_score_sum = None

    def _recalculate_undetermined(self, data_by_run: Dict[str, RunSummary]):
        # We have to calculate "corrected" unknown read counts when parsing more than
        # one demux file that belong to the same run. To do this: add up all the reads
        # in a lane that were assigned to samples, then take the total reads in a lane
        # (which is taken from the sum of all reads in a single file), subtract the former
        # from the latter, and use that as "undetermined samples in lane."
        total_reads_per_lane: Dict[str, int] = defaultdict(int)
        for _, cnt_by_lane in self.total_reads_in_lane_per_demuxfile.items():
            for lane_id, cnt in cnt_by_lane.items():
                if total_reads_per_lane[lane_id] != 0 and cnt != total_reads_per_lane[lane_id]:
                    log.error(
                        "Warning: different amounts of reads per lane across input files! "
                        "Cannot expect calculations to be accurate!"
                    )
                total_reads_per_lane[lane_id] = cnt

        for _, run in data_by_run.items():
            for lane_id, lane in run.lanes.items():
                determined_reads = 0
                for _, sample in lane.samples.items():
                    determined_reads += sample.clusters
                if self.undetermined_reads_per_lane:
                    self.undetermined_reads_per_lane[lane_id] = total_reads_per_lane[lane_id] - determined_reads

    def parse_demux_data(
        self,
        demux_file: RunInfo,
        data_by_run: Dict[str, RunSummary],
        data_by_sample: Dict[str, SampleSummary],
        num_demux_files: int,
    ):
        """
        Parse a bclconvert output stats csv, populate variables appropriately
        """
        total_reads_in_lane: Dict[str, int] = defaultdict(int)

        run_id = demux_file.run_id
        with demux_file.path.open() as fh:
            reader: csv.DictReader[str] = csv.DictReader(fh, delimiter=",")
            for row in reader:
                sname = row["SampleID"]
                if self.is_ignore_sample(sname):
                    continue
                lane_id = f"L{row['Lane']}"

                # Add up number of reads, regardless of undetermined or not
                total_reads_in_lane[lane_id] += int(row["# Reads"])

                # Adn don't include undetermined reads at all in any of the further calculations
                if sname == "Undetermined":
                    if num_demux_files == 1:
                        self.undetermined_reads_per_lane[lane_id] += int(row["# Reads"])
                    continue

                run = data_by_run.get(run_id)
                if run is None:
                    run = RunSummary(cluster_length=demux_file.cluster_length, run_id=run_id)
                    data_by_run[run_id] = run

                lane = run.lanes.get(lane_id)
                if lane is None:
                    lane = LaneSummary(cluster_length=demux_file.cluster_length, run_id=run_id)
                    run.lanes[lane_id] = lane

                sample = data_by_sample.get(sname)
                if sample is None:
                    sample = SampleSummary(cluster_length=demux_file.cluster_length, run_id=run_id)
                    data_by_sample[sname] = sample

                self.add_data_source(
                    path=demux_file.path,
                    s_name=sname,
                    module="bclconvert",
                    section="bclconvert-runinfo-demux-csv",
                )

                chunk = ChunkMetrics(
                    run_id=run_id,
                    cluster_length=demux_file.cluster_length,
                    clusters=int(row["# Reads"]),
                    calculated_yield=int(row["# Reads"]) * demux_file.cluster_length,
                    perfect_index_reads=int(row["# Perfect Index Reads"]),
                    one_mismatch_index_reads=int(row["# One Mismatch Index Reads"]),
                    index=str(row["Index"]),
                )
                lane.samples[sname] = chunk
                sample.lanes[lane_id] = chunk

                # Columns only present pre v3.9.3, after they moved to quality_metrics
                if (yield_q30 := row.get("# of >= Q30 Bases (PF)")) is not None:
                    chunk.yield_q30 = int(yield_q30)
                    run.yield_q30 = (run.yield_q30 or 0) + chunk.yield_q30
                    lane.yield_q30 = (lane.yield_q30 or 0) + chunk.yield_q30
                    sample.yield_q30 = (sample.yield_q30 or 0) + chunk.yield_q30
                if (qscore := row.get("Mean Quality Score (PF)")) is not None:
                    calc_qscore_sum = float(qscore) * chunk.calculated_yield
                    chunk.calculated_qscore_sum = calc_qscore_sum
                    run.calculated_qscore_sum = (run.calculated_qscore_sum or 0) + calc_qscore_sum
                    lane.calculated_qscore_sum = (lane.calculated_qscore_sum or 0) + calc_qscore_sum
                    sample.calculated_qscore_sum = (sample.calculated_qscore_sum or 0) + calc_qscore_sum

                # Not all demux files have Sample_Project column
                if (sproj := row.get("Sample_Project")) is not None:
                    chunk.sample_project = str(sproj)

                # Total run stats
                run.clusters += chunk.clusters
                run.calculated_yield += chunk.calculated_yield
                run.perfect_index_reads += chunk.perfect_index_reads
                run.one_mismatch_index_reads += chunk.one_mismatch_index_reads

                # Total lane stats
                lane.clusters += chunk.clusters
                lane.calculated_yield += chunk.calculated_yield
                lane.perfect_index_reads += chunk.perfect_index_reads
                lane.one_mismatch_index_reads += chunk.one_mismatch_index_reads

                # Total sample stats
                sample.clusters += chunk.clusters
                sample.calculated_yield += chunk.calculated_yield
                sample.perfect_index_reads += chunk.perfect_index_reads
                sample.one_mismatch_index_reads += chunk.one_mismatch_index_reads
                if chunk.sample_project != sample.sample_project:
                    log.warning(
                        f"Sample {sname} has different project names on different lanes: "
                        f"{chunk.sample_project} != {sample.sample_project}, overriding"
                    )
                    sample.sample_project = chunk.sample_project
                if sample.index is not None and chunk.index != sample.index:
                    log.warning(
                        f"Sample {sname} has different indices on different lanes: "
                        f"{chunk.index} != {sample.index}, overriding"
                    )
                    sample.index = chunk.index

        self.total_reads_in_lane_per_demuxfile[demux_file.path] = total_reads_in_lane

    def parse_qmetrics_data(
        self,
        data_by_run: Dict[str, RunSummary],
        data_by_sample: Dict[str, SampleSummary],
        qmetrics_file: RunInfo,
    ):
        """
        Parse a bclconvert output stats CSV, populate variables appropriately
        """
        self.total_reads_in_lane_per_demuxfile[qmetrics_file.path] = dict()

        reader: csv.DictReader[str] = csv.DictReader(qmetrics_file.path.open(), delimiter=",")
        for row in reader:
            s_name = row["SampleID"]
            if s_name == "Undetermined":  # don't include undetermined reads at all in any of the calculations
                continue
            if self.is_ignore_sample(s_name):
                continue

            # data_by_lane: Dict[str, RunLaneSummary] = data_by_lane_by_run[qmetrics_file.run_id]
            # data_by_sample: Dict[str, SampleSummary] = data_by_sample_by_run[qmetrics_file.run_id]
            lane_id = f"L{row['Lane']}"

            if qmetrics_file.run_id not in data_by_run:
                log.warning(f"Found unrecognised run {qmetrics_file.run_id} in Quality Metrics file, skipping")
                continue
            run = data_by_run[qmetrics_file.run_id]

            if lane_id not in run.lanes:
                log.warning(
                    f"Found unrecognised lane {lane_id} in Quality Metrics file for run {qmetrics_file.run_id}, skipping"
                )
                continue
            lane = run.lanes[lane_id]

            if s_name not in lane.samples or s_name not in data_by_sample:
                log.warning(f"Found unrecognised sample {s_name} in Quality Metrics file, skipping")
                continue
            sample = data_by_sample[s_name]

            self.add_data_source(
                path=qmetrics_file.path,
                s_name=s_name,
                module="bclconvert",
                section="bclconvert-runinfo-quality-metrics-csv",
            )

            # Parse the stats that moved to this file in v3.9.3
            run.yield_ = (run.yield_ or 0) + int(row["Yield"])
            run.yield_q30 = (run.yield_q30 or 0) + int(row["YieldQ30"])

            lane.yield_ = (lane.yield_ or 0) + int(row["Yield"])
            lane.yield_q30 = (lane.yield_q30 or 0) + int(row["YieldQ30"])

            sample.yield_ = (sample.yield_ or 0) + int(row["Yield"])
            sample.yield_q30 = (sample.yield_q30 or 0) + int(row["YieldQ30"])

            chunk = lane.samples[s_name]  # this sample in this lane
            chunk.yield_ = (chunk.yield_ or 0) + int(row["Yield"])
            chunk.yield_q30 = (chunk.yield_q30 or 0) + int(row["YieldQ30"])

            # Collecting to re-calculate mean_quality:
            chunk.quality_score_sum = (chunk.quality_score_sum or 0) + float(row["QualityScoreSum"])
            lane.quality_score_sum = (lane.quality_score_sum or 0) + float(row["QualityScoreSum"])
            sample.quality_score_sum = (sample.quality_score_sum or 0) + float(row["QualityScoreSum"])
            run.quality_score_sum = (run.quality_score_sum or 0) + float(row["QualityScoreSum"])

    def _parse_top_unknown_barcodes(self, run: RunSummary):
        for unknown_barcode_file in self.find_log_files("bclconvert/unknown_barcodes", filehandles=True):
            barcode_reader = csv.DictReader(unknown_barcode_file["f"], delimiter=",")
            for unknown_barcode_row in barcode_reader:
                lane_id = "L" + str(unknown_barcode_row["Lane"])
                if lane_id not in run.lanes:
                    log.warning(f"Found unrecognised lane {lane_id} in Top Unknown Barcode file, skipping")
                    continue
                barcode = str(unknown_barcode_row["index"]) + "-" + str(unknown_barcode_row["index2"])
                run.lanes[lane_id].top_unknown_barcodes[barcode] = int(unknown_barcode_row["# Reads"])

    # @staticmethod
    # def _total_reads_for_run(data_by_lane_by_run: Dict[str, Dict[str, RunLaneSummary]], run_id: str) -> int:
    #     totalreads = 0
    #     for _, lane in data_by_lane_by_run[run_id].items():
    #         totalreads += lane.clusters
    #     return totalreads

    @staticmethod
    def _total_reads_all_runs(data_by_run: Dict[str, RunSummary]) -> int:
        totalreads = 0
        for _, run in data_by_run.items():
            totalreads += run.clusters
        return totalreads

    def _set_lane_percentage_stats(self, data: LaneSummary, cluster_length: int):
        if data.yield_q30 is not None and (data.clusters * cluster_length):
            data.percent_yield_q30 = (float(data.yield_q30) / (data.clusters * cluster_length)) * 100.0
        if data.clusters:
            data.percent_perfect_index_reads = (float(data.perfect_index_reads) / data.clusters) * 100.0
            data.percent_one_mismatch_index_reads = float(data.one_mismatch_index_reads) / data.clusters * 100.0
        if data.yield_q30 is not None and (gs := self._get_genome_size()) is not None:
            data.depth = float(data.yield_q30) / gs

    def sample_stats_table(self, data_by_sample: Dict[str, SampleSummary]) -> ViolinPlot:
        depth_available = any(sample.depth is not None for sample in data_by_sample.values())

        rows_by_sample: Dict[SampleGroup, List[InputRow]] = {}
        for sname, sample in data_by_sample.items():
            rows = [
                InputRow(
                    sample=SampleName(sname),
                    data={ColumnKey(k.strip("_")): v for k, v in sample.__dict__.items()},
                )
            ]
            if len(sample.lanes) > 1:
                for lane_id, lane in sample.lanes.items():
                    rows.append(
                        InputRow(
                            sample=SampleName(sname + " - " + lane_id),
                            data={ColumnKey(k.strip("_")): v for k, v in lane.__dict__.items()},
                        )
                    )
            else:
                rows[0].sample = SampleName(sname + " (" + list(sample.lanes.keys())[0] + ")")
            rows_by_sample[SampleGroup(sname)] = rows

        headers: Dict[str, ColumnDict] = {}
        if depth_available:
            headers["depth"] = {
                "title": "Coverage",
                "description": (
                    "Estimated sequencing depth based on the number of bases with quality score greater or equal to "
                    "Q30, "
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
            "title": "Clusters",
            "description": f"Total number of clusters (read pairs) for this sample as determined by bclconvert demultiplexing ({config.read_count_desc})",
            "scale": "Blues",
            "shared_key": "read_count",
        }
        headers["yield"] = {
            "title": "Yield",
            "description": f"Total number of bases for this sample as determined by bclconvert demultiplexing ({config.base_count_desc})",
            "scale": "Greens",
            "shared_key": "base_count",
        }
        headers["percent_clusters"] = {
            "title": "Clusters/run",
            "description": "Percentage of clusters (read pairs) for this sample in this run, as determined by bclconvert demultiplexing",
            "scale": "Blues",
            "max": 100,
            "min": 0,
            "suffix": "%",
        }
        headers["percent_yield"] = {
            "title": "Yield/run",
            "description": "Percentage of sequenced bases for this sample in this run",
            "scale": "Greens",
            "max": 100,
            "min": 0,
            "suffix": "%",
        }
        headers["yield_q30"] = {
            "title": "Bases ≥ Q30 (PF)",
            "description": f"Number of bases with a Phred score of 30 or higher, passing filter ({config.base_count_desc})",
            "scale": "Blues",
            "shared_key": "base_count",
        }
        headers["percent_yield_q30"] = {
            "title": "Bases ≥ Q30 (PF)",
            "description": f"Percent of bases with a Phred score of 30 or higher, passing filter ({config.base_count_desc})",
            "scale": "Greens",
            "max": 100,
            "min": 0,
            "suffix": "%",
        }
        headers["percent_perfect_index_reads"] = {
            "title": "Perfect index",
            "description": "Percent of reads with perfect index (0 mismatches)",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["percent_one_mismatch_index_reads"] = {
            "title": "One mismatch index",
            "description": "Percent of reads with one mismatch index",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["mean_quality"] = {
            "title": "Mean quality sscore",
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

        # GroupT = Union[Mapping[ColumnKeyT, Optional[ValueT]], InputRow, Sequence[InputRow]]
        # SectionT = Mapping[GroupKeyT, GroupT]
        return table.plot(
            rows_by_sample,
            headers,
            table_config,
        )

    def lane_stats_table(self, data_by_run: Dict[str, RunSummary]) -> ViolinPlot:
        depth_available = any(run.depth is not None for run in data_by_run.values())

        rows_by_sample: Dict[SampleGroup, List[InputRow]] = {}
        for run_id, run in data_by_run.items():
            for lane_id, lane in run.lanes.items():
                runlane_id = lane_id if len(data_by_run) == 1 else run_id + " - " + lane_id
                rows = [
                    InputRow(
                        sample=SampleName(runlane_id),
                        data={ColumnKey(k.strip("_")): v for k, v in lane.__dict__.items()},
                    )
                ]
                if len(lane.samples) > 1:
                    for sname, sample in lane.samples.items():
                        rows.append(
                            InputRow(
                                sample=SampleName(runlane_id + " - " + sname),
                                data={ColumnKey(k.strip("_")): v for k, v in sample.__dict__.items()},
                            )
                        )
                else:
                    rows[0].sample = SampleName(runlane_id + " (" + list(lane.samples.keys())[0] + ")")
                rows_by_sample[SampleGroup(runlane_id)] = rows

        headers: Dict[str, ColumnDict] = {}
        if depth_available:
            headers["depth"] = {
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

        headers["clusters"] = {
            "title": "Clusters",
            "description": f"Total number of clusters (read pairs) for this sample as determined by bclconvert "
            f"demultiplexing ({config.read_count_desc})",
            "scale": "Blues",
            "shared_key": "read_count",
        }
        headers["yield"] = {
            "title": "Yield",
            "description": f"Total number of bases for this sample as determined by bclconvert demultiplexing ("
            f"{config.base_count_desc})",
            "scale": "Greens",
            "shared_key": "base_count",
        }
        headers["bases_q30"] = {
            "title": "Bases ≥ Q30 (PF)",
            "description": f"Number of bases with a Phred score of 30 or higher, p"
            f"assing filter ({config.base_count_desc})",
            "scale": "Blues",
            "shared_key": "base_count",
        }
        headers["percent_yield_q30"] = {
            "title": "Bases ≥ Q30 (PF)",
            "description": "Percent of bases with a Phred score of 30 or higher, passing filter",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "Greens",
        }
        headers["perfect_index_reads"] = {
            "title": "Perfect index",
            "description": f"Reads with perfect index - 0 mismatches ({config.read_count_desc})",
            "scale": "Blues",
            "shared_key": "read_count",
        }

        headers["one_mismatch_index_reads"] = {
            "title": "One mismatch",
            "description": f"Reads with one mismatch index ({config.read_count_desc})",
            "scale": "Spectral",
            "shared_key": "read_count",
        }
        headers["percent_perfect_index_reads"] = {
            "title": "Perfect index",
            "description": "Percent of reads with perfect index - 0 mismatches",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["percent_one_mismatch_index_reads"] = {
            "title": "One mismatch",
            "description": "Percent of reads with one mismatch",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["mean_quality"] = {
            "title": "Mean quality score",
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

        return table.plot(rows_by_sample, headers, table_config)

    @staticmethod
    def prepend_runid(runid: str, rest: str) -> str:
        return str(runid) + " - " + str(rest)

    def get_bar_data_from_counts(
        self,
        data_by_lane: Mapping[str, LaneSummary],
        total_runs: Optional[int] = None,
    ) -> Dict[str, Dict[str, int]]:
        # For per-lane stats we fetch undetermined reads, too.
        bar_data: Dict[str, Dict[str, int]] = {}
        for key, metrics in data_by_lane.items():
            if total_runs == 1:
                key = key.split(" - ")[1]

            bar_data[key] = {
                "perfect": metrics.perfect_index_reads,
                "imperfect": metrics.clusters - metrics.perfect_index_reads,
            }
            if total_runs == 1 and key in self.undetermined_reads_per_lane:
                bar_data[key]["undetermined"] = self.undetermined_reads_per_lane[key]

        return bar_data

    @staticmethod
    def get_bar_data_from_undetermined(
        data_by_run: Dict[str, RunSummary], top_n: int = 20
    ) -> Dict[str, Dict[str, int]]:
        """
        Get data to plot for undetermined barcodes.
        """

        bar_data: Dict[str, Dict[str, int]] = defaultdict(dict)
        # get undetermined barcodes for each lanes
        for _, run in data_by_run.items():
            for lane_id, lane in run.lanes.items():
                try:
                    for barcode, count in islice(lane.top_unknown_barcodes.items(), top_n):
                        bar_data[barcode][lane_id] = count
                except AttributeError:
                    pass
                except KeyError:
                    pass

        # Sort by value
        bar_data = dict(sorted(bar_data.items(), key=lambda item: sum(item[1].values()), reverse=True))
        return {key: value for key, value in islice(bar_data.items(), top_n)}
