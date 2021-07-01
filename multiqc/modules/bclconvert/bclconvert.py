import json
import logging
import operator
import os
from collections import OrderedDict, defaultdict
from itertools import islice
import csv
import decimal
from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule
import xml.etree.ElementTree as ET

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="bclconvert",
            anchor="bclconvert",
            href="https://support.illumina.com/",
            info="can be used to both demultiplex data and convert BCL files"
            " to FASTQ file formats for downstream analysis.",
        )

        # set up and collate bclconvert run and demux files
        self.last_run_id = None
        self.cluster_length = None
        self.multiple_sequencing_runs = False
        bclconvert_demuxes = self._collate_log_files()
        self.num_demux_files = len(bclconvert_demuxes)

        self.bclconvert_data = dict()
        for demux in bclconvert_demuxes:
            self.bclconvert_data[demux["run_id"]] = dict()

        # variables to store reads for undetermined read recalcuation
        self.per_lane_undetermined_reads = dict()
        self.total_reads_in_lane_per_file = dict()

        for demux in bclconvert_demuxes:
            self.parse_demux_data(demux)

        if self.num_demux_files == 0:
            raise UserWarning
        elif self.num_demux_files > 1 and not self.multiple_sequencing_runs:
            log.warning("Found multiple runs from the same sequencer output")
            self.intro += """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Detected multiple bclconvert runs from the same sequencer output.
                    Runs were merged and undetermined stats were recalculated.
                </div>
            """
            self._recalculate_undetermined()
        elif self.multiple_sequencing_runs:
            # if we have data from multple sequencing runs, the recalculation in _recalculate_undetermined() wont work. in this case we suppress/hide the info.
            log.warning("Found multiple sequencer runs")
            self.intro += """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Detected multiple sequencer runs.
                    Sample stats were merged.
                </div>
            """
            self.per_lane_undetermined_reads = None
        elif self.num_demux_files == 1:
            # only possible to calculate barcodes per lane when parsing a single bclconvert run
            self._parse_top_unknown_barcodes()

        # Collect counts by lane and sample
        self.bclconvert_bylane = dict()
        self.bclconvert_bysample = dict()
        self.bclconvert_bysample_lane = dict()
        self.source_files = dict()
        self._split_data_by_lane_and_sample()

        # Filter to strip out ignored sample names
        self.bclconvert_bylane = self.ignore_samples(self.bclconvert_bylane)
        self.bclconvert_bysample = self.ignore_samples(self.bclconvert_bysample)
        self.bclconvert_bysample_lane = self.ignore_samples(self.bclconvert_bysample_lane)

        # Return with Warning if no files are found
        if len(self.bclconvert_bylane) == 0 and len(self.bclconvert_bysample) == 0:
            raise UserWarning
        log.info("{} lanes and {} samples found".format(len(self.bclconvert_bylane), len(self.bclconvert_bysample)))

        # Print source files
        for s in self.source_files.keys():
            self.add_data_source(
                s_name=s,
                source=",".join(list(set(self.source_files[s]))),
                module="bclconvert",
                section="bclconvert-bysample",
            )

        self.write_data_file(
            {str(k): self.bclconvert_bylane[k] for k in self.bclconvert_bylane.keys()}, "multiqc_bclconvert_bylane"
        )
        self.write_data_file(self.bclconvert_bysample, "multiqc_bclconvert_bysample")

        # Add sections for summary stats per flow cell
        self.add_section(
            name="Sample Statistics",
            anchor="bclconvert-samplestats",
            description="Statistics about each sample for each flowcell",
            plot=self.sample_stats_table(),
        )

        self.add_section(
            name="Lane Statistics",
            anchor="bclconvert-lanestats",
            description="Statistics about each lane for each flowcell",
            plot=self.lane_stats_table(),
        )

        # Add section for counts by lane
        cats = OrderedDict()
        cats["perfect"] = {"name": "Perfect Index Reads"}
        cats["imperfect"] = {"name": "Mismatched Index Reads"}
        cats["undetermined"] = {"name": "Undetermined Reads"}
        extra = ""
        if self.num_demux_files > 1 and not self.multiple_sequencing_runs:
            extra = """
                <div class="alert alert-warning">
                    <strong>Warning:</strong> Found multiple runs from the same sequencer output.
                    Runs were merged and <em>Undetermined Reads</em> were recalculated.
                </div>
            """
        elif self.multiple_sequencing_runs:
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
                self.get_bar_data_from_counts(self.bclconvert_bylane),
                cats,
                {
                    "id": "bclconvert_lane_counts",
                    "title": "bclconvert: Clusters by lane",
                    "ylab": "Number of clusters",
                    "hide_zero_cats": False,
                },
            ),
        )

        # Add section for counts by sample
        # get cats for per-lane tab
        lcats = set()
        for s_name in self.bclconvert_bysample_lane:
            lcats.update(self.bclconvert_bysample_lane[s_name].keys())
        lcats = sorted(list(lcats))
        self.add_section(
            name="Clusters by sample",
            anchor="bclconvert-bysample",
            description="Number of reads per sample.",
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combinned. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot=bargraph.plot(
                [self.get_bar_data_from_counts(self.bclconvert_bysample), self.bclconvert_bysample_lane],
                [cats, lcats],
                {
                    "id": "bclconvert_sample_counts",
                    "title": "bclconvert: Clusters by sample",
                    "hide_zero_cats": False,
                    "ylab": "Number of clusters",
                    "data_labels": ["Index mismatches", "Counts per lane"],
                },
            ),
        )

        # Add section with undetermined barcodes
        if self.num_demux_files == 1:
            self.add_section(
                name="Undetermined barcodes by lane",
                anchor="undetermine_by_lane",
                description="Undetermined barcodes by lanes",
                plot=bargraph.plot(
                    self.get_bar_data_from_undetermined(self.bclconvert_bylane),
                    None,
                    {
                        "id": "bclconvert_undetermined",
                        "title": "bclconvert: Undetermined barcodes by lane",
                        "ylab": "Count",
                        "tt_percentages": False,
                        "use_legend": True,
                        "tt_suffix": "reads",
                    },
                ),
            )

    def _get_genome_size(self):
        gs = getattr(config, "bclconvert", {}).get("genome_size")
        if gs:
            try:
                gs = float(gs)
            except ValueError:
                presets = {"hg19_genome": 2897310462, "hg38_genome": 3049315783, "mm10_genome": 2652783500}
                if gs in presets:
                    gs = presets[gs]
                else:
                    log.warning(
                        "The size for genome "
                        + gs
                        + " is unknown to MultiQC, "
                        + "please specify it explicitly or choose one of the following: "
                        + ", ".join(presets.keys())
                        + "."
                    )
                    gs = None
        log.debug("Determied genome size as " + str(gs))
        return gs

    def _set_up_reads_dictionary(self, dict):
        dict["reads"] = 0
        dict["yield"] = 0
        dict["perfect_index_reads"] = 0
        dict["one_mismatch_index_reads"] = 0
        dict["basesQ30"] = 0
        dict["mean_quality"] = 0

    def _parse_single_runinfo_file(self, runinfo_file):
        # get run id and readlength from RunInfo.xml
        try:
            root = ET.fromstring(runinfo_file["f"])
            run_id = root.find("Run").get("Id")
            read_length = root.find("./Run/Reads/Read[1]").get(
                "NumCycles"
            )  # ET indexes first element at 1, so here we're gettign the first NumCycles
        except:
            log.error(f"Could not parse RunInfo.xml to get RunID and read length in '{runinfo_file['root']}'")
            raise UserWarning
        return {"run_id": run_id, "read_length": int(read_length), "cluster_length": int(read_length) * 2}

    def _find_log_files_and_sort(self, search_pattern, sort_field):
        logs = list()
        for m in self.find_log_files(search_pattern):
            logs.append(m)
        return sorted(logs, key=lambda i: i[sort_field])  # sort on root directory

    def _collate_log_files(self):
        # this function returns a list of self.find_log_files('bclconvert/demux') dicts, with the run_id added on, sorted by root directory
        #
        # to get all that we must parse runinfo files and match them with demux files, because demux files dont contain run-ids
        # we need to match demux and runinfo logs from the same directory, but find_log_files() does not guarantee order; however it provides root dir, so we use that

        demuxes = self._find_log_files_and_sort("bclconvert/demux", "root")
        runinfos = self._find_log_files_and_sort("bclconvert/runinfo", "root")

        if not len(demuxes) == len(runinfos):
            log.error(
                f"Different amount of Demux Stats files and RunInfoXML files found: '{len(demuxes)}' and '{len(runinfos)}'"
            )
            raise UserWarning

        for idx, runinfo in enumerate(runinfos):
            rundata = self._parse_single_runinfo_file(runinfo)

            if runinfo["root"] != demuxes[idx]["root"]:
                log.error(
                    f"Expected RunInfo.xml file in '{runinfo['root']}' and Demultiplex_Stats.csv in '{demuxes[idx]['root']}' to be in the same directory"
                )
                raise UserWarning
            else:
                demuxes[idx]["run_id"] = rundata["run_id"]

            if self.last_run_id and rundata["run_id"] != self.last_run_id:
                self.multiple_sequencing_runs = (
                    True  # this will mean we supress unknown reads, since we can't do a recaculation
                )

            self.last_run_id = rundata["run_id"]

            if self.cluster_length and rundata["cluster_length"] != self.cluster_length:
                log.error(
                    "Detected different read lengths across the RunXml files. This should not be - read length across all input directories must be the same."
                )
                raise UserWarning
            self.cluster_length = rundata["cluster_length"]

        return demuxes

    def _recalculate_undetermined(self):
        # We have to calculate "corrected" unknown read counts when parsing more than one bclconvert run. To do this:
        # Add up all the reads in a lane that were assigned to samples, then take the total reads in a lane (which is taken _from the sum of all reads in a single file_),
        # subtract the former from the latter, and use that as "undetermined samples in lane."
        total_reads_per_lane = dict()
        for filename, lanedata in self.total_reads_in_lane_per_file.items():
            for lane_id, reads in lanedata.items():
                if lane_id not in total_reads_per_lane:
                    total_reads_per_lane[lane_id] = int()
                if total_reads_per_lane[lane_id] and reads != total_reads_per_lane[lane_id]:
                    log.error(
                        "Warning: different amounts of reads per lane across input files! Cannot expect calculatons to be accurate!"
                    )
                total_reads_per_lane[lane_id] = reads

        run_data = self.bclconvert_data[self.last_run_id]  # in this situaiton we have only one run id
        for lane_id, lane in run_data.items():
            determined_reads = 0
            for sample_id, sample in lane["samples"].items():
                determined_reads += sample["reads"]
            self.per_lane_undetermined_reads[lane_id] = total_reads_per_lane[lane_id] - determined_reads

    def parse_demux_data(self, myfile):
        # parse a bclconvert output stats csv, populate variables appropriately
        filename = str(os.path.join(myfile["root"], myfile["fn"]))
        self.total_reads_in_lane_per_file[filename] = dict()

        # reader = csv.DictReader(myfile["f"], delimiter=",")
        reader = csv.DictReader(open(filename), delimiter=",")
        for row in reader:
            run_data = self.bclconvert_data[myfile["run_id"]]
            lane_id = "L{}".format(row["Lane"])
            if lane_id not in run_data:
                run_data[lane_id] = {}
                run_data[lane_id]["samples"] = {}
                self._set_up_reads_dictionary(run_data[lane_id])
                self.per_lane_undetermined_reads[lane_id] = 0
            lane = run_data[lane_id]

            sample = row["SampleID"]

            if sample != "Undetermined":  # dont include undetermined reads at all in any of the calculations...
                if sample not in run_data[lane_id]["samples"]:
                    run_data[lane_id]["samples"][sample] = {}
                    self._set_up_reads_dictionary(run_data[lane_id]["samples"][sample])
                    run_data[lane_id]["samples"][sample]["filename"] = os.path.join(myfile["root"], myfile["fn"])
                    run_data[lane_id]["samples"][sample]["samples"] = {}

                lane_sample = run_data[lane_id]["samples"][sample]  # this sample in this lane

                # total lane stats
                lane["reads"] += int(row["# Reads"])
                lane["yield"] += int(row["# Reads"]) * (self.cluster_length)
                lane["perfect_index_reads"] += int(row["# Perfect Index Reads"])
                lane["one_mismatch_index_reads"] += int(row["# One Mismatch Index Reads"])
                lane["basesQ30"] += int(row["# of >= Q30 Bases (PF)"])

                # stats for this sample in this lane
                lane_sample["reads"] += int(row["# Reads"])
                lane_sample["yield"] += int(row["# Reads"]) * (self.cluster_length)
                lane_sample["perfect_index_reads"] += int(row["# Perfect Index Reads"])
                lane_sample["one_mismatch_index_reads"] += int(row["# One Mismatch Index Reads"])
                lane_sample["basesQ30"] += int(row["# of >= Q30 Bases (PF)"])
                lane_sample["mean_quality"] += float(row["Mean Quality Score (PF)"])

            if lane_id not in self.total_reads_in_lane_per_file[filename]:
                self.total_reads_in_lane_per_file[filename][lane_id] = 0

            self.total_reads_in_lane_per_file[filename][lane_id] += int(
                row["# Reads"]
            )  # add up number of reads, regardless of undetermined or not

            if self.num_demux_files == 1 and sample == "Undetermined":
                self.per_lane_undetermined_reads[lane_id] += int(row["# Reads"])

    def _parse_top_unknown_barcodes(self):
        run_data = self.bclconvert_data[self.last_run_id]

        for unknown_barcode_file in self.find_log_files("bclconvert/unknown_barcodes", filehandles=True):
            barcode_reader = csv.DictReader(unknown_barcode_file["f"], delimiter=",")
            for unknown_barcode_row in barcode_reader:
                thislane = "L" + str(unknown_barcode_row["Lane"])
                if "top_unknown_barcodes" not in run_data[thislane]:
                    run_data[thislane]["top_unknown_barcodes"] = {}
                thisbarcode = str(unknown_barcode_row["index"]) + "-" + str(unknown_barcode_row["index2"])
                run_data[thislane]["top_unknown_barcodes"][thisbarcode] = int(unknown_barcode_row["# Reads"])

    def _total_reads_for_run(self, run_id):
        totalreads = 0
        for lane_id, lane in self.bclconvert_data[run_id].items():
            totalreads += lane["reads"]
        return totalreads

    def _total_reads_all_runs(self):
        totalreads = 0
        for key, run_data in self.bclconvert_data.items():
            for lane_id, lane in run_data.items():
                totalreads += lane["reads"]
        return totalreads

    def _set_lane_percentage_stats(self, data):
        try:
            data["percent_Q30"] = (float(data["basesQ30"]) / float(data["reads"] * (self.cluster_length))) * 100.0
        except ZeroDivisionError:
            data["percent_Q30"] = "NA"
        try:
            data["percent_perfectIndex"] = (float(data["perfect_index_reads"]) / float(data["reads"])) * 100.0
        except ZeroDivisionError:
            data["percent_perfectIndex"] = "NA"
        try:
            data["percent_oneMismatch"] = float(data["one_mismatch_index_reads"]) / float(data["reads"]) * 100.0
        except ZeroDivisionError:
            data["percent_oneMismatch"] = "NA"
        try:
            data["depth"] = float(data["basesQ30"]) / float(self._get_genome_size())
        except ZeroDivisionError:
            data["depth"] = "NA"
        except TypeError:
            data["depth"] = "NA"

    def _split_data_by_lane_and_sample(self):
        # populate a collection of "stats across all lanes" and "stats across all sampeles"

        for run_id, r in self.bclconvert_data.items():
            # set stats for each lane (across all samples) in bclconvert_bylane dictionary
            for lane_id, lane in r.items():
                self._set_lane_percentage_stats(lane)

                lane_key_name = self.prepend_runid(run_id, lane_id)
                self.bclconvert_bylane[lane_key_name] = {
                    "depth": lane["depth"],
                    "reads": lane["reads"],
                    "yield": lane["yield"],
                    "perfect_index_reads": lane["perfect_index_reads"],
                    "one_mismatch_index_reads": lane["one_mismatch_index_reads"],
                    "basesQ30": lane["basesQ30"],
                    "percent_Q30": lane["percent_Q30"],
                    "percent_perfectIndex": lane["percent_perfectIndex"],
                    "percent_oneMismatch": lane["percent_oneMismatch"],
                    "top_unknown_barcodes": self.get_unknown_barcodes(lane["top_unknown_barcodes"])
                    if "top_unknown_barcodes" in lane
                    else {},
                }

                # now set stats for each sample (across all lanes) in bclconvert_bysample dictionary
                for sample_id, sample in lane["samples"].items():
                    if sample_id not in self.bclconvert_bysample:
                        self.bclconvert_bysample[sample_id] = {}
                        self._set_up_reads_dictionary(self.bclconvert_bysample[sample_id])

                    s = self.bclconvert_bysample[sample_id]

                    s["reads"] += int(sample["reads"])
                    s["yield"] += int(sample["yield"])
                    s["perfect_index_reads"] += int(sample["perfect_index_reads"])
                    s["one_mismatch_index_reads"] += int(sample["one_mismatch_index_reads"])
                    s["basesQ30"] += int(sample["basesQ30"])
                    s["mean_quality"] += float(sample["mean_quality"])

                    try:
                        if "depth" not in s:
                            s["depth"] = 0
                        s["depth"] += float(sample["basesQ30"]) / float(self._get_genome_size())
                    except ZeroDivisionError:
                        s["depth"] = "NA"
                    except TypeError:
                        s["depth"] = "NA"

                    if sample_id not in ["top_unknown_barcodes"]:
                        if sample_id not in self.source_files:
                            self.source_files[sample_id] = []
                        self.source_files[sample_id].append(sample["filename"])

    def get_unknown_barcodes(self, lane_unknown_barcode):
        """Python 2.* dictionaries are not sorted.
        This function return an `OrderedDict` sorted by barcode count.
        """
        try:
            sorted_barcodes = OrderedDict(
                sorted(lane_unknown_barcode.items(), key=operator.itemgetter(1), reverse=True)
            )
        except AttributeError:
            sorted_barcodes = None
        return sorted_barcodes

    def sample_stats_table(self):
        sample_stats_data = dict()
        total_reads = self._total_reads_all_runs()

        for sample_id, sample in self.bclconvert_bysample.items():
            # percent stats for bclconvert-bysample i.e. stats for sample across all lanes
            try:
                perfect_percent = "{0:.1f}".format(float(100.0 * sample["perfect_index_reads"] / sample["reads"]))
            except ZeroDivisionError:
                perfect_percent = "0.0"
            try:
                one_mismatch_pecent = "{0:.1f}".format(
                    float(100.0 * sample["one_mismatch_index_reads"] / sample["reads"])
                )
            except ZeroDivisionError:
                one_mismatch_pecent = "0.0"

            try:
                yield_q30_percent = "{0:.1f}".format(float(100.0 * (sample["basesQ30"] / sample["yield"])))
            except ZeroDivisionError:
                yield_q30_percent = "0.0"  #

            try:
                percent_yield = (float(sample["yield"]) / float((total_reads) * (self.cluster_length))) * 100.0
            except ZeroDivisionError:
                percent_yield = "NA"

            try:
                percent_reads = (float(sample["reads"]) / float(total_reads)) * 100.0
            except ZeroDivisionError:
                percent_reads = "NA"

            sample_stats_data[sample_id] = {
                "depth": sample["depth"],
                "basesQ30": sample["basesQ30"],
                "reads": sample["reads"],
                "percent_reads": percent_reads,
                "yield": sample["yield"],
                "percent_yield": percent_yield,
                "yield_q30_percent": yield_q30_percent,
                # "perfect_index": samle['perfect_index_reads'], # don't need these
                # "one_mismatch_index_reads": sample['one_mismatch_index_reads'],
                "perfect_pecent": perfect_percent,
                "one_mismatch_pecent": one_mismatch_pecent,
            }

        headers = OrderedDict()
        if sample["depth"] != "NA":
            headers["depth"] = {
                "title": "Coverage",
                "description": (
                    "Estimated sequencing depth based on the number of bases with quality score greater or equal to Q30, "
                    "assuming the genome size is {}, as provided in config".format(self._get_genome_size())
                ),
                "min": 0,
                "suffix": "X",
                "scale": "BuPu",
            }

        headers["reads"] = {
            "title": "{} Clusters".format(config.read_count_prefix),
            "description": "Total number of clusters (read pairs) for this sample as determined by bclconvert demultiplexing ({})".format(
                config.read_count_desc
            ),
            "scale": "Blues",
            "shared_key": "read_count",
        }
        headers["yield"] = {
            "title": "Yield ({})".format(config.base_count_prefix),
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
            "title": "Bases ({}) &ge; Q30 (PF)".format(config.base_count_prefix),
            "description": "Number of bases with a Phred score of 30 or higher, passing filter ({})".format(
                config.base_count_desc
            ),
            "scale": "Blues",
            "shared_key": "base_count",
        }
        headers["yield_q30_percent"] = {
            "title": "% Bases &ge; Q30 (PF)",
            "description": "Percent of bases with a Phred score of 30 or higher, passing filter ({})".format(
                config.base_count_desc
            ),
            "scale": "Greens",
            "max": 100,
            "min": 0,
            "suffix": "%",
        }
        headers["perfect_pecent"] = {
            "title": "% Perfect Index",
            "description": "Percent of reads with perfect index (0 mismatches)",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        headers["one_mismatch_pecent"] = {
            "title": "% One Mismatch Index",
            "description": "Percent of reads with one mismatch index",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }

        # Table config
        table_config = {
            "namespace": "bclconvert",
            "id": "bclconvert-sample-stats-table",
            "table_title": "bclconvert Sample Statistics",
            "no_beeswarm": True,
        }

        return table.plot(sample_stats_data, headers, table_config)

    def lane_stats_table(self):
        for lane_id, lane in self.bclconvert_bylane.items():
            try:
                yield_q30_percent = "{0:.1f}".format(float(100.0 * (lane["basesQ30"] / lane["yield"])))
            except ZeroDivisionError:
                yield_q30_percent = "0.0"
            self.bclconvert_bylane[lane_id]["yield_q30_percent"] = yield_q30_percent

        headers = OrderedDict()
        if lane["depth"] != "NA":
            headers["depth-lane"] = {
                "title": "Coverage",
                "description": (
                    "Estimated sequencing depth based on the number of bases with quality score greater or equal to Q30, "
                    "assuming the genome size is {}, as provided in config".format(self._get_genome_size())
                ),
                "suffix": "X",
                "scale": "BuPu",
            }

        headers["reads-lane"] = {
            "title": "{} Clusters".format(config.read_count_prefix),
            "description": "Total number of clusters (read pairs) for this sample as determined by bclconvert demultiplexing ({})".format(
                config.read_count_desc
            ),
            "scale": "Blues",
            "shared_key": "read_count",
        }
        headers["yield-lane"] = {
            "title": "Yield ({})".format(config.base_count_prefix),
            "description": "Total number of bases for this sample as determined by bclconvert demultiplexing ({})".format(
                config.base_count_desc
            ),
            "scale": "Greens",
            "shared_key": "base_count",
        }
        headers["basesQ30-lane"] = {
            "title": "Bases ({}) &ge; Q30 (PF)".format(config.base_count_prefix),
            "description": "Number of bases with a Phred score of 30 or higher, passing filter ({})".format(
                config.base_count_desc
            ),
            "scale": "Blues",
            "shared_key": "base_count",
        }
        headers["yield_q30_percent-lane"] = {
            "title": "% Bases &ge; Q30 (PF)",
            "description": "Percent of bases with a Phred score of 30 or higher, passing filter",
            "max": 100,
            "min": 0,
            "scale": "Greens",
        }
        headers["perfect_index_reads-lane"] = {
            "title": "{} Perfect Index".format(config.read_count_prefix),
            "description": "Reads with perfect index - 0 mismatches ({})".format(config.read_count_desc),
            "scale": "Blues",
            "shared_key": "read_count",
        }

        headers["one_mismatch_index_reads-lane"] = {
            "title": "{} One Mismatch".format(config.read_count_prefix),
            "description": "Reads with one mismatch index ({})".format(config.read_count_desc),
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
            "title": "% One Mismach",
            "description": "Percent of reads with one mismatch",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
            "suffix": "%",
        }

        # Table config
        table_config = {
            "namespace": "bclconvert-lane",
            "id": "bclconvert-lane-stats-table",
            "table_title": "bclconvert Lane Statistics",
            "col1_header": "Run ID - Lane",
            "no_beeswarm": True,
        }

        # new dict with matching keys for plotting (this avoids duplicate html id linting errors)
        bclconvert_bylane_foroutput = dict()
        for laneid, lanestats in self.bclconvert_bylane.items():
            if laneid not in bclconvert_bylane_foroutput:
                bclconvert_bylane_foroutput[laneid] = dict()
            for key, value in lanestats.items():
                bclconvert_bylane_foroutput[laneid][key + "-lane"] = value

        return table.plot(bclconvert_bylane_foroutput, headers, table_config)

    def prepend_runid(self, runId, rest):
        return str(runId) + " - " + str(rest)

    def get_bar_data_from_counts(self, counts):
        # for per-lane stats we fetch undetermied reads, too
        bar_data = {}
        for key, value in counts.items():
            bar_data[key] = {
                "perfect": value["perfect_index_reads"],
                "imperfect": value["reads"] - value["perfect_index_reads"],
            }
            try:
                if key.startswith(
                    self.prepend_runid(self.last_run_id, "")  # this wont run in multiple sequencing run situations
                ):  # per-lane stats start with a prepended run id, this is a per-lane entry
                    this_lane_id = key.replace(self.prepend_runid(self.last_run_id, ""), "")
                    rundata = self.bclconvert_data[self.last_run_id]
                    if this_lane_id in rundata:  # this is definitely a lane
                        bar_data[key]["undetermined"] = self.per_lane_undetermined_reads[this_lane_id]
            except TypeError:
                # do nothing, there is no Undetermined - this will happen in case of multiple run ids
                pass

        return bar_data

    def get_bar_data_from_undetermined(self, flowcells):
        # Get data to plot for undetermined barcodes.

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

        # sort results
        bar_data = OrderedDict(sorted(bar_data.items(), key=lambda x: sum(x[1].values()), reverse=True))
        return OrderedDict((key, value) for key, value in islice(bar_data.items(), 20))
