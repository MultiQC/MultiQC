"""MultiQC module to parse output from vg stats"""

import logging
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import violin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="VG",
            anchor="vg",
            href="https://github.com/vgteam/vg",
            info="is a tool to manipulate and analyze graphical genomes",
            doi="10.1038/nbt.4227",
        )

        # Parse vg stats data
        data_by_sample: Dict = dict()
        for f in self.find_log_files("vg/stats", filehandles=True):
            data = self.parse_vg_stats_file(f)
            if f["s_name"] in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

            # Calculate additional metrics
            try:
                data["Percent Aligned"] = (data["Total aligned"] / data["Total alignments"]) * 100.0
                data["Percent Properly Paired"] = (data["Total properly paired"] / data["Total alignments"]) * 100.0
            except ZeroDivisionError:
                data["Percent Aligned"] = 0.0
                data["Percent Properly Paired"] = 0.0
            except KeyError:
                pass

            data_by_sample[f["s_name"]] = data

            self.add_data_source(f)
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        log.info(f"Found {len(data_by_sample)} report{'s' if len(data_by_sample) > 1 else ''}")

        # Write parsed report data to a file
        self.write_data_file(data_by_sample, "multiqc_vg_stats")

        # Add data to summary table
        headers = {
            "Percent Aligned": {
                "title": "Aligned",
                "description": "Percentage of total reads aligned by vg giraffe to pangenomic reference",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "Percent Properly Paired": {
                "title": "Properly paired",
                "description": "Percentage of graph-aligned reads in a GAM file that are properly paired",
                "max": 100,
                "min": 0,
                "scale": "Blues",
                "suffix": "%",
            },
            "Mapping quality": {
                "title": "MQ",
                "description": "Average mapping quality of graph-aligned reads in a GAM file",
                "max": 60,
                "min": 0,
                "scale": "BuGn",
            },
        }

        # Add % Aligned, % Properly Paired, and average MQ to summary table
        self.general_stats_addcols(data_by_sample, headers)

        # Add violin plot section to report, plotting major vg stats metrics
        self.make_violin_plot(data_by_sample)

    @staticmethod
    def parse_vg_stats_file(f):
        """
        Load the vg stats file.
        The stats reports assume the following structure:
        ---
        Total alignments: 786889820
        Total primary: 786889820
        Total secondary: 0
        Total aligned: 777651532
        Total perfect: 382250293
        Total gapless (softclips allowed): 772411557
        Total paired: 786889820
        Total properly paired: 774424970
        Alignment score: mean 141.081, median 156, stdev 30.4529, max 161 (382250293 reads)
        Mapping quality: mean 52.9098, median 60, stdev 17.7611, max 60 (643610302 reads)
        Insertions: 6503589 bp in 2310323 read events
        Deletions: 10842444 bp in 4336335 read events
        Substitutions: 544551412 bp in 544551412 read events
        Softclips: 11476128752 bp in 248401284 read events
        Total time: 771240 seconds
        Speed: 1020.29 reads/second
        """
        # Load the file
        data = {}
        for line in f["f"]:
            var = list()
            val = list()
            s = line.strip().split(": ")
            if s[1].isnumeric():
                var.append(s[0])
                val.append(s[1])
            elif "mean" in line:
                var.append(s[0])
                val.append(s[1].split(",")[0].split(" ")[1])
            elif "read events" in line:
                var.append(s[0] + " (bp)")
                val.append(s[1].split(" bp ")[0])
                var.append(s[0] + " (reads)")
                val.append(s[1].split(" bp in ")[1].split(" ")[0])
            elif s[1].split(" ")[0].isnumeric():
                unit = s[1].split(" ")[1]
                var.append(s[0] + " (" + unit + ")")
                val.append(s[1].split(" ")[0])
            else:
                continue
            for k, v in zip(var, val):
                assert v.replace(".", "", 1).isdigit()
                data[k] = float(v)
        return data

    def make_violin_plot(self, data_by_sample: Dict):
        # Make dot plot of counts
        table_column_metadata = {}
        reads = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": config.read_count_prefix,
            "tt_decimals": 2,
            "shared_key": "read_count",
        }
        bases = {
            "min": 0,
            "modify": lambda x: float(x) * config.base_count_multiplier,
            "suffix": config.base_count_prefix,
            "tt_decimals": 2,
            "shared_key": "base_count",
        }
        hours = {
            "min": 0,
            "modify": lambda x: float(x) / 3600,
            "suffix": "h",
            "tt_decimals": 2,
            "shared_key": "cpu_hours",
        }

        table_column_metadata["Total alignments"] = dict(reads, **{"title": "Total reads"})
        table_column_metadata["Total aligned"] = dict(reads, **{"title": "Reads aligned"})
        table_column_metadata["Total perfect"] = dict(
            reads,
            **{
                "title": "Total reads perfectly aligned",
                "description": "Properly paired reads aligned with no indels, soft clips, or mismatches",
            },
        )
        table_column_metadata["Total properly paired"] = dict(
            reads,
            **{
                "title": "Properly paired",
                "description": "Read pair is mapped to same contig, properly oriented, and within 6 SD of mean insert length.",
            },
        )
        table_column_metadata["Insertions (reads)"] = dict(reads, **{"title": "Insertions (reads)"})
        table_column_metadata["Deletions (reads)"] = dict(reads, **{"title": "Deletions (reads)"})
        table_column_metadata["Substitutions (reads)"] = dict(reads, **{"title": "Substitutions (reads)"})
        table_column_metadata["Softclips (reads)"] = dict(reads, **{"title": "Softclips (reads)"})
        table_column_metadata["Insertions (bp)"] = dict(bases, **{"title": "Insertions (bases)"})
        table_column_metadata["Deletions (bp)"] = dict(bases, **{"title": "Deletions (bases)"})
        table_column_metadata["Substitutions (bp)"] = dict(bases, **{"title": "Substitutions (bases)"})
        table_column_metadata["Softclips (bp)"] = dict(bases, **{"title": "Bases Softclipped"})
        table_column_metadata["Percent Aligned"] = dict(
            **{"title": "% reads aligned", "suffix": "%", "shared_key": "vg_read_percentage"}
        )
        table_column_metadata["Percent Properly Paired"] = dict(
            **{"title": "'% reads properly paired'", "suffix": "%", "shared_key": "vg_read_percentage"}
        )
        table_column_metadata["Alignment score"] = dict(
            **{"title": "Mean alignment score", "description": "Smith-Waterman alignment score", "min": 0, "max": 150}
        )
        table_column_metadata["Mapping quality"] = dict(
            **{"title": "Mean MQ", "description": "Read MQ (Li and Durbin method)", "min": 0, "max": 60}
        )

        ## Data to report out in table format, but not showed by default in report
        table_column_metadata["Total primary"] = dict(reads, **{"title": "Total primary alignments", "hidden": True})
        table_column_metadata["Total secondary"] = dict(
            reads, **{"title": "Total secondary alignments", "hidden": True}
        )
        table_column_metadata["Total paired"] = dict(
            reads, **{"title": "Total paired", "description": "Total number of reads in a read pair", "hidden": True}
        )
        table_column_metadata["Total gapless (softclips allowed)"] = dict(
            reads,
            **{
                "title": "Total gapless (softclips allowed)",
                "description": "Total reads without indels or substitutions",
                "hidden": True,
            },
        )
        table_column_metadata["Total time (seconds)"] = dict(
            hours,
            **{"title": "Total running time (cpu-hours)", "description": "Total cpu-hours per sample", "hidden": True},
        )

        self.add_section(
            name="Alignment stats",
            # anchor="vg_stats",
            description="This module parses the output from <code>vg stats</code>. All numbers in millions.",
            helptext="These data represent alignment statistics calculated from a gam file of reads aligned to a graphical genome. The process of surjecting these alignments to a linear reference to produce a sam/bam aligned file results in the loss of aligned reads. (For example, a read that aligns to a non-reference insertion has no aligned location in the linear reference genome.) As result, these statistics will differ from those produced by a samtools stats report of the surjected bam file.",
            plot=violin.plot(
                data=data_by_sample,
                headers=table_column_metadata,
                pconfig={
                    "id": "vg_stats",
                    "title": "VG stats: Graphical Genome Alignment Stats",
                },
            ),
        )
