import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import violin, bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses the
    [vg stats](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe#evaluating-with-vg-stats)
    reports that summarize the stats of read alignment to a graphical genome in a GAM file.

    `vg stats` is capable of producing many reports summarizing many aspects of graphical
    genomes, including specific aspects of aligned GAM files such as node coverage.
    This module is not meant to gather those data. Rather, this module is designed to
    summarize the alignment performance of GAM files produced by `vg giraffe` created
    from the stdout of the `vg stats` command:

    ```sh
    $ vg stats -a mapped.gam > sample-stats.txt
    $ cat sample-stats.txt
    Total alignments: 727413268
    Total primary: 727413268
    Total secondary: 0
    Total aligned: 717826332
    Total perfect: 375143620
    Total gapless (softclips allowed): 714388968
    Total paired: 727413268
    Total properly paired: 715400510
    Alignment score: mean 129.012, median 132, stdev 31.5973, max 161 (244205781 reads)
    Mapping quality: mean 52.8552, median 60, stdev 17.7742, max 60 (589259353 reads)
    Insertions: 3901467 bp in 1466045 read events
    Deletions: 6759252 bp in 2795331 read events
    Substitutions: 281648245 bp in 281648245 read events
    Softclips: 11480269152 bp in 252773804 read events
    Total time: 291465 seconds
    Speed: 2495.71 reads/second
    ```

    It is not guaranteed that output created using any other parameter combination can
    be parsed using this module.

    The graphical reports are designed to mimic a samtools stats report, including:

    1. A bar chart showing the breakdown of aligned, perfectly aligned, and unaligned reads.
    2. A violin plot for all metrics.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="VG",
            anchor="vg",
            href="https://github.com/vgteam/vg",
            info="Toolkit to manipulate and analyze graphical genomes, including read alignment",
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
        self.general_stats_addcols(
            data_by_sample,
            {
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
            },
        )

        self.make_barplot(data_by_sample)

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

    def make_barplot(self, data_by_sample: Dict):
        cats = {
            "Perfect": {"color": "#437bb1", "name": "Perfect"},
            "Imperfect": {"color": "#FF9933", "name": "Imperfect"},
            "Unaligned": {"color": "#b1084c", "name": "Unaligned"},
        }

        bar_data_by_sample: Dict = {
            s: {
                "Perfect": d["Total perfect"],
                "Imperfect": d["Total aligned"] - d["Total perfect"],
                "Unaligned": d["Total alignments"] - d["Total aligned"],
            }
            for s, d in data_by_sample.items()
        }
        self.add_section(
            name="Reads aligned",
            anchor="vg_stats_aligned",
            description="Read alignment breakdown from <code>vg stats</code> as calculated from the GAM file produced by `vg giraffe`",
            helptext="""
"Imperfect" is calculated as the difference between "Total aligned" and "Total perfect",
and "Unaligned" is calculated as the difference between "Total alignments" and "Total aligned" reads.
""",
            plot=bargraph.plot(
                data=bar_data_by_sample,
                cats=cats,
                pconfig={
                    "id": "vg_stats_aligned_barplot",
                    "title": "VG stats: alignment scores",
                },
            ),
        )

    def make_violin_plot(self, data_by_sample: Dict):
        # Make dot plot of counts
        headers = {
            "Total alignments": {
                "title": "Total reads",
                "shared_key": "read_count",
            },
            "Total aligned": {
                "title": "Reads aligned",
                "shared_key": "read_count",
            },
            "Total perfect": {
                "title": "Total reads perfectly aligned",
                "description": "Properly paired reads aligned with no indels, soft clips, or mismatches",
                "shared_key": "read_count",
            },
            "Total properly paired": {
                "title": "Properly paired",
                "description": "Read pair is mapped to same contig, properly oriented, and within 6 SD of mean insert length.",
                "shared_key": "read_count",
            },
            "Insertions (reads)": {"title": "Insertions (reads)", "shared_key": "read_count"},
            "Deletions (reads)": {"title": "Deletions (reads)", "shared_key": "read_count"},
            "Substitutions (reads)": {"title": "Substitutions (reads)", "shared_key": "read_count"},
            "Softclips (reads)": {"title": "Softclips (reads)", "shared_key": "read_count"},
            "Insertions (bp)": {"title": "Insertions (bases)", "shared_key": "base_count"},
            "Deletions (bp)": {"title": "Deletions (bases)", "shared_key": "base_count"},
            "Substitutions (bp)": {"title": "Substitutions (bases)", "shared_key": "base_count"},
            "Softclips (bp)": {"title": "Bases soft-clipped", "shared_key": "base_count"},
            "Percent Aligned": {
                "title": "% reads aligned",
                "suffix": "%",
                "shared_key": "vg_read_percentage",
            },
            "Percent Properly Paired": {
                "title": "'% reads properly paired'",
                "suffix": "%",
                "shared_key": "vg_read_percentage",
            },
            "Alignment score": {
                "title": "Mean alignment score",
                "description": "Smith-Waterman alignment score",
                "min": 0,
                "max": 150,
            },
            "Mapping quality": {
                "title": "Mean MQ",
                "description": "Read MQ (Li and Durbin method)",
                "min": 0,
                "max": 60,
            },
            "Total primary": {
                "title": "Total primary alignments",
                "hidden": True,
                "shared_key": "read_count",
            },
            "Total secondary": {
                "title": "Total secondary alignments",
                "hidden": True,
                "shared_key": "read_count",
            },
            "Total paired": {
                "title": "Total paired",
                "description": "Total number of reads in a read pair",
                "hidden": True,
                "shared_key": "read_count",
            },
            "Total gapless (softclips allowed)": {
                "title": "Total gapless (softclips allowed)",
                "description": "Total reads without indels or substitutions",
                "hidden": True,
                "shared_key": "read_count",
            },
            "Total time (seconds)": {
                "min": 0,
                "modify": lambda x: float(x) / 3600,
                "suffix": "h",
                "tt_decimals": 2,
                "title": "Total running time (cpu-hours)",
                "description": "Total cpu-hours per sample",
                "hidden": True,
            },
        }

        # Data to report out in table format, but not showed by default in report

        self.add_section(
            name="Alignment stats",
            anchor="vg_stats",
            description="Alignment metrics from <code>vg stats</code> as calculated from the GAM file produced by `vg giraffe`",
            helptext="These data represent alignment statistics calculated from a GAM file of reads aligned to a graphical genome. The process of mapping these alignments to a linear reference to produce a SAM/BAM aligned file results in the loss of aligned reads. (For example, a read that aligns to a non-reference insertion has no aligned location in the linear reference genome.) As result, these statistics will differ from those produced by a samtools stats report of the mapping BAM file.",
            plot=violin.plot(
                data=data_by_sample,
                headers=headers,
                pconfig={
                    "id": "vg_stats_violin",
                    "title": "VG stats: Graphical Genome Alignment Stats",
                },
            ),
        )
