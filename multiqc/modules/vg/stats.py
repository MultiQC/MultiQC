"""MultiQC module to parse output from vg stats"""

import logging

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import violin

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
          name='VG',
          anchor='vg',
          href="https://github.com/vgteam/vg",
          info="is a tool to manipulate and analyze graphical genomes, including align reads to graphical genomes.",
          doi="10.1126/science.abg8871",
        #   comment="Summarizes the output of vg stats -a aligned.gam"
        )

        # Parse vg stats data
        self.vg_stats_map = dict()
        for f in self.find_log_files("vg/vg_stats", filehandles=True):
            self.vg_stats_map.update(self.get_vg_stats_data(f))
            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        self.vg_stats_map = self.ignore_samples(self.vg_stats_map)

        if len(self.vg_stats_map) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.vg_stats_map)} report{'s' if len(self.vg_stats_map) > 1 else ''}")

        # Write parsed report data to a file
        self.write_data_file(self.vg_stats_map, "multiqc_vg_stats")

        # Add data to summary table
        headers = {
            'Percent Aligned': {
                'title': '% Aligned (vg gam)',
                'description': 'Percentage of total reads aligned by vg giraffe to pangenomic reference',
                'max': 100,
                'min': 0,
                "scale": "RdYlGn",
                'suffix': '%'
            },
            'Percent Properly Paired': {
                'title': '% Properly Paired (vg gam)',
                'description': 'Percentage of graph-aligned reads in gam file that are properly paired',
                'max': 100,
                'min': 0,
                "scale": "RdYlGn",
                'suffix': '%'
            },
            'Mapping quality': {
                'title': 'Mapping quality (vg gam)',
                'description': 'Average mapping quality of graph-aligned reads in gam file',
                'max': 60,
                'min': 0,
                'scale': 'RdYlGn'
            }
        }

        # Add % Aligned, % Properly Paired, and average MQ to summary table
        self.general_stats_addcols(self.vg_stats_map, headers)

        # Add violin plot section to report, plotting major vg stats metrics
        self.make_violin_plot()


    def get_vg_stats_data(self, f):
        data = self.parse_vg_stats_file(f)

        # Calculate additional metrics
        try:
            data["Percent Aligned"] = (data["Total aligned"] / data["Total alignments"]) * 100.0
            data["Percent Properly Paired"] = (data["Total properly paired"] / data["Total alignments"]) * 100.0
        except ZeroDivisionError:
            data["Percent Aligned"] = 0.0
            data["Percent Properly Paired"] = 0.0
        except KeyError:
            pass

        # Add the data to the main vg_stats_map dict
        if f['s_name'] in self.vg_stats_map:
            log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

        return {f['s_name']: data}


    def parse_vg_stats_file(self, f):
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
        for line in f['f']:
            var=list()
            val=list()
            s = line.strip().split(': ')
            if s[1].isnumeric():
                var.append(s[0])
                val.append(s[1])
            elif 'mean' in line:
                var.append(s[0])
                val.append(s[1].split(',')[0].split(' ')[1])
            elif 'read events' in line:
                var.append(s[0]+' (bp)')
                val.append(s[1].split(' bp ')[0])
                var.append(s[0]+' (reads)')
                val.append(s[1].split(' bp in ')[1].split(' ')[0])
            elif s[1].split(' ')[0].isnumeric():
                unit = s[1].split(' ')[1]
                var.append(s[0] + ' (' + unit + ')')
                val.append(s[1].split(' ')[0])
            else:
                continue
            for k, v in zip(var, val):
                assert v.replace('.','',1).isdigit()
                data[k] = float(v)
        return data


    def make_violin_plot(self):
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
            **{"title": "Total reads perfectly aligned", "description": "Properly paired reads aligned with no indels, soft clips, or mismatches"},
        )
        table_column_metadata["Total properly paired"] = dict(reads, **{"title": "Properly paired", "description": "Read pair is mapped to same contig, properly oriented, and within 6 SD of mean insert length."})
        table_column_metadata["Insertions (reads)"] = dict(reads, **{"title": "Insertions (reads)"})
        table_column_metadata["Deletions (reads)"] = dict(reads, **{"title": "Deletions (reads)"})
        table_column_metadata["Substitutions (reads)"] = dict(reads, **{"title": "Substitutions (reads)"})
        table_column_metadata["Softclips (reads)"] = dict(reads, **{"title": "Softclips (reads)"})
        table_column_metadata["Insertions (bp)"] = dict(
            bases, **{"title": "Insertions (bases)"}
        )
        table_column_metadata["Deletions (bp)"] = dict(bases, **{"title": "Deletions (bases)"})
        table_column_metadata["Substitutions (bp)"] = dict(bases, **{"title": "Substitutions (bases)"})
        table_column_metadata["Softclips (bp)"] = dict(bases, **{"title": "Bases Softclipped"})
        table_column_metadata["Percent Aligned"] = dict(**{"title": '% reads aligned', "suffix":'%', 'shared_key':'vg_read_percentage'})
        table_column_metadata["Percent Properly Paired"] = dict(**{"title": "'% reads properly paired'", "suffix":'%', 'shared_key':'vg_read_percentage'})
        table_column_metadata["Alignment score"] = dict(**{"title": "Mean alignment score", "description": "Smith-Waterman alignment score",'min':0, 'max':150})
        table_column_metadata["Mapping quality"] = dict(**{"title": "Mean MQ", "description": "Read MQ (Li and Durbin method)", 'min':0, 'max':60})

        ## Data to report out in table format, but not showed by default in report
        table_column_metadata["Total primary"] = dict(reads, **{"title": "Total primary alignments", 'hidden':True})
        table_column_metadata["Total secondary"] = dict(reads, **{"title": "Total secondary alignments", 'hidden':True})
        table_column_metadata["Total paired"] = dict(reads, **{"title": "Total paired", "description": "Total number of reads in a read pair", 'hidden':True})
        table_column_metadata["Total gapless (softclips allowed)"] = dict(reads, **{"title": "Total gapless (softclips allowed)", "description": "Total reads without indels or substitutions", 'hidden':True})
        table_column_metadata["Total time (seconds)"] = dict(hours, **{"title": "Total running time (cpu-hours)", "description":"Total cpu-hours per sample", 'hidden':True})
        # print (self.vg_stats_map.keys())
        print ([x for x in self.vg_stats_map['HSB8153'].keys() if x not in table_column_metadata.keys()])
        print ([x for x in table_column_metadata.keys() if x not in self.vg_stats_map['HSB8153'].keys()])
        self.add_section(
            name="Alignment stats",
            anchor="vg stats",
            description="This module parses the output from <code>vg stats</code>. All numbers in millions.",
            helptext="These data represent alignment statistics calculated from a gam file of reads aligned to a graphical genome. The process of surjecting these alignments to a linear reference to produce a sam/bam aligned file results in the loss of aligned reads. (For example, a read that aligns to a non-reference insertion has no aligned location in the linear reference genome.) As result, these statistics will differ from those produced by a samtools stats report of the surjected bam file.",
            plot=violin.plot(
                data=self.vg_stats_map,
                headers=table_column_metadata,
                pconfig={
                    "id": "vg_stats",
                    "title": "VG stats: Graphical Genome Alignment Stats",
                    "save_file": False
                },
            )
        )