""" MultiQC module to parse output from RSEM/rsem-calculate-expression """


import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph
from multiqc.utils import config

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    RSEM module class, parses .cnt file .
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Rsem",
            anchor="rsem",
            href="https://deweylab.github.io/RSEM/",
            info="RSEM (RNA-Seq by Expectation-Maximization) is a software package for"
            "estimating gene and isoform expression levels from RNA-Seq data.",
            doi="10.1186/1471-2105-12-323",
        )

        self.rsem_mapped_data = dict()
        self.rsem_multimapping_data = dict()

        # Find and load any count file
        for f in self.find_log_files("rsem"):
            self.parse_rsem_report(f)

        # Filter to strip out ignored sample names
        self.rsem_mapped_data = self.ignore_samples(self.rsem_mapped_data)

        if len(self.rsem_mapped_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.rsem_mapped_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.rsem_mapped_data, "multiqc_rsem")

        # Basic Stats Table
        self.rsem_stats_table()

        # Assignment bar plot
        self.rsem_mapped_reads_plot()

        # Multimapping line plot
        self.rsem_multimapping_plot()

    def parse_rsem_report(self, f):
        """Parse the rsem cnt stat file.
        Description of cnt file found : https://github.com/deweylab/RSEM/blob/master/cnt_file_description.txt
        """
        data = dict()
        multimapping_hist = dict()
        in_hist = False
        for line in f["f"].splitlines():
            s = line.split()
            if len(s) > 3:
                # Line: N0 N1 N2 N_tot
                # N0, number of unalignable reads;
                # N1, number of alignable reads; (nUnique + nMulti)
                # N2, number of filtered reads due to too many alignments
                # N_tot = N0 + N1 + N2
                data["Unalignable"] = int(s[0])
                data["Alignable"] = int(s[1])
                data["Filtered"] = int(s[2])
                data["Total"] = int(s[3])
                try:
                    data["alignable_percent"] = (float(s[1]) / float(s[3])) * 100.0
                except ZeroDivisionError:
                    data["alignable_percent"] = 0
            elif len(s) == 3:
                # Line: nUnique nMulti nUncertain
                # nUnique, number of reads aligned uniquely to a gene
                # nMulti, number of reads aligned to multiple genes; nUnique + nMulti = N1;
                # nUncertain, number of reads aligned to multiple locations in the given reference sequences, which include isoform-level multi-mapping reads
                data["Unique"] = int(s[0])
                data["Multi"] = int(s[1])
                data["Uncertain"] = int(s[2])
            elif len(s) == 2:
                # Number of gene alignments and read count
                if in_hist or int(s[0]) == 0:
                    in_hist = True
                    try:
                        multimapping_hist[int(s[0])] = int(s[1])
                    except ValueError:
                        pass
            else:
                break
        try:
            assert data["Unique"] + data["Multi"] == data["Alignable"]
        except AssertionError:
            log.warning(f"Unique + Multimapping read counts != alignable reads! '{f['fn']}'")
            return None
        except KeyError:
            log.warning(f"Error parsing RSEM counts file '{f['fn']}'")
            return None

        # Save parsed data
        if len(data) > 0:
            if f["s_name"] in self.rsem_mapped_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.rsem_mapped_data[f["s_name"]] = data
            self.add_data_source(f)
        if len(multimapping_hist) > 0:
            self.rsem_multimapping_data[f["s_name"]] = multimapping_hist

    def rsem_stats_table(self):
        """Take the parsed stats from the rsem report and add them to the
        basic stats table at the top of the report"""
        headers = {
            "alignable_percent": {
                "title": f"% Alignable, {config.read_count_prefix}",
                "description": f"% Alignable reads, {config.read_count_desc}",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            }
        }
        self.general_stats_addcols(self.rsem_mapped_data, headers)

    def rsem_mapped_reads_plot(self):
        """Make the rsem assignment rates plot"""

        # Plot categories
        keys = {
            "Unique": {"color": "#437bb1", "name": "Aligned uniquely to a gene"},
            "Multi": {"color": "#e63491", "name": "Aligned to multiple genes"},
            "Filtered": {"color": "#b1084c", "name": "Filtered due to too many alignments"},
            "Unalignable": {"color": "#7f0000", "name": "Unalignable reads"},
        }

        # Config for the plot
        config = {
            "id": "rsem_assignment_plot",
            "title": "RSEM: Mapped reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_zero_cats": False,
        }

        self.add_section(
            name="Mapped Reads",
            anchor="rsem_mapped_reads",
            description="A breakdown of how all reads were aligned for each sample.",
            plot=bargraph.plot(self.rsem_mapped_data, keys, config),
        )

    def rsem_multimapping_plot(self):
        """Make a line plot showing the multimapping levels"""

        pconfig = {
            "id": "rsem_multimapping_rates",
            "title": "RSEM: Multimapping Rates",
            "ylab": "Counts",
            "xlab": "Number of alignments",
            "xDecimals": False,
            "ymin": 0,
            "tt_label": "<b>{point.x} alignments</b>: {point.y:.0f}",
        }

        self.add_section(
            name="Multimapping rates",
            anchor="rsem_multimapping",
            description="A frequency histogram showing how many reads were aligned to `n` reference regions.",
            helptext="""In an ideal world, every sequence reads would align uniquely to a single location in the
                reference. However, due to factors such as repeititve sequences, short reads and sequencing errors,
                reads can be align to the reference 0, 1 or more times. This plot shows the frequency of each factor
                of multimapping. Good samples should have the majority of reads aligning once.""",
            plot=linegraph.plot(self.rsem_multimapping_data, pconfig),
        )
