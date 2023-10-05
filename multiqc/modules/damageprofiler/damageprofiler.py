""" MultiQC module to parse output from DamageProfiler """


import json
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        """
        damageprofiler module class
        """

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="DamageProfiler",
            anchor="damageprofiler",
            href="https://github.com/Integrative-Transcriptomics/DamageProfiler",
            info="a tool to determine damage patterns on ancient DNA.",
            doi="10.1093/bioinformatics/btab190",
        )

        # Init empty dictionaries

        self.threepGtoAfreq_data = dict()
        self.fivepCtoTfreq_data = dict()
        self.lgdist_fw_data = dict()
        self.lgdist_rv_data = dict()
        self.summary_metrics_data = dict()

        # Find and load JSON file
        for f in self.find_log_files("damageprofiler", filehandles=True):
            self.parseJSON(f)

        # Filter to strip out ignored sample names
        self.threepGtoAfreq_data = self.ignore_samples(self.threepGtoAfreq_data)
        self.fivepCtoTfreq_data = self.ignore_samples(self.fivepCtoTfreq_data)
        self.lgdist_fw_data = self.ignore_samples(self.lgdist_fw_data)
        self.lgdist_rv_data = self.ignore_samples(self.lgdist_rv_data)
        self.summary_metrics_data = self.ignore_samples(self.summary_metrics_data)

        if len(self.summary_metrics_data) == 0:
            raise ModuleNoSamplesFound

        # Write parsed report data to a file
        self.write_data_file(self.summary_metrics_data, "multiqc_damageprofiler_metrics")

        # Basic Stats Table, use generic function to add data to general table
        self.dmgprof_misinc_stats(self.threepGtoAfreq_data, "3 Prime", "G>A")
        self.dmgprof_misinc_stats(self.fivepCtoTfreq_data, "5 Prime", "C>T")
        self.addSummaryMetrics(self.summary_metrics_data)

        # Add plots
        if len(self.threepGtoAfreq_data) > 0:
            self.add_section(
                name="3P misincorporation plot",
                description="3' misincorporation plot for G>A substitutions",
                helptext="""
                This plot shows the frequency of G>A substitutions at the 3' read ends. Typically, one would observe high substitution percentages for ancient DNA, whereas modern DNA does not show these in higher extents.
                """,
                plot=self.threeprime_plot(),
            )
        if len(self.fivepCtoTfreq_data) > 0:
            self.add_section(
                name="5P misincorporation plot",
                description="5' misincorporation plot for C>T substitutions",
                helptext="""
                This plot shows the frequency of C>T substitutions at the 5' read ends. Typically, one would observe high substitution percentages for ancient DNA, whereas modern DNA does not show these in higher extents.
                """,
                plot=self.fiveprime_plot(),
            )
        if len(self.lgdist_fw_data) > 0 and len(self.lgdist_rv_data) > 0:
            self.add_section(
                name="Forward read length distribution",
                description="Read length distribution for forward strand (+) reads.",
                helptext="""
                This plot shows the read length distribution of the forward reads in the investigated sample. Reads below lengths of 30bp are typically filtered, so the plot doesn't show these in many cases. A shifted distribution of read lengths towards smaller read lengths (e.g around 30-50bp) is also an indicator of ancient DNA.
                """,
                plot=self.lgdistplot(self.lgdist_fw_data, "Forward"),
            )
            self.add_section(
                name="Reverse read length distribution",
                description="Read length distribution for reverse strand (-) reads.",
                helptext="""
                This plot shows the read length distribution of the reverse reads in the investigated sample. Reads below lengths of 30bp are typically filtered, so the plot doesn't show these in many cases. A shifted distribution of read lengths towards smaller read lengths (e.g around 30-50bp) is also an indicator of ancient DNA.
                """,
                plot=self.lgdistplot(self.lgdist_rv_data, "Reverse"),
            )

    # Parse our nice little JSON file
    def parseJSON(self, f):
        """Parse the JSON output from DamageProfiler and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except Exception as e:
            print(e)
            log.warning("Could not parse DamageProfiler JSON: '{}'".format(f["fn"]))
            return None

        # Get sample name from JSON first
        s_name = self.clean_s_name(parsed_json["metadata"]["sample_name"], f)
        self.add_data_source(f, s_name)

        # Add version info
        version = parsed_json["metadata"]["version"]
        self.add_software_version(version, s_name)

        # Add 3' G to A data
        self.threepGtoAfreq_data[s_name] = parsed_json["dmg_3p"]

        # Add 5' C to T data
        self.fivepCtoTfreq_data[s_name] = parsed_json["dmg_5p"]

        # Add lendist forward
        self.lgdist_fw_data[s_name] = parsed_json["lendist_fw"]

        # Add lendist reverse
        self.lgdist_rv_data[s_name] = parsed_json["lendist_rv"]

        # Add summary metrics
        self.summary_metrics_data[s_name] = parsed_json["summary_stats"]

    #### Tables from here on

    def dmgprof_misinc_stats(self, dict_to_plot, readend, substitution):
        """Take the parsed stats from the DamageProfiler and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["{}1".format(readend)] = {
            "id": "misinc-stats-1st-{}-{}".format(readend, substitution),
            "title": "{} {} 1st base".format(readend, substitution),
            "description": "{} 1st base substitution frequency for {}".format(readend, substitution),
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGnBu",
            "modify": lambda x: x * 100.0,
        }
        headers["{}2".format(readend)] = {
            "id": "misinc-stats-2nd-{}-{}".format(readend, substitution),
            "title": "{} {} 2nd base".format(readend, substitution),
            "description": "{} 2nd base substitution frequency for {}".format(readend, substitution),
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "BuGn",
            "hidden": True,
            "modify": lambda x: x * 100.0,
        }

        # Create new small subset dictionary for entries (we need just the first two data (k,v) pairs from each report)
        # Only the first two parts are informative from both 3' and 5' ends of reads (1st, 2nd base damage pattern)
        data = OrderedDict()
        dict_to_add = dict()

        for key in dict_to_plot.keys():
            tmp = dict_to_plot[key]
            pos = [readend + "1", readend + "2"]
            strlist = tmp[:2]
            tuples = list(zip(pos, strlist))
            data = dict((x, y) for x, y in tuples)
            # Extract first two elements from list
            dict_to_add[key] = data

        self.general_stats_addcols(dict_to_add, headers)

    def addSummaryMetrics(self, dict_to_plot):
        """Take the parsed stats from the DamageProfiler and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()

        headers["std"] = {
            "title": "Read length std. dev.",
            "description": "Read length std. dev.",
            "suffix": "bp",
            "scale": "PuBu",
            "format": "{:,.2f}",
            "shared_key": "read_length",
            "hidden": True,
        }

        headers["median"] = {
            "title": "Median read length",
            "description": "Median read length",
            "suffix": "bp",
            "scale": "YlGnBu",
            "format": "{:,.2f}",
            "shared_key": "read_length",
        }
        headers["mean_readlength"] = {
            "title": "Mean read length",
            "description": "Mean read length",
            "suffix": "bp",
            "scale": "PuBuGn",
            "format": "{:,.2f}",
            "shared_key": "read_length",
            "hidden": True,
        }

        self.general_stats_addcols(dict_to_plot, headers)

    #### Plotting from here on
    # Nice Linegraph plot for lgdist data

    def lgdistplot(self, dict_to_use, orientation):
        """Generate a read length distribution plot"""

        data = dict()
        for s_name in dict_to_use:
            try:
                data[s_name] = {int(d): int(dict_to_use[s_name][d]) for d in dict_to_use[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug("No valid data for forward read lgdist input!")
            return None

        config = {
            "id": "length-distribution-{}".format(orientation),
            "title": "DamageProfiler: Read length distribution - {} ".format(orientation),
            "ylab": "Number of reads",
            "xlab": "Readlength (bp)",
            "xDecimals": False,
            "tt_label": "{point.y} reads of length {point.x}",
            "ymin": 0,
            "xmin": 0,
        }
        return linegraph.plot(data, config)

    # Linegraph plot for 3pGtoA
    def threeprime_plot(self):
        """Generate a 3' G>A linegraph plot"""

        data = dict()
        dict_to_add = dict()
        # Create tuples out of entries
        for key in self.threepGtoAfreq_data:
            pos = list(range(1, len(self.threepGtoAfreq_data.get(key))))
            # Multiply values by 100 to get %
            tmp = [i * 100.0 for i in self.threepGtoAfreq_data.get(key)]
            tuples = list(zip(pos, tmp))
            # Get a dictionary out of it
            data = dict((x, y) for x, y in tuples)
            dict_to_add[key] = data

        config = {
            "id": "threeprime_misinc_plot",
            "title": "DamageProfiler: 3' G>A misincorporation plot",
            "ylab": "% G to A substituted",
            "xlab": "Nucleotide position from 3'",
            "tt_label": "{point.y:.2f} % G>A misincorporations at nucleotide position {point.x}",
            "ymin": 0,
            "xmin": 1,
        }

        return linegraph.plot(dict_to_add, config)

    # Linegraph plot for 5pCtoT
    def fiveprime_plot(self):
        """Generate a 5' C>T linegraph plot"""

        data = dict()
        dict_to_add = dict()
        # Create tuples out of entries
        for key in self.fivepCtoTfreq_data:
            pos = list(range(1, len(self.fivepCtoTfreq_data.get(key))))
            tmp = [i * 100.0 for i in self.fivepCtoTfreq_data.get(key)]
            tuples = list(zip(pos, tmp))
            # Get a dictionary out of it
            data = dict((x, y) for x, y in tuples)
            dict_to_add[key] = data

        config = {
            "id": "fiveprime_misinc_plot",
            "title": "DamageProfiler: 5' C>T misincorporation plot",
            "ylab": "% C to T substituted",
            "xlab": "Nucleotide position from 5'",
            "tt_label": "{point.y:.2f} % C>T misincorporations at nucleotide position {point.x}",
            "ymin": 0,
            "xmin": 1,
        }

        return linegraph.plot(dict_to_add, config)
