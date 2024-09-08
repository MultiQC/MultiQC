import logging
import os

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module parses the base `misincorporation` output.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="mapDamage",
            anchor="mapdamage",
            href="https://github.com/ginolhac/mapDamage",
            info="Tracks and quantifies damage patterns in ancient DNA sequences.",
            doi="https://doi.org/10.1093/bioinformatics/btt193",
        )

        # Init empty dictionaries
        self.threepGtoAfreq_data = dict()
        self.fivepCtoTfreq_data = dict()
        self.lgdist_fw_data = dict()
        self.lgdist_rv_data = dict()

        # Find and load log files
        for f in self.find_log_files("mapdamage", filehandles=True):
            self.parse_logs(f)

        # Filter to strip out ignored sample names
        self.threepGtoAfreq_data = self.ignore_samples(self.threepGtoAfreq_data)
        self.fivepCtoTfreq_data = self.ignore_samples(self.fivepCtoTfreq_data)
        self.lgdist_fw_data = self.ignore_samples(self.lgdist_fw_data)
        self.lgdist_rv_data = self.ignore_samples(self.lgdist_rv_data)

        ## Stop computations if there is no data after ignoring samples
        if (
            len(self.threepGtoAfreq_data) == 0
            and len(self.fivepCtoTfreq_data) == 0
            and len(self.lgdist_fw_data) == 0
            and len(self.lgdist_rv_data) == 0
        ):
            raise ModuleNoSamplesFound

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        ## No need to run self.write_data_file as all the data imported is already in a TSV format and can be (almost) directly used for plotting.

        # Basic Stats Table, use generic function to add data to general table
        self.dmgprof_misinc_stats(self.threepGtoAfreq_data, "3 Prime", "G>A")
        self.dmgprof_misinc_stats(self.fivepCtoTfreq_data, "5 Prime", "C>T")

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
                description="Read length distribution for reverse strand (−) reads.",
                helptext="""
                This plot shows the read length distribution of the reverse reads in the investigated sample. Reads below lengths of 30bp are typically filtered, so the plot doesn't show these in many cases. A shifted distribution of read lengths towards smaller read lengths (e.g around 30-50bp) is also an indicator of ancient DNA.
                """,
                plot=self.lgdistplot(self.lgdist_rv_data, "Reverse"),
            )

    ## Parse misincorporation file into a dict with 1 key ('dmg_5p' or 'dmg_3p') and 1 value (list of 2nd column values)
    def parse_misincorporation(self, f):
        """Parse the misincorporation data from mapDamage output files"""
        misincorporation_dict = {}
        values = []
        try:
            for line in f["f"]:
                if line.startswith("#"):
                    continue
                else:
                    line = line.rstrip().split("\t")
                    if line[0] == "pos":
                        readend = line[1][0:2]
                    else:
                        values.append(float(line[1]))
            misincorporation_dict["dmg_" + readend] = values
            return misincorporation_dict
        except Exception as e:
            print(e)
            log.warning(f"Could not parse mapDamage misincorporation file: '{f['fn']}'")
            return None

    ## Parse length distribution file into a dict with 2 keys ('lendist_fw' and 'lendist_rv') and 1 value each (dict of {length : count})
    def parse_length_distribution(self, f):
        """Parse the length distribution data from mapDamage output files"""
        length_distribution_dict = {"lendist_fw": {}, "lendist_rv": {}}
        try:
            for line in f["f"]:
                ## Skip commented lines and header
                if line.startswith("#") or line.startswith("Std"):
                    continue
                else:
                    ## Lines are tab-separated, first column is strand, second column is read length, third column is count
                    line = line.rstrip().split("\t")
                    if line[0] == "+":
                        strand = "fw"
                    else:
                        strand = "rv"

                    length_distribution_dict["lendist_" + strand][line[1]] = int(line[2])
            return length_distribution_dict
        except Exception as e:
            print(e)
            log.warning(f"Could not parse mapDamage length distribution file: '{f['fn']}'")
            return None

    # Parse input files
    def parse_logs(self, f):
        """Parse the output from mapDamage"""

        # Get sample name from result directory name
        ## Remove the *results_ prefix from the sample name if it is there.
        s_name = self.clean_s_name(f["root"], f, root=os.path.dirname(f["root"])).lstrip("*results_")
        self.add_data_source(f, s_name)

        if f["fn"].endswith("_freq.txt") and f["fn"].startswith("3p"):
            # Add 3' G to A data
            self.threepGtoAfreq_data[s_name] = self.parse_misincorporation(f)["dmg_3p"]

        elif f["fn"].endswith("_freq.txt") and f["fn"].startswith("5p"):
            # Add 5' C to T data
            self.fivepCtoTfreq_data[s_name] = self.parse_misincorporation(f)["dmg_5p"]

        elif f["fn"] == "lgdistribution.txt":
            # Add lendist forward & reverse
            lgdist_data = self.parse_length_distribution(f)
            self.lgdist_fw_data[s_name] = lgdist_data["lendist_fw"]
            self.lgdist_rv_data[s_name] = lgdist_data["lendist_rv"]

    #### Tables from here on

    def dmgprof_misinc_stats(self, dict_to_plot, readend, substitution):
        """Take the parsed stats from the mapDamage and add it to the
        basic stats table at the top of the report"""

        headers = {
            f"mapdamage-{readend}1": {
                "rid": f"misinc-stats-1st-{readend}-{substitution}",
                "title": f"{readend} {substitution} 1st base",
                "description": f"{readend} 1st base substitution frequency for {substitution}",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGnBu",
                "modify": lambda x: x * 100.0,
            },
            f"mapdamage-{readend}2": {
                "rid": f"misinc-stats-2nd-{readend}-{substitution}",
                "title": f"{readend} {substitution} 2nd base",
                "description": f"{readend} 2nd base substitution frequency for {substitution}",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "BuGn",
                "hidden": True,
                "modify": lambda x: x * 100.0,
            },
        }

        # Create new small subset dictionary for entries (we need just the first two data (k,v) pairs from each report)
        # Only the first two parts are informative from both 3' and 5' ends of reads (1st, 2nd base damage pattern)
        dict_to_add = dict()

        for key in dict_to_plot.keys():
            tmp = dict_to_plot[key]
            pos = ["mapdamage-" + readend + "1", "mapdamage-" + readend + "2"]
            strlist = tmp[:2]
            tuples = list(zip(pos, strlist))
            data = dict((x, y) for x, y in tuples)
            # Extract first two elements from list
            dict_to_add[key] = data

        self.general_stats_addcols(dict_to_add, headers)

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
            "id": f"mapdamage-length-distribution-{orientation}",
            "title": f"mapDamage: Read length distribution - {orientation} ",
            "ylab": "Number of reads",
            "xlab": "Readlength (bp)",
            "x_decimals": False,
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
            "id": "mapdamage-threeprime_misinc_plot",
            "title": "mapDamage: 3' G>A misincorporation plot",
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
            "id": "mapdamage-fiveprime_misinc_plot",
            "title": "mapDamage: 5' C>T misincorporation plot",
            "ylab": "% C to T substituted",
            "xlab": "Nucleotide position from 5'",
            "tt_label": "{point.y:.2f} % C>T misincorporations at nucleotide position {point.x}",
            "ymin": 0,
            "xmin": 1,
        }

        return linegraph.plot(dict_to_add, config)
