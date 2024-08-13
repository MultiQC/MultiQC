import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses the pychopper stats file. Pychopper needs to be run with the `-S stats_output` option to create the file. The name of the output file defines the sample name.

    The stats file is a three column `tsv` file with the format `category name value`.

    Currently only two stats are displayed in MultiQC. Two bargraphs are created for the read classication and the strand orientation of the identified full length transcripts. Additional stats could be included on further request.

    The general stats table contains a value that displays the percentage of full length transcripts. This value is calculated from the cumulative length of reads where Pychopper found primers at both ends.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Pychopper",
            anchor="pychopper",
            href="https://github.com/nanoporetech/pychopper",
            info="Identifies, orients, trims and rescues full length Nanopore cDNA reads. Can also rescue fused reads.",
            # Can't find a DOI // doi=
        )

        # Parse stats file
        self.pychopper_data = {}
        for f in self.find_log_files("pychopper"):
            sample = f["s_name"]
            self.pychopper_data[sample] = {}
            lines = f["f"].splitlines()
            for line in lines[1:]:
                category, name, value = line.split()
                if category not in self.pychopper_data[sample]:
                    self.pychopper_data[sample][category] = {}
                self.pychopper_data[sample][category][name] = float(value)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.pychopper_data = self.ignore_samples(self.pychopper_data)

        # Raise user warning if no data found
        if len(self.pychopper_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.pychopper_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add to general statistics table:
        # Percentage of full length transcripts
        data_general_stats = {}
        for sample in self.pychopper_data.keys():
            data_general_stats[sample] = {}
            c = self.pychopper_data[sample]["Classification"]
            ftp = c["Primers_found"] * 100 / (c["Primers_found"] + c["Rescue"] + c["Unusable"])
            data_general_stats[sample]["ftp"] = ftp

        headers = {
            "ftp": {
                "title": "Full-Length cDNA",
                "description": "Percentage of full length cDNA reads with correct primers at both ends",
                "suffix": "%",
                "max": 100,
                "min": 0,
            }
        }

        self.general_stats_addcols(data_general_stats, headers)

        # Write data file
        self.write_data_file(self.pychopper_data, "multiqc_pychopper")

        # Report sections
        self.add_section(
            name="cDNA Read Classification",
            description=(
                """
                This plot shows the cDNA read categories identified by Pychopper </br>
                """
            ),
            helptext=(
                """
                There are three possible cases:

                * **Primers found**: Full length cDNA reads with correct primers at both ends.
                * **Rescued reads**: Split fusion reads.
                * **Unusable**: Reads without correct primer combinations.
                """
            ),
            anchor="pychopper_classification",
            plot=self.plot_classification(),
        )
        self.add_section(
            name="cDNA Strand Orientation",
            description=(
                """
                This plot shows the strand orientation of full length cDNA reads
                """
            ),
            helptext=(
                """
                Nanopore cDNA reads are always read forward. To estimate their original strand,
                Pychopper searches for the location of the start and end primers and assigns the reads accordingly.
                """
            ),
            anchor="pychopper_orientation",
            plot=self.plot_orientation(),
        )

    # Plotting functions
    def plot_classification(self):
        """Generate the cDNA read classification plot"""

        pconfig = {
            "id": "pychopper_classification_plot",
            "title": "Pychopper: Read classification",
            "ylab": "",
            "x_decimals": False,
            "ymin": 0,
        }

        data_classification = {}
        for sample in self.pychopper_data.keys():
            data_classification[sample] = {}
            data_classification[sample] = self.pychopper_data[sample]["Classification"]

        cats = ["Primers_found", "Rescue", "Unusable"]
        return bargraph.plot(data_classification, cats, pconfig)

    def plot_orientation(self):
        """Generate the read strand orientation plot"""

        pconfig = {
            "id": "pychopper_orientation_plot",
            "title": "Pychopper: Strand Orientation",
            "ylab": "",
            "cpswitch_c_active": False,
            "x_decimals": False,
            "ymin": 0,
        }

        data_orientation = {}
        for sample in self.pychopper_data.keys():
            data_orientation[sample] = {}
            data_orientation[sample] = self.pychopper_data[sample]["Strand"]

        cats = ["+", "-"]
        return bargraph.plot(data_orientation, cats, pconfig)
