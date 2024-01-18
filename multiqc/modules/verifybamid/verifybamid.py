""" MultiQC module to parse output from VerifyBAMID """


import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    module class, parses stderr logs.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="VerifyBAMID",
            anchor="verifybamid",
            href="https://genome.sph.umich.edu/wiki/VerifyBamID",
            info="detects sample contamination and/or sample swaps.",
            doi="10.1016/j.ajhg.2012.09.004",
        )

        # flag to hide columns if no chip data
        self.hide_chip_columns = True

        # default values for columns
        self.col_config_defaults = {
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.3f}",
            "modify": lambda x: x * 100.0 if x != "NA" else x,
            "scale": "OrRd",
        }

        # dictionary to hold all data for each sample
        self.verifybamid_data = dict()

        # for each file ending in self.SM
        for f in self.find_log_files("verifybamid/selfsm"):
            # pass the file to function self.parse_selfsm to parse file
            parsed_data = self.parse_selfsm(f)
            # if a result was returned
            if parsed_data is not None:
                # for each sample extracted from the file
                for s_name in parsed_data:
                    # if there are duplicate sample names
                    if s_name in self.verifybamid_data:
                        # write this to log
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    # add the sample as a key to the verifybamid_data dictionary and the dictionary of values as the value
                    self.verifybamid_data[s_name] = parsed_data[s_name]

                    # add data source to multiqc_sources.txt
                    self.add_data_source(f, s_name)

        # Filter to strip out ignored sample names as per config.yaml
        self.verifybamid_data = self.ignore_samples(self.verifybamid_data)

        if len(self.verifybamid_data) == 0:
            raise ModuleNoSamplesFound

        # print number of verifyBAMID reports found and parsed
        log.info(f"Found {len(self.verifybamid_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.verifybamid_data, "multiqc_verifybamid")

        # add to General Stats Table
        self.verifybamid_general_stats_table()

        # add section with the values from the verify BAM ID output
        self.verifybamid_table()

    def parse_selfsm(self, f):
        """Go through selfSM file and create a dictionary with the sample name as a key,"""
        # create a dictionary to populate from this sample's file
        parsed_data = dict()
        # set a empty variable which denotes if the headers have been read
        headers = None
        # for each line in the file
        for line in f["f"].splitlines():
            # split the line on tab
            s = line.split("\t")
            # if we haven't already read the header line
            if headers is None:
                # assign this list to headers variable
                headers = s
            # for all rows after the first
            else:
                # clean the sample name (first column) and assign to s_name
                s_name = self.clean_s_name(s[0], f)
                # create a dictionary entry with the first column as a key (sample name) and empty dictionary as a value
                parsed_data[s_name] = {}
                # for each item in list of items in the row
                for i, v in enumerate(s):
                    # if it's not the first element (if it's not the name)
                    if i != 0:
                        # see if CHIP is in the column header and the value is not NA
                        if headers[i].startswith("CHIP") and v != "NA":
                            # set hide_chip_columns = False so they are not hidden
                            self.hide_chip_columns = False
                        # try and convert the value into a float
                        try:
                            # and add to the dictionary the key as the corrsponding item from the header and the value from the list
                            parsed_data[s_name][headers[i]] = float(v)
                        # if can't convert to float...
                        except ValueError:
                            # add to the dictionary the key as the corrsponding item from the header and the value from the list
                            parsed_data[s_name][headers[i]] = v

                # Rename some column titles that can be variably named (See issue #1112)
                if "FREEMIX(alpha)" in parsed_data[s_name]:
                    parsed_data[s_name]["FREEMIX"] = parsed_data[s_name]["FREEMIX(alpha)"]

        # else return the dictionary
        return parsed_data

    def verifybamid_general_stats_table(self):
        """Take the percentage of contamination from all the parsed *.SELFSM files and add it to the basic stats table at the top of the report"""

        # create a dictionary to hold the columns to add to the general stats table
        headers = dict()
        # available columns are:
        # SEQ_ID RG  CHIP_ID #SNPS   #READS  AVG_DP  FREEMIX FREELK1 FREELK0 FREE_RH FREE_RA CHIPMIX CHIPLK1 CHIPLK0 CHIP_RH CHIP_RA DPREF   RDPHET  RDPALT
        # see https://genome.sph.umich.edu/wiki/VerifyBamID#Interpreting_output_files

        # add the CHIPMIX column if we have data
        if not self.hide_chip_columns:
            headers["CHIPMIX"] = dict(
                self.col_config_defaults,
                **{
                    "title": "Contamination (S+A)",
                    "description": "VerifyBamID: CHIPMIX -   Sequence+array estimate of contamination (NA if the external genotype is unavailable)",
                },
            )

        # add the FREEMIX column. set the title and description
        headers["FREEMIX"] = dict(
            self.col_config_defaults,
            **{
                "title": "Contamination (S)",
                "description": "VerifyBamID: FREEMIX -   Sequence-only estimate of contamination.",
            },
        )

        # pass the data dictionary and header dictionary to function to add to table.
        self.general_stats_addcols(self.verifybamid_data, headers)

    def verifybamid_table(self):
        """
        Create a table with all the columns from verify BAM ID
        """

        # create an ordered dictionary to preserve the order of columns
        headers = dict()
        # add each column and the title and description (taken from verifyBAMID website)
        headers["RG"] = {
            "title": "Read Group",
            "description": "ReadGroup ID of sequenced lane.",
            "hidden": all([s["RG"] == "ALL" for s in self.verifybamid_data.values()]),
        }
        if not self.hide_chip_columns:
            headers["CHIP_ID"] = {"title": "Chip ID", "description": "ReadGroup ID of sequenced lane."}
        headers["#SNPS"] = {
            "title": "SNPS",
            "description": "# SNPs passing the criteria from the VCF file",
            "format": "{:,.0f}",
            "min": 0,
            "scale": "BuPu",
        }
        headers["#READS"] = {
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Number of reads loaded from the BAM file ({config.read_count_desc})",
            "format": "{:,.1f}",
            "modify": lambda x: x * config.read_count_multiplier if x != "NA" else x,
            "shared_key": "read_count",
            "min": 0,
            "scale": "GnBu",
        }
        headers["AVG_DP"] = {
            "title": "Average Depth",
            "description": "Average sequencing depth at the sites in the VCF file",
            "suffix": " X",
            "min": 0,
            "scale": "YlGn",
        }
        # use default columns
        headers["FREEMIX"] = dict(
            self.col_config_defaults,
            **{
                "title": "Contamination (Seq)",
                "description": "VerifyBamID: FREEMIX -   Sequence-only estimate of contamination.",
            },
        )
        headers["FREELK1"] = {
            "title": "FREEELK1",
            "format": "{:,.0f}",
            "description": "Maximum log-likelihood of the sequence reads given estimated contamination under sequence-only method",
            "min": 0,
            "scale": "RdYlGn",
        }
        headers["FREELK0"] = {
            "title": "FREELK0",
            "format": "{:,.0f}",
            "description": "Log-likelihood of the sequence reads given no contamination under sequence-only method",
            "min": 0,
            "scale": "RdYlGn",
        }
        headers["FREE_RH"] = {
            "title": "FREE_RH",
            "description": "Estimated reference bias parameter Pr(refBase|HET) (when --free-refBias or --free-full is used)",
            "hidden": all([s["FREE_RH"] == "NA" for s in self.verifybamid_data.values()]),
        }
        headers["FREE_RA"] = {
            "title": "FREE_RA",
            "description": "Estimated reference bias parameter Pr(refBase|HOMALT) (when --free-refBias or --free-full is used)",
            "hidden": all([s["FREE_RA"] == "NA" for s in self.verifybamid_data.values()]),
        }

        # Only print Chip columns to the report if we have data
        if not self.hide_chip_columns:
            headers["CHIPMIX"] = dict(
                self.col_config_defaults,
                **{
                    "title": "Contamination S+A",
                    "description": "VerifyBamID: CHIPMIX -   Sequence+array estimate of contamination (NA if the external genotype is unavailable)",
                },
            )
            headers["CHIPLK1"] = {
                "title": "CHIPLK1",
                "description": "Maximum log-likelihood of the sequence reads given estimated contamination under sequence+array method (NA if the external genotypes are unavailable)",
            }
            headers["CHIPLK0"] = {
                "title": "CHIPLK0",
                "description": " Log-likelihood of the sequence reads given no contamination under sequence+array method (NA if the external genotypes are unavailable)",
            }
            headers["CHIP_RH"] = {
                "title": "CHIP_RH",
                "description": "Estimated reference bias parameter Pr(refBase|HET) (when --chip-refBias or --chip-full is used)",
            }
            headers["CHIP_RA"] = {
                "title": "CHIP_RA",
                "description": "Estimated reference bias parameter Pr(refBase|HOMALT) (when --chip-refBias or --chip-full is used)",
            }

        headers["DPREF"] = {
            "title": "DPREF",
            "description": "Depth (Coverage) of HomRef site (based on the genotypes of (SELF_SM/BEST_SM), passing mapQ, baseQual, maxDepth thresholds.",
            "hidden": all([s["DPREF"] == "NA" for s in self.verifybamid_data.values()]),
        }
        headers["RDPHET"] = {
            "title": "RDPHET",
            "description": "DPHET/DPREF, Relative depth to HomRef site at Heterozygous site.",
            "hidden": all([s["RDPHET"] == "NA" for s in self.verifybamid_data.values()]),
        }
        headers["RDPALT"] = {
            "title": "RDPALT",
            "description": "DPHET/DPREF, Relative depth to HomRef site at HomAlt site.",
            "hidden": all([s["RDPALT"] == "NA" for s in self.verifybamid_data.values()]),
        }

        tconfig = {
            "namespace": "VerifyBAMID",
            "id": "verifybamid-results",
            "title": "VerifyBAMID: Results",
        }

        # send the plot to add section function with data dict and headers
        self.add_section(
            anchor="verifybamid-table",
            description="The following values provide estimates of sample contamination. Click help for more information.",
            helptext="""
            **Please note that `FREEMIX` is named _Contamination (Seq)_ and `CHIPMIX`
            is named _Contamination (S+A)_ in this MultiQC report.**

            VerifyBamID provides a series of information that is informative to determine
            whether the sample is possibly contaminated or swapped, but there is no single
            criteria that works for every circumstances. There are a few unmodeled factor
            in the estimation of `[SELF-IBD]/[BEST-IBD]` and `[%MIX]`, so please note that the
            MLE estimation may not always exactly match to the true amount of contamination.
            Here we provide a guideline to flag potentially contaminated/swapped samples:

            * Each sample or lane can be checked in this way.
              When `[CHIPMIX] >> 0.02` and/or `[FREEMIX] >> 0.02`, meaning 2% or more of
              non-reference bases are observed in reference sites, we recommend to examine
              the data more carefully for the possibility of contamination.
            * We recommend to check each lane for the possibility of sample swaps.
              When `[CHIPMIX] ~ 1` AND `[FREEMIX] ~ 0`, then it is possible that the sample
              is swapped with another sample. When `[CHIPMIX] ~ 0` in `.bestSM` file,
              `[CHIP_ID]` might be actually the swapped sample. Otherwise, the swapped
              sample may not exist in the genotype data you have compared.
            * When genotype data is not available but allele-frequency-based estimates of
              `[FREEMIX] >= 0.03` and `[FREELK1]-[FREELK0]` is large, then it is possible
              that the sample is contaminated with other sample. We recommend to use
              per-sample data rather than per-lane data for checking this for low coverage
              data, because the inference will be more confident when there are large number
              of bases with depth 2 or higher.

            _Copied from the [VerifyBAMID documentation](https://genome.sph.umich.edu/wiki/VerifyBamID) - see the link for more details._
            """,
            plot=table.plot(self.verifybamid_data, headers, tconfig),
        )
