""" MultiQC module to parse output from HiCExplorer """


import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="HiCExplorer",
            anchor="hicexplorer",
            href="https://hicexplorer.readthedocs.io",
            info=" addresses the common tasks of Hi-C analysis from processing to visualization.",
            doi=["10.1093/nar/gky504", "10.1093/nar/gkaa220"],
        )

        self.hicexplorer_data = dict()
        for f in self.find_log_files("hicexplorer"):
            if f["fn"] != "QC_table.txt":
                # Parse the log file
                parsed_data = self.parse_logs(f["f"])
                # Build the sample ID
                s_name = f"{f['root']}_{f['s_name']}_{parsed_data['File']}"
                s_name = self.clean_s_name(s_name, f)
                # Save results
                if s_name in self.hicexplorer_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.hicexplorer_data[s_name] = parsed_data
                self.add_data_source(f, s_name=s_name)

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, s_name)

        self.hicexplorer_data = self.ignore_samples(self.hicexplorer_data)

        if len(self.hicexplorer_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.hicexplorer_data)} reports")

        self.write_data_file(self.hicexplorer_data, "multiqc_hicexplorer")

        self.colors = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
            "#D2691E",
        ]

        # detect version of QC file
        # no 'Pairs mappable, unique and high quality' --> version 1.7
        # Contains 'pairs used' --> 1.8 - 3.1
        # Contains 'Hi-C contacts' --> since 3.2
        hicexplorer_versions = set()
        for s_name in self.hicexplorer_data:
            # compatibility to HiCExplorer <= 1.7 version QC files
            if (
                "Pairs mappable, unique and high quality" not in self.hicexplorer_data[s_name]
                and "Pairs considered" in self.hicexplorer_data[s_name]
            ):
                hicexplorer_versions.add("1.7")
                self.hicexplorer_data[s_name]["Pairs mappable, unique and high quality"] = self.hicexplorer_data[
                    s_name
                ]["Pairs considered"]
                for key in ["One mate unmapped", "One mate not unique", "One mate low quality"]:
                    self.hicexplorer_data[s_name]["Pairs mappable, unique and high quality"] -= self.hicexplorer_data[
                        s_name
                    ][key]

            # compatibility to HiCExplorer <= 3.1 version QC files
            if "Pairs considered" in self.hicexplorer_data[s_name]:
                self.hicexplorer_data[s_name]["Sequenced reads"] = self.hicexplorer_data[s_name]["Pairs considered"]
                self.hicexplorer_data[s_name]["Hi-c contacts"] = self.hicexplorer_data[s_name]["Pairs used"]
                self.hicexplorer_data[s_name]["Low mapping quality"] = self.hicexplorer_data[s_name][
                    "One mate low quality"
                ]
                self.hicexplorer_data[s_name]["Intra short range (< 20kb)"] = self.hicexplorer_data[s_name][
                    "Short range"
                ]
                self.hicexplorer_data[s_name]["Intra long range (>= 20kb)"] = self.hicexplorer_data[s_name][
                    "Long range"
                ]
                self.hicexplorer_data[s_name]["Read pair type: inward pairs"] = self.hicexplorer_data[s_name][
                    "Inward pairs"
                ]
                self.hicexplorer_data[s_name]["Read pair type: outward pairs"] = self.hicexplorer_data[s_name][
                    "Outward pairs"
                ]
                self.hicexplorer_data[s_name]["Read pair type: left pairs"] = self.hicexplorer_data[s_name][
                    "Left pairs"
                ]
                self.hicexplorer_data[s_name]["Read pair type: right pairs"] = self.hicexplorer_data[s_name][
                    "Right pairs"
                ]

            elif "Sequenced reads" in self.hicexplorer_data[s_name]:
                hicexplorer_versions.add("3.2")

        log.debug(f"HiCExplorer versions: {', '.join(hicexplorer_versions)}")

        # prepare the basic statistics for hicexplorer
        self.hicexplorer_basic_statistics()

        # key lists for plotting
        keys_categorization_of_reads_considered = [
            "Pairs mappable, unique and high quality",
            "One mate unmapped",
            "One mate not unique",
            "Low mapping quality",
        ]
        keys_mappable_unique_and_high_quality = [
            "Hi-c contacts",
            "Self ligation (removed)",
            "Same fragment",
            "Self circle",
            "Dangling end",
            "One mate not close to rest site",
            "Duplicated pairs",
        ]
        keys_list_contact_distance = ["Intra short range (< 20kb)", "Intra long range (>= 20kb)", "Inter chromosomal"]
        keys_list_read_orientation = [
            "Read pair type: inward pairs",
            "Read pair type: outward pairs",
            "Read pair type: left pairs",
            "Read pair type: right pairs",
            "Inter chromosomal",
        ]

        # prepare the detail report section
        self.add_section(
            name="Mapping statistics",
            anchor="hicexplorer_categorization_of_considered_reads",
            plot=self.hicexplorer_create_plot(
                keys_categorization_of_reads_considered,
                "HiCExplorer: Categorization of considered reads",
                "categorization",
            ),
            description="This shows how the sequenced read pairs were mapped and those filtered due to mapping problems.",
            helptext="""
                * **Pairs mappable, unique and high quality**
                    * The count of reads that were considered as valid reads and were not one of the following:
                * **One mate unmapped**
                    * Filtered out read because one mate was not mapped.
                * **One mate not unique**
                    * Filtered out read because one mate was not unique.
                * **Low mapping quality**
                    * Filtered out because one mate was having a low quality.
            """,
        )

        self.add_section(
            name="Read filtering",
            anchor="hicexplorer_pairs_categorized",
            plot=self.hicexplorer_create_plot(
                keys_mappable_unique_and_high_quality,
                "HiCExplorer: Categorization of reads - Pairs mappable, unique and high quality",
                "mapping",
            ),
            description="This figure contains the number of reads that were finally used to build the "
            "Hi-C matrix along with the reads that where filtered out.",
            helptext="""
            * **Dangling ends**
                * These are reads that start with the restriction site and constitute reads that were digested but no ligated.
            * **Same fragment**
                * These are read mates, facing inward, separated by up to 800 bp that do not have a restriction enzyme in between.
                * These read pairs are not valid Hi-C pairs.
            * **Self circle**
                * Self circles are defined as pairs within 25kb with 'outward' read orientation
            * **Self ligation**
                * These are read pairs with a restriction site in between that are within 800 bp.""",
        )

        self.add_section(
            name="Contact distance",
            anchor="hicexplorer_contact_distance",
            plot=self.hicexplorer_create_plot(
                keys_list_contact_distance, "HiCExplorer: Contact distance", "contact_distance"
            ),
            description="This figure contains information about the distance and location of the valid pairs used.",
            helptext="""
            * **Intra long range**
                * Pairs with a distance greater than 20 kilobases
            * **Intra short range**
                * Pairs with a distance less than 20 kilobases
            * **Inter chromosomal**
                * Interchromosomal pairs.
            """,
        )

        self.add_section(
            name="Read orientation",
            anchor="hicexplorer_read_orientation",
            plot=self.hicexplorer_create_plot(
                keys_list_read_orientation, "HiCExplorer: Read orientation", "orientation"
            ),
            description="This figure contains information about the orientation of the read pairs.",
            helptext="""
                * **Inward pairs**
                    * First mate is a forward read, second is reverse.
                    * `--------------->              <----------------`
                * **Outward pairs**
                    * First mate is a reverse read, second is forward.
                    * `--------------->              <----------------`
                * **Left pairs**
                    * Both are reverse reads.
                    * `<---------------              <----------------`
                * **Right pairs**
                    * Both are forward reads.
                    * `--------------->              ---------------->`
            """,
        )

    def parse_logs(self, f):
        """Parse a given HiCExplorer log file from hicBuildMatrix."""
        data = {}
        for line in f.splitlines():
            # catch empty lines
            if len(line) == 0:
                continue
            s = line.split("\t")
            # Skip header
            if s[0].startswith("#"):
                continue
            data_ = []
            # catch lines with descriptive content: "Of Hi-C contacts:"
            for i in s[1:]:
                if len(i) == 0:
                    continue
                try:
                    i = i.replace("(", "")
                    i = i.replace(")", "")
                    i = i.replace(",", "")
                    data_.append(float(i))
                except ValueError:
                    data_.append(i)
            if len(data_) == 0:
                continue
            if s[0].startswith("Intra short range (< 20kb)"):
                s[0] = "Intra short range (< 20kb)"
            elif s[0].startswith("same fragment"):
                s[0] = "same fragment"
            elif s[0].startswith("short range"):
                s[0] = "short range"
            s[0] = s[0].capitalize()
            data[s[0]] = data_[0]
            try:
                data[s[0] + "_pct"] = data_[1]
            except IndexError:
                pass
        return data

    def hicexplorer_basic_statistics(self):
        """Create the general statistics for HiCExplorer."""
        data = {}
        for s_name in self.hicexplorer_data:
            max_distance_key = "Max rest. site distance"
            total_pairs = self.hicexplorer_data[s_name]["Sequenced reads"]
            try:
                self.hicexplorer_data[s_name][max_distance_key]
            except KeyError:
                max_distance_key = "Max library insert size"
            data_ = {
                "Sequenced reads": self.hicexplorer_data[s_name]["Sequenced reads"],
                "Hi-c contacts": self.hicexplorer_data[s_name]["Hi-c contacts"] / total_pairs,
                "Mapped": self.hicexplorer_data[s_name]["One mate unmapped"] / total_pairs,
                "Min rest. site distance": self.hicexplorer_data[s_name]["Min rest. site distance"],
                max_distance_key: self.hicexplorer_data[s_name][max_distance_key],
            }
            data[s_name] = data_
        headers = {
            "Sequenced reads": {
                "title": f"{config.read_count_prefix} Pairs",
                "description": f"Total number of read pairs ({config.read_count_desc})",
                "shared_key": "read_count",
            },
            "Hi-c contacts": {
                "title": "% Used pairs",
                "max": 100,
                "min": 0,
                "modify": lambda x: x * 100,
                "suffix": "%",
            },
            "Mapped": {
                "title": "% Mapped",
                "max": 100,
                "min": 0,
                "modify": lambda x: (1 - x) * 100,
                "scale": "RdYlGn",
                "suffix": "%",
            },
            "Min rest. site distance": {
                "title": "Min RE dist",
                "description": "Minimum restriction site distance (bp)",
                "format": "{:.0f}",
                "suffix": " bp",
            },
            max_distance_key: {
                "title": "Max RE dist",
                "description": f"{max_distance_key} (bp)",
                "format": "{:.0f}",
                "suffix": " bp",
            },
        }

        self.general_stats_addcols(data, headers)

    def hicexplorer_create_plot(self, pKeyList, pTitle, pId):
        """Create the graphics containing information about the read quality."""

        keys = dict()
        for i, key_ in enumerate(pKeyList):
            keys[key_] = {"color": self.colors[i]}

        config = {
            "id": "hicexplorer_" + pId + "_plot",
            "title": pTitle,
            "ylab": "Number of Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.hicexplorer_data, keys, config)
