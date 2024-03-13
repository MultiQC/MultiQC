""" MultiQC module to parse output from FastQ Screen """


import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"Fastq_screen version: ([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FastQ Screen",
            anchor="fastq_screen",
            href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/",
            info="allows you to screen a library of sequences in FastQ format against"
            " a set of sequence databases so you can see if the composition of the"
            " library matches with what you expect.",
            doi="10.12688/f1000research.15931.2",
        )

        # Find and load any FastQ Screen reports
        self.fq_screen_data = dict()
        self.num_orgs = 0
        for f in self.find_log_files("fastq_screen", filehandles=True):
            parsed_data = self.parse_fqscreen(f)
            if parsed_data is not None:
                if f["s_name"] in self.fq_screen_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f)
                self.fq_screen_data[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.fq_screen_data = self.ignore_samples(self.fq_screen_data)

        if len(self.fq_screen_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.fq_screen_data)} reports")

        # Section 1 - Alignment Profiles
        self.add_section(name="Mapped Reads", anchor="fastq_screen_mapped_reads", plot=self.fqscreen_simple_plot())

        # Section 2 - Optional bisfulfite plot
        self.fqscreen_bisulfite_plot()

        # Write the total counts and percentages to files
        self.write_data_file(self.parse_csv(), "multiqc_fastq_screen")

    def parse_fqscreen(self, f):
        """Parse the FastQ Screen output into a 3D dict"""
        parsed_data = {}
        nohits_pct = None
        headers = None
        bs_headers = None
        for line in f["f"]:
            # Skip comment lines
            if line.startswith("#"):
                version_match = re.search(VERSION_REGEX, line)
                if version_match:
                    self.add_software_version(version_match.group(1), f["s_name"])
                continue
            if line.startswith("%Hit_no_genomes:") or line.startswith("%Hit_no_libraries:"):
                nohits_pct = float(line.split(":", 1)[1])
                parsed_data["No hits"] = {"percentages": {"one_hit_one_library": nohits_pct}}
            else:
                s = line.strip().split("\t")

                # Regular FastQ Screen table section
                if len(s) == 12:
                    if headers is None:
                        headers = s
                    else:
                        # Can't use #Reads in subset as varies. #Reads_processed should be same for all orgs in a sample
                        parsed_data["total_reads"] = int(s[1])
                        # Loop through all columns
                        parsed_data[s[0]] = {"percentages": {}, "counts": {}}
                        for idx, h in enumerate(headers):
                            if idx == 0:
                                continue
                            dtype = "percentages" if h.startswith("%") else "counts"
                            field = (
                                h.replace("%", "")
                                .replace("#", "")
                                .replace("genomes", "libraries")
                                .replace("genome", "library")
                                .lower()
                            )
                            parsed_data[s[0]][dtype][field] = float(s[idx])

                # Optional Bisulfite table section
                elif len(s) == 9:
                    if bs_headers is None:
                        bs_headers = s
                    else:
                        # Loop through all columns
                        parsed_data[s[0]]["bisulfite_percentages"] = {}
                        parsed_data[s[0]]["bisulfite_counts"] = {}
                        for idx, h in enumerate(bs_headers):
                            if idx == 0:
                                continue
                            dtype = "bisulfite_percentages" if h.startswith("%") else "bisulfite_counts"
                            field = h.replace("%", "").replace("#", "").lower()
                            parsed_data[s[0]][dtype][field] = float(s[idx])

        if len(parsed_data) == 0:
            return None

        # Calculate no hits counts
        if parsed_data.get("total_reads") is not None and nohits_pct is not None:
            parsed_data["No hits"]["counts"] = {
                "one_hit_one_library": int((nohits_pct / 100.0) * float(parsed_data["total_reads"]))
            }
        else:
            log.warning(f"Couldn't find number of reads with no hits for '{f['s_name']}'")

        self.num_orgs = max(len(parsed_data), self.num_orgs)
        return parsed_data

    def parse_csv(self):
        totals = dict()
        for s in sorted(self.fq_screen_data.keys()):
            totals[s] = dict()
            for org in self.fq_screen_data[s]:
                if org == "total_reads":
                    totals[s]["total_reads"] = self.fq_screen_data[s][org]
                    continue
                try:
                    k = f"{org} counts"
                    totals[s][k] = self.fq_screen_data[s][org]["counts"]["one_hit_one_library"]
                    totals[s][k] += self.fq_screen_data[s][org]["counts"].get("multiple_hits_one_library", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["counts"].get("one_hit_multiple_libraries", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["counts"].get("multiple_hits_multiple_libraries", 0)
                except KeyError:
                    pass
                try:
                    k = f"{org} percentage"
                    totals[s][k] = self.fq_screen_data[s][org]["percentages"]["one_hit_one_library"]
                    totals[s][k] += self.fq_screen_data[s][org]["percentages"].get("multiple_hits_one_library", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["percentages"].get("one_hit_multiple_libraries", 0)
                    totals[s][k] += self.fq_screen_data[s][org]["percentages"].get(
                        "multiple_hits_multiple_libraries", 0
                    )
                except KeyError:
                    pass
        return totals

    def fqscreen_simple_plot(self):
        """Makes a simple bar plot with summed alignment counts for
        each species, stacked."""

        # First, sum the different types of alignment counts
        data = {}
        org_counts = {}
        for s_name in sorted(self.fq_screen_data):
            data[s_name] = dict()
            sum_alignments = 0
            for org in self.fq_screen_data[s_name]:
                if org == "total_reads":
                    continue
                try:
                    data[s_name][org] = self.fq_screen_data[s_name][org]["counts"]["one_hit_one_library"]
                except KeyError:
                    log.error(
                        "No counts found for '{}' ('{}'). Could be malformed or very old FastQ Screen results. Skipping sample".format(
                            org, s_name
                        )
                    )
                    continue
                try:
                    data[s_name][org] += self.fq_screen_data[s_name][org]["counts"]["multiple_hits_one_library"]
                except KeyError:
                    pass
                sum_alignments += data[s_name][org]
                org_counts[org] = org_counts.get(org, 0) + data[s_name][org]

            # Calculate hits in multiple genomes
            if "total_reads" in self.fq_screen_data[s_name]:
                data[s_name]["Multiple Genomes"] = self.fq_screen_data[s_name]["total_reads"] - sum_alignments

        # Sort the categories by the total read counts
        cats = dict()
        for org in sorted(org_counts, key=org_counts.get, reverse=True):
            if org not in cats and org != "No hits":
                cats[org] = {"name": org}

        # Strip empty dicts
        for s_name in list(data.keys()):
            if len(data[s_name]) == 0:
                del data[s_name]

        pconfig = {
            "id": "fastq_screen_plot",
            "title": "FastQ Screen: Mapped Reads",
            "cpswitch_c_active": False,
            "ylab": "Mapped reads",
        }
        cats["Multiple Genomes"] = {"name": "Multiple Genomes", "color": "#820000"}
        cats["No hits"] = {"name": "No hits", "color": "#cccccc"}

        return bargraph.plot(data, cats, pconfig)

    def fqscreen_bisulfite_plot(self):
        """Make a stacked barplot for the bisulfite data, if we have any"""

        pconfig = {
            "id": "fastq_screen_bisulfite_plot",
            "title": "FastQ Screen: Bisulfite Mapping Strand Orientation",
            "hide_zero_cats": False,
            "ylab": "Reads",
            "data_labels": [],
        }

        cats = {
            "original_top_strand": {"name": "Original top strand", "color": "#80cdc1"},
            "complementary_to_original_top_strand": {
                "name": "Complementary to original top strand",
                "color": "#018571",
            },
            "complementary_to_original_bottom_strand": {
                "name": "Complementary to original bottom strand",
                "color": "#a6611a",
            },
            "original_bottom_strand": {"name": "Original bottom strand", "color": "#dfc27d"},
        }

        # Pull out the data that we want
        pdata_unsorted = {}
        org_counts = {}
        max_count = 0
        for s_name, d in self.fq_screen_data.items():
            for org, dd in d.items():
                if org not in pdata_unsorted:
                    pdata_unsorted[org] = {}
                    org_counts[org] = 0
                try:
                    pdata_unsorted[org][s_name] = dd["bisulfite_counts"]
                    total_counts = sum([c for c in dd["bisulfite_counts"].values()])
                    org_counts[org] += total_counts
                    max_count = max(max_count, total_counts)
                except (TypeError, KeyError):
                    pass

        # Consistent y-max across all genomes
        pconfig["ymax"] = max_count

        # Show tabbed bar plots, sorted by total read count
        pdata = []
        pcats = []
        for org in sorted(org_counts, key=org_counts.get, reverse=True):
            if org_counts[org] > 0:
                pconfig["data_labels"].append(org)
                pdata.append(pdata_unsorted[org])
                pcats.append(cats)

        if len(pdata) == 0:
            return None

        self.add_section(
            name="Bisulfite Reads", anchor="fastq_screen_bisulfite", plot=bargraph.plot(pdata, pcats, pconfig)
        )
