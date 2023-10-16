""" MultiQC module to parse output from Lima """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="Lima",
            anchor="lima",
            href="https://github.com/PacificBiosciences/barcoding",
            info=" is used to demultiplex PacBio single-molecule sequencing reads.",
            # No publication / DOI // doi=
        )

        # To store the summary data
        self.lima_summary = dict()
        self.lima_counts = dict()

        # Parse the output files
        self.parse_summary_files()
        self.parse_counts_files()

        # Remove filtered samples
        self.lima_summary = self.ignore_samples(self.lima_summary)
        self.lima_counts = self.ignore_samples(self.lima_counts)

        # Write the data files to disk
        if not self.lima_summary and not self.lima_counts:
            raise ModuleNoSamplesFound
        log.info("Found {} reports".format(len(self.lima_summary)))
        log.info("Found {} samples".format(len(self.lima_counts)))

        if self.lima_summary:
            self.write_data_file(self.lima_summary, "multiqc_lima_summary")
        if self.lima_counts:
            self.write_data_file(self.lima_counts, "multiqc_lima_counts")

        # Add a graph of all filtered ZMWs
        self.plot_filter_data()

        # Here, we make a difference between samples that have been renamed,
        # and samples that have their name derived from lima itself. Since lima
        # uses barcode1--barcode2 as sample names, these will never match the
        # other samples in the MultiQC report.
        #
        # Therefore, we do two different things here.
        # 1. All samples that have been renamed will be added to the general
        #       statistics table.
        # 2. All samples that have not been renamed will be added to their own
        #       table in the Lima section of the report.
        lima_renamed_count = dict()
        lima_original_count = dict()

        for sample in self.lima_counts:
            if sample in config.sample_names_replace.values():
                lima_renamed_count[sample] = self.lima_counts[sample]
            else:
                lima_original_count[sample] = self.lima_counts[sample]

        # Get the headers and tconfig for the table
        headers, tconfig = self.make_headers_config()

        # Add a graph for the data values in the counts file, for the samples
        # that haven't been renamed
        if lima_original_count:
            self.make_counts_table(lima_original_count, headers, tconfig)

        # Add renamed samples to the general statistics table, since we assume
        # they are named consistenly now
        if lima_renamed_count:
            self.add_general_stats(lima_renamed_count, headers)

    def parse_summary_files(self):
        for f in self.find_log_files("lima/summary", filehandles=True):
            data = parse_PacBio_log(f["f"])
            if len(data) > 0:
                if f["s_name"] in self.lima_summary:
                    log.debug(f"Duplicate summary sample name found! Overwriting: {f['s_name']}")
                self.lima_summary[f["s_name"]] = data
                self.add_data_source(f)

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, f["s_name"])

    def parse_counts_files(self):
        for f in self.find_log_files("lima/counts", filehandles=True):
            data = self.parse_lima_counts(f["f"], f)
            # Check for duplicate samples
            for sample in data:
                if sample in self.lima_counts:
                    log.debug(f"Duplicate counts sample name found! Overwriting: {sample}")
            # After warning the user, overwrite the samples
            self.lima_counts.update(data)
            self.add_data_source(f)

    def parse_lima_counts(self, file_content, f):
        """Parse lima counts file"""

        # The first line is the header
        header = next(file_content).strip().split()

        # A dictionary to store the results
        lima_counts = dict()
        for line in file_content:
            spline = line.strip().split()
            data = {field: value for field, value in zip(header, spline)}

            first_barcode = data["IdxFirstNamed"]
            second_barcode = data["IdxCombinedNamed"]
            # The format barcode1--barcode2 is also used in
            sample = self.clean_s_name(f"{first_barcode}--{second_barcode}", f)
            counts = data["Counts"]
            mean_score = data["MeanScore"]
            lima_counts[sample] = {"Counts": int(counts), "MeanScore": float(mean_score)}

        return lima_counts

    def plot_filter_data(self):
        # First, we gather all filter results for each lima summary
        all_filters = set()
        for data in self.lima_summary.values():
            for reason in self.filter_and_pass(data):
                all_filters.add(reason)
        all_filters = list(all_filters)
        # Move passing reads to the start of the list
        all_filters.insert(0, all_filters.pop(all_filters.index("ZMWs above all thresholds")))

        # Flatten the data into a structure that can be plotted
        plot_data = dict()
        for s_name, data in self.lima_summary.items():
            plot_data[s_name] = dict()
            for reason in all_filters:
                plot_data[s_name][reason] = dict()
                # This is where the filter reasons are stored
                zmw_marginals = data["ZMWs below any threshold"]["ZMW marginals"]
                for filter_reason in zmw_marginals:
                    if filter_reason == reason:
                        plot_data[s_name][reason] = zmw_marginals[filter_reason]["count"]
                        break
                # If we didn't find filter reason for this filename, set it to zero
                else:
                    plot_data[s_name][reason] = 0
                # Unless 'reason' is actually the special case of the number of unfiltered reads
                if reason == "ZMWs above all thresholds":
                    plot_data[s_name][reason] = data["ZMWs above all thresholds"]["count"]

        # Plot configuration
        config = {
            "id": "lima-filter-graph",
            "title": "Lima: ZMW filtering results",
            "ylab": "Number of ZMWs",
            "xlab": "Lima summary file",
        }

        self.add_section(
            name="ZMW filtering",
            anchor="lima-filter",
            description=("The number of ZMWs that failed or passed all Lima filters"),
            helptext="""
                The number of ZMWs that passed all filters is shown as **ZMWs above all thresholds**,
                all other categories that are shown in the graph represent the number of ZMWs that
                were dropped for the specified reason.
            """,
            plot=bargraph.plot(plot_data, all_filters, config),
        )

    def filter_and_pass(self, data):
        """Get the reasons why each ZMW was filtered"""
        reasons = dict()

        # Add why ZMWs were filtered
        filter_data = data["ZMWs below any threshold"]["ZMW marginals"]
        for reason in filter_data:
            reasons[reason] = filter_data[reason]["count"]
        # Add the ZMWs that passed
        reasons["ZMWs above all thresholds"] = data["ZMWs above all thresholds"]["count"]

        return reasons

    def make_headers_config(self):
        """Prepare the headers and config for the lima counts table"""
        tconfig = {
            "id": "multiqc_lima_counts",
            "namespace": "Lima",
            "title": "Lima: Number of Reads",
            "anchor": "multiqc_lima_counts",
            "ylab": "# Reads",
        }

        headers = OrderedDict()
        headers["Counts"] = {
            "title": f"Read Count ({config.long_read_count_prefix})",
            "description": f"Number of reads for each sample or barcode pair ({config.long_read_count_desc})",
            "modify": lambda x: x * config.long_read_count_multiplier,
            "shared_key": "long_read_counts",
            "format": "{:,.2f}",
            "scale": "PuBuGn",
        }
        headers["MeanScore"] = {
            "title": "Quality Score",
            "description": "The mean quality score of the reads for each sample or barcode pair",
            "scale": "Spectral",
        }
        return headers, tconfig

    def make_counts_table(self, counts, headers, tconfig):
        self.add_section(
            name="Per sample count data",
            anchor="multiqc_lima_count",
            description="""
                Per sample or per barcode statistics from Lima.
                For instructions on how to display sample names instead of `barcode--barcode` pairs,
                please see the [MultiQC Lima documentation](https://multiqc.info/docs/#lima).
            """,
            plot=table.plot(counts, headers, tconfig),
        )

    def add_general_stats(self, counts, headers):
        """Add (renamed) samples to the general statistics table"""
        self.general_stats_addcols(counts, headers)


def parse_PacBio_log(file_content):
    """Parse PacBio log file"""
    data = dict()
    # This is a local dictionary to store which annotations belong to which
    # result dictionary. This will be used to structure the output, but will
    # not be part of the output itself
    annotations = dict()
    current_annotation = None

    for line in file_content:
        # Get rid of trailing newlines
        line = line.strip("\n")
        # Did we enter a new section with annotations for an earlier result?
        # If so, we will only add an empty dictionary with the correct name
        # These field are of the format "something something for (A):"
        section_header_pattern = " for [(][A-Z][)][:]$"
        if re.search(section_header_pattern, line):
            linedata = parse_line(line)
            ann = linedata["annotation"]
            # Cut off the ' for (B):' part
            name = line[:-9]
            # We make a new heading with the current name under the data that
            # matches the current annotation
            current_annotation = dict()
            # We keep the dictonary accessible under 'current_annotation',
            # so we can keep adding new data to it without having to keep track
            # of where it belongs
            annotations[ann][name] = current_annotation
            continue

        linedata = parse_line(line)

        # If we got no data, we reached the end of the annotated section
        if not linedata:
            current_annotation = None
            continue

        # Lets get the name of the data
        name = linedata.pop("name")
        # If we are in an annotated section, we add the data to the current
        # annotation
        if current_annotation is not None:
            current_annotation[name] = linedata
        # Otherwise, we add the newfound annotation to the dictionary in case
        # we find a corresponding section later on.
        # The data from the current line we add directly to the output data
        else:
            annotation = linedata.pop("annotation")
            annotations[annotation] = linedata
            data[name] = linedata

    return data


def parse_line(line):
    """Parse a line from the Lima log file"""
    data = dict()

    # If we got an empty line to parse
    if not line.strip():
        return data

    # Split the line on the colon character
    key, values = line.strip().split(":")

    # The key can have multiple parts
    keys = key.strip().split()

    # Does the key have an annotation (A), (B) etc
    if re.fullmatch("[(][A-Z][)]", keys[-1]):
        # We store the annotation without the bracets
        data["annotation"] = keys[-1][1:-1]
        # And we add the rest of the key as name
        data["name"] = " ".join(keys[:-1])
    # Otherwise, we just store the name
    else:
        data["name"] = " ".join(keys)

    # Parsing the values
    values = values.strip().split()
    # Are there are no values we are done
    if not values:
        return data

    # If there is a single value
    if len(values) == 1:
        value = values.pop()
        # Is the value a percentage
        if value.endswith("%"):
            data["percentage"] = float(value[:-1])
        # Otherwise, it is an integer
        else:
            data["count"] = int(value)
    elif len(values) == 2:
        # If there are multiple values, they are in the format
        # count (percentage%)
        count = values[0]
        percentage = values[1]

        # The percentage should be in the format:
        # (12.34%) or (100%) or (0%) or (-nan%)
        # So we can remove the bracets first
        percentage = percentage[1:-1]
        # Then we make sure it is one of the supported patterns
        assert re.fullmatch(r"(\d+\.\d+%|\d+%|-nan%)", percentage)

        # If the percentage is this weird nan, set it to 0
        if percentage == "-nan%":
            data["percentage"] = 0.0
        else:  # Otherwise, we cut of the % sign and convert to float
            data["percentage"] = float(percentage[:-1])
        # Add the count as an integer
        data["count"] = int(count)

    return data
