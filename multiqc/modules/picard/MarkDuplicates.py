"""MultiQC submodule to parse output from Picard MarkDuplicates"""

import logging
from typing import Dict, List

import math
from collections import defaultdict

from multiqc import config
from multiqc.modules.picard import util
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(
    module,
    sp_key="picard/markdups",
):
    """
    Find Picard MarkDuplicates reports and parse their data.
    Note that this function is also used by the biobambam2 module, that's why
    the parameter.
    """

    data_by_sample: Dict = dict()

    # Get custom config value
    merge_multiple_libraries: bool = getattr(config, "picard_config", {}).get("markdups_merge_multiple_libraries", True)

    # Function to save results at end of table
    def save_table_results(s_name, keys, parsed_data, recompute_merged_metrics):
        # No data
        if len(keys) == 0 or len(parsed_data) == 0:
            return

        # User has requested each library is kept separate
        # Update the sample name to append the library name
        if not merge_multiple_libraries and len(parsed_data) > 0:
            s_name = f"{s_name} - {parsed_data['LIBRARY']}"

        # Skip - No reads
        try:
            if parsed_data["READ_PAIRS_EXAMINED"] == 0 and parsed_data["UNPAIRED_READS_EXAMINED"] == 0:
                log.warning(f"Skipping MarkDuplicates sample '{s_name}' as log contained no reads")
                return
        # Skip - Missing essential fields
        except KeyError:
            log.warning(f"Skipping MarkDuplicates sample '{s_name}' as missing essential fields")
            return

        # Recompute PERCENT_DUPLICATION and ESTIMATED_LIBRARY_SIZE
        if recompute_merged_metrics:
            parsed_data["PERCENT_DUPLICATION"] = calculate_percentage_duplication(parsed_data)
            parsed_data["ESTIMATED_LIBRARY_SIZE"] = estimate_library_size(parsed_data)

        # Save the data
        if s_name in data_by_sample:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
        data_by_sample[s_name] = parsed_data
        module.add_data_source(f, s_name, section="DuplicationMetrics")

        # End of metrics table - reset for next sample
        if len(vals) < 6:
            return True

        # On to the next library if not merging
        else:
            return False

    # Go through logs and find Metrics
    for f in module.find_log_files(sp_key, filehandles=True):
        s_name = f["s_name"]
        parsed_lists: Dict[str, List] = defaultdict(list)
        keys = None
        in_stats_block = False
        recompute_merged_metrics = False

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="MarkDuplicates",
                sentieon_algo="Dedup",
            )
            if maybe_s_name:
                s_name = maybe_s_name

            if s_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class="DuplicationMetric", sentieon_algo="Dedup"):
                keys = f["f"].readline().strip("\n").split("\t")
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                in_stats_block = True

            # Currently parsing the METRICS table
            elif in_stats_block:
                vals = line.rstrip("\n").split("\t")

                # End of the METRICS table, or multiple libraries, and we're not merging them
                if len(vals) < 6 or (not merge_multiple_libraries and len(parsed_lists) > 0):
                    parsed_data = {k: parsed_list[0] for k, parsed_list in parsed_lists.items()}
                    if save_table_results(s_name, keys, parsed_data, recompute_merged_metrics):
                        # Reset for next file if returned True
                        s_name = None
                        parsed_lists = defaultdict(list)
                        keys = None
                        in_stats_block = False
                        recompute_merged_metrics = False
                    continue

                # Parse the column values
                if keys is not None:
                    assert len(vals) == len(keys), (keys, vals, f)
                    for k, v in zip(keys, vals):
                        # More than one library present and merging stats
                        if k in parsed_lists:
                            recompute_merged_metrics = True

                        v = v.strip()
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                        parsed_lists[k].append(v)

        parsed_data = {}
        for k in parsed_lists:
            # Sometimes a numerical column will an empty string, so converting "" to 0.0
            if all(isinstance(x, float) or x == "" for x in parsed_lists[k]):
                parsed_data[k] = sum(0.0 if x == "" else x for x in parsed_lists[k])
            else:
                parsed_data[k] = "/".join(str(x) for x in parsed_lists[k])

        # Files with no extra lines after last library
        if in_stats_block:
            save_table_results(s_name, keys, parsed_data, recompute_merged_metrics)

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_dups")

    # Add to general stats table
    headers = {
        "PERCENT_DUPLICATION": {
            "title": "Duplication",
            "description": "Mark Duplicates - Percent Duplication",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "OrRd",
            "modify": lambda x: util.multiply_hundred(x),
        }
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="Mark Duplicates")

    # Make the bar plot and add to the MarkDuplicates section
    #
    # The table in the Picard metrics file contains some columns referring
    # read pairs and some referring to single reads.
    for s_name, metr in data_by_sample.items():
        metr["READS_IN_DUPLICATE_PAIRS"] = 2.0 * metr["READ_PAIR_DUPLICATES"]
        metr["READS_IN_UNIQUE_PAIRS"] = 2.0 * (metr["READ_PAIRS_EXAMINED"] - metr["READ_PAIR_DUPLICATES"])
        metr["READS_IN_UNIQUE_UNPAIRED"] = metr["UNPAIRED_READS_EXAMINED"] - metr["UNPAIRED_READ_DUPLICATES"]
        metr["READS_IN_DUPLICATE_PAIRS_OPTICAL"] = 2.0 * metr["READ_PAIR_OPTICAL_DUPLICATES"]
        metr["READS_IN_DUPLICATE_PAIRS_NONOPTICAL"] = (
            metr["READS_IN_DUPLICATE_PAIRS"] - metr["READS_IN_DUPLICATE_PAIRS_OPTICAL"]
        )
        metr["READS_IN_DUPLICATE_UNPAIRED"] = metr["UNPAIRED_READ_DUPLICATES"]
        metr["READS_UNMAPPED"] = metr["UNMAPPED_READS"]

    keys = dict()
    keys_r = [
        "READS_IN_UNIQUE_PAIRS",
        "READS_IN_UNIQUE_UNPAIRED",
        "READS_IN_DUPLICATE_PAIRS_OPTICAL",
        "READS_IN_DUPLICATE_PAIRS_NONOPTICAL",
        "READS_IN_DUPLICATE_UNPAIRED",
        "READS_UNMAPPED",
    ]
    for k in keys_r:
        keys[k] = {"name": k.replace("READS_", "").replace("IN_", "").replace("_", " ").title()}

    # Config for the plot
    pconfig = {
        "id": f"{module.anchor}_deduplication",
        "title": f"{module.name}: Deduplication Stats",
        "ylab": "# Reads",
        "cpswitch_counts_label": "Number of Reads",
        "cpswitch_c_active": False,
    }

    module.add_section(
        name="Mark Duplicates",
        anchor=f"{module.anchor}-markduplicates",
        description="Number of reads, categorised by duplication state. **Pair counts "
        "are doubled** - see help text for details.",
        helptext="""
        The table in the Picard metrics file contains some columns referring
        read pairs and some referring to single reads.

        To make the numbers in this plot sum correctly, values referring to pairs are 
        doubled
        according to the scheme below:

        * `READS_IN_DUPLICATE_PAIRS = 2 * READ_PAIR_DUPLICATES`
        * `READS_IN_UNIQUE_PAIRS = 2 * (READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES)`
        * `READS_IN_UNIQUE_UNPAIRED = UNPAIRED_READS_EXAMINED - 
        UNPAIRED_READ_DUPLICATES`
        * `READS_IN_DUPLICATE_PAIRS_OPTICAL = 2 * READ_PAIR_OPTICAL_DUPLICATES`
        * `READS_IN_DUPLICATE_PAIRS_NONOPTICAL = READS_IN_DUPLICATE_PAIRS - 
        READS_IN_DUPLICATE_PAIRS_OPTICAL`
        * `READS_IN_DUPLICATE_UNPAIRED = UNPAIRED_READ_DUPLICATES`
        * `READS_UNMAPPED = UNMAPPED_READS`
        """,
        plot=bargraph.plot(data_by_sample, keys, pconfig),
    )

    # Return the number of detected samples to the parent module
    return data_by_sample.keys()


def calculate_percentage_duplication(d):
    """
    Picard calculation to compute percentage duplication

    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob
    /78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam
    /DuplicationMetrics.java#L106-L110
    """
    try:
        dups = d["UNPAIRED_READ_DUPLICATES"] + d["READ_PAIR_DUPLICATES"] * 2
        examined = d["UNPAIRED_READS_EXAMINED"] + d["READ_PAIRS_EXAMINED"] * 2
        if examined == 0:
            return 0
        return float(dups) / float(examined)
    except KeyError:
        return 0


def estimate_library_size(d):
    """
    Picard calculation to estimate library size

    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob
    /78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam
    /DuplicationMetrics.java#L153-L164

    Note: Optical duplicates are contained in duplicates and therefore do not enter
    the calculation here.
    See also the computation of READ_PAIR_NOT_OPTICAL_DUPLICATES.

     * Estimates the size of a library based on the number of paired end molecules
     observed
     * and the number of unique pairs observed.
     * <p>
     * Based on the Lander-Waterman equation that states:
     * C/X = 1 - exp( -N/X )
     * where
     * X = number of distinct molecules in library
     * N = number of read pairs
     * C = number of distinct fragments observed in read pairs
    """

    try:
        readPairs = d["READ_PAIRS_EXAMINED"] - d["READ_PAIR_OPTICAL_DUPLICATES"]
        uniqueReadPairs = d["READ_PAIRS_EXAMINED"] - d["READ_PAIR_DUPLICATES"]
    except KeyError:
        return None

    readPairDuplicates = readPairs - uniqueReadPairs

    if readPairs > 0 and readPairDuplicates > 0:
        m = 1.0
        M = 100.0

        if uniqueReadPairs >= readPairs or f(m * uniqueReadPairs, uniqueReadPairs, readPairs) < 0:
            logging.warning("Picard recalculation of ESTIMATED_LIBRARY_SIZE skipped - metrics " "look wrong")
            return None

        # find value of M, large enough to act as other side for bisection method
        while f(M * uniqueReadPairs, uniqueReadPairs, readPairs) > 0:
            M *= 10.0

        # use bisection method (no more than 40 times) to find solution
        for i in range(40):
            r = (m + M) / 2.0
            u = f(r * uniqueReadPairs, uniqueReadPairs, readPairs)
            if u == 0:
                break
            elif u > 0:
                m = r
            elif u < 0:
                M = r

        return uniqueReadPairs * (m + M) / 2.0
    else:
        return None


def f(x, c, n):
    """
    Picard calculation used when estimating library size

    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob
    /78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam
    /DuplicationMetrics.java#L172-L177

    * Method that is used in the computation of estimated library size.
    """
    return c / x - 1 + math.exp(-n / x)
