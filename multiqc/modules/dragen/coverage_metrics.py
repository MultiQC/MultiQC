#!/usr/bin/env python
# coding=utf-8

'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers coverage metrics data and prepares it for the output report.
Mostly implemented for the DRAGEN v3.7. Further improvements can/will be made.

It relies on the following official sources:
https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/CoverageMetricsReport_fDG.htm
https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QCMetricsCoverageReports_fDG_dtSW.htm

Section - 'DRAGEN DNA Pipeline', Subsection - 'QC Metrics and Coverage/Callability Reports':
https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf

The following sentence from the User Guide deserves attention:
"Additional coverage metrics can be enabled, and additional coverage regions can be specified."
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''

import re
import logging
from collections import defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils.util_functions import write_data_file

"""
The common dragen cleaner can be written later, for now it is placed locally,
but is disabled to provide the compatibility with other dragen modules.
from .utils import clean_sample
"""

NAMESPACE = "DRAGEN coverage"


'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The following code section provides a clean and simple tool for setting configurations:
https://github.com/ewels/MultiQC/blob/master/docs/plots.md#creating-a-table

The SINGLE_HEADER is used to define the common settings for all coverage headers.
'''
SINGLE_HEADER = {
    "namespace": NAMESPACE,        # Name for grouping. Prepends desc and is in Config Columns modal
    "title": "No title",           # Short title, table column title
    "description": "No",           # Longer description, goes in mouse hover text
    "max": None,                   # Minimum value in range, for bar / colour coding
    "min": 0,                      # Maximum value in range, for bar / colour coding
    "ceiling": None,               # Maximum value for automatic bar limit
    "floor": None,                 # Minimum value for automatic bar limit
    "minRange": None,              # Minimum range for automatic bar
    "scale": "GnBu",               # Colour scale for colour coding. False to disable.
    "bgcols": {},                  # Dict with values: background colours for categorical data.
    "colour": "15, 150, 255",      # Colour for column grouping
    "suffix": "",                  # Suffix for value (eg. "%")
    "format": "{:,.1f}",           # Value format string - default 1 decimal place
    "cond_formatting_rules": {},   # Rules for conditional formatting table cell values.
    "cond_formatting_colours": [], # Styles for conditional formatting of table cell values
    "shared_key": None,            # See the link for description
    "modify": None,                # Lambda function to modify values, special case, see below
    "hidden": True,                # Set to True to hide the column on page load
}
# The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER.
# These are just extensions, but can be used as if they were a part of the SINGLE_HEADER.
# These are not set by default. Only the SINGLE_HEADER is used as a basis for a header.
EXTRA_HEADER = {
    "hidden2": False,               # For non-general plots in the coverage section.
}


'''-----------------------------------------------------
The composition of the region strings is very sensitive.
It has to coincide with the standard formats of regions.
But the letter case is actually irrelevant.
'''
WGS = "genome"
QC = "QC coverage region"

# The desired regions can be added:
# TRGT = "target region"
BED  = "target bed" # May be incorrect. No real data was available.

# The REGIONS holds the variables. You can add common region's settings here.
# These configs have higher priority than those in the SINGLE_HEADER.
# They also overwrite the automatically created ones.
REGIONS = {
WGS:{
    "hidden": False,
    "scale":"Reds",
    "hidden2": False,
},
QC:{
    "scale":"Greens",
    "hidden2": False,
},
BED:{},
}

'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The METRICS is the container for "high level" representations of the standard metrics.
They are designed to be seen as simple identifiers of the well known metrics.
You can set any possible real/virtual configuration defined in the SINGLE_HEADER.

Some rules:

- You can use the single regions defined above as keys to set
  region-specific configurations for a certain metric.

- "modify" must be (lambda x: x if isinstance(x, str) else expression)
  Or generally speaking your function must check for string and return it.
  If not string (eg int or float), then do something with x and return the result.

- The descending precedence of configurations:
    1. Region-specific metric
    2. General metric
    3. General region
    4. Automatically defined settings
    5. Basic configs in the SINGLE_HEADER

- If you want to add new metrics, then you could try the following:
  - All characters must be in lower case. No exclusions.
  - No extra padded space on the left and right sides.
  - If tokens (eg words, operators) are separated from each other by a sequnce
    of spaces, replace each such sequnce with just one single-space string " ".
  - If the previous steps did not help, then you might need to take a look
    at the "make_metric_id" and/or "make_user_configs".


General info:

- All the configs defined below are for demonstration purposes only.
  You can adapt them according to your personal taste.

- You have the final word in defining settings.
  The code trusts you and does not check the chosen values.
'''

auto_configs_enabled = True
user_configs_enabled = True

V2 = "(second value)" # Used to create a unique ID if a metric has two values.

METRICS = {

    "aligned bases":{
        "title": "Aln bases", "scale": "RdYlGn", "colour": "0, 0, 255",
        "description": "Total number of aligned bases.", "hidden2": False,
    },
    "aligned bases in region":{
        "title": "Bases on target", "scale": "Reds", "colour": "0, 0, 255",

        WGS:{
            "scale": "Purples",
            "description": "Number of uniquely mapped bases to genome.",
        },
        QC:{"scale": "Greens",},
    },
    "aligned bases in region" + V2:{
        "title": "Bases on trg pct", "max": 100, "scale": "RdGy", "suffix": " %",
        "colour": "0, 0, 255", "bgcols": {"NA":"#00FFFF",},

        WGS:{
            "scale": "Purples",
            "description": "Percentage of uniquely mapped bases to genome.",
        },
        QC:{"scale": "Greens",},
    },


    "aligned reads":{
        "colour": "0, 255, 0", "scale": "PuBu", "title": "Aln reads",
        "description": "Total number of aligned reads.",
    },
    "aligned reads in region":{
        "title": "Reads on target", "scale": "RdGy", "colour": "0, 255, 0",
    },
    "aligned reads in region" + V2:{
        "title": "Reads on trg pct", "max": 100, "suffix": " %",
        "scale": "RdGy", "colour": "0, 255, 0",
    },


    "average alignment coverage over region":{
        "title": "Depth", "scale": "BrBG",
    },
    "average autosomal coverage over region":{
        "title": "Mean aut cov", "suffix": " x", "bgcols": {"NA": "#FF00FF"},
        WGS:{
            "description":(
                "Average autosomal coverage over genome. Calculated as the number of bases " +
                "that aligned to the autosomal loci in genome divided by the total number of loci " +
                "in the autosomal loci in genome. If there is no autosomes in the reference genome, " +
                "or the region does not intersect autosomes, this metric shows as NA."),
        },
    },
    "average mitochondrial coverage over region":{
        "title": "MT cov", "scale": "PRGn", "suffix": " x",
    },
    "average chr x coverage over region":{
        "title": "X cov", "scale": "Blues", "suffix": " x",
    },
    "average chr y coverage over region":{
        "title": "Y cov", "scale": "Blues", "suffix": " x",
    },

    # This metric can handle any float. '0.0' is just a tecnical trick.
    # The same configs are set for all possible floats.
    # Only the '>' operator and the '*mean' are supported.
    # Title and description shall not be set here, they are automated.
    "uniformity of coverage (pct > 0.0*mean) over region":{
        "scale": "RdBu", "suffix": " %", "colour": "255, 0, 0",
    },

    # The following two metrics are the most tedious. The creation of the title
    # and description are automated. If you want to target some specific case,
    # then you can specify a certain combination of (ix, inf) or (ix, jx).
    # Just write a group in the "extra", you can consider it as a separate
    # metric, so you can also add general or region-specific configs there.
    # '0.x' are just a tecnical trick. Any natural number is meant instead.

    # PCT of region with coverage [ix, inf)
    "pct of region with coverage [0x:inf)":{
        "max": 100, "scale": "Blues", "suffix": " %", "colour": "255, 50, 25",
        WGS:{
            "scale": "Purples",
        },
        QC:{
            "scale": "Greens",
        },
        "extra":{
            ("0x", "inf"):{
                "hidden2": False, "title": "⩾0x&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", "scale": "BrBG",
                "cond_formatting_rules":{
                    "red": [{"s_ne": 100.0}],
                    "green": [{"s_eq": 100.0}],
                },
                "cond_formatting_colours":[
                    {"red": "#FF0000"},
                    {"green": "#00FF00"},
                ],
                WGS:{
                    "scale": "Oranges",
                    "description": "Percentage of sites in genome with any coverage.",
                },
                QC:{
                    "description": "Percentage of sites in QC coverage region with any coverage.",
                },
                BED:{
                    "description": "Percentage of sites in target bed with any coverage.",
                },
            },
            ("1x", "inf"):{
                "hidden2": False, "title": "⩾1x&nbsp;&nbsp;",
            },
        },
    },
    # PCT of region with coverage [ix, jx)
    "pct of region with coverage [0x:0x)":{
        "max": 100, "scale": "Blues", "suffix": " %",
        "extra":{
            ("0x", "1x"):{
                "title": "0x&nbsp;&nbsp;&nbsp;",
                WGS:{
                    "scale": "Oranges",
                    "description": "Percentage of sites in genome with no coverage.",
                },
                QC:{
                    "description": "Percentage of sites in QC coverage region with no coverage.",
                },
                BED:{
                    "description": "Percentage of sites in target bed with no coverage.",
                },
            },
        },
    },


    "median autosomal coverage over region":{
        "title": "Med aut cov", "scale": "PuBuGn", "suffix": " x",
    },
    "mean/median autosomal coverage ratio over region":{
        "title": "Mean/med aut cov", "scale": "Blues", "suffix": " x",
    },
    "xavgcov/yavgcov ratio over region":{
        "title": "XAvgCov/YAvgCov", "scale": "Blues",
    },
    "xavgcov/autosomalavgcov ratio over region":{
        "title": "XAvgCov/AutAvgCov",
    },
    "yavgcov/autosomalavgcov ratio over region":{
        "title": "YAvgCov/AutAvgCov",
    },
}



##*****************************************************************************************##

##                               The actual code starts here

##*****************************************************************************************##

"""
The cov_data dictionary has the following structure:
  samples x regions x subregions x metrics x values, where:
  - sample is the <output prefix> as defined in the standard:
    <output prefix>.<coverage region prefix>_coverage_metrics.csv
    Samples are chosen to be the first layer because of how "base_module.ignore_samples" works.
  - region can be the genome, a target region, or a QC coverage region, maybe something else.
    Used only to automatically create a plot or set of plots for each region.
  - subregion is the <coverage region prefix> as defined in the standard.

The cov_headers dictionary stores all found metrics with associated configurations.
"""
cov_data = defaultdict(lambda: defaultdict(dict))
cov_headers = defaultdict(dict)

"""
The dragen.py provides a common interface for all dragen modules.
To avoid cluttering up the global dragen space, only one method is defined.
Other methods can be added as well to provide modularity.
Collected coverage data can be also transferred to other dragen modules as an attribute.
"""
class DragenCoverageMetrics(BaseMultiqcModule):
    def add_coverage_metrics(self):
        """Called from dragen.py. Does different abstract tasks."""

        global cov_data

        for file in self.find_log_files("dragen/coverage_metrics"):
            if parse(file): self.add_data_source(file, section = "stats")

        if not cov_data: return set()

        '''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Coverage data is collected now. The headers are also complete and available for use.
        If other files (eg json) are needed, then these can be found with find_log_files,
        processed and used afterwards. See afterqc.py or supernova.py for example.
        -------------------------------------------------------------------------------------
        CSS/JS can be imported, see fastqc.py for details.
        -------------------------------------------------------------------------------------
        The ignore_samples and clean_s_name use the Multiqc's config module.

        Assumption: some/most users are not or dont want to be familiar with the cleaning
        details of the code and most will set configurations for the original sample names.
        If it is true, then the proper sequence of execution steps would be:
        1. Ignore samples.
        2. Use the global clean_s_name.
        3. Use the local dragen cleaner.

        If not, then the sequence can be always rearranged.
        --------------------------------------------------------------------------------------
        In theory some combinations of sample names, configurations and implementation aspects
        of both cleaners may produce duplicates. Cleaned names are checked for duplicates, if
        presented then their original names are used instead. The issue will be reported.
        """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''
        cov_data = self.ignore_samples(cov_data)
        check_cleaned_names(
        {s_name: clean_sample(self.clean_s_name(s_name)) for s_name in cov_data.keys()})

        # Write data into the general table.
        gen_data, gen_headers = make_general_stats()
        self.general_stats_addcols(gen_data, gen_headers, namespace = NAMESPACE)

        # Write general data to file.
        self.write_data_file(gen_data, "dragen_cov_metrics")

        # Report all found errors and warnings.
        make_log_report()

        self.add_section(
            name = "Coverage metrics",
            anchor = "dragen-cov-metrics",
            description = """
            Coverage metrics over a region (where the region can be a target region,
            a QC coverage region, or the whole genome). Press the `Help` button for details.
            """,
            helptext = """
            The following criteria are used when calculating coverage:

            * Duplicate reads and clipped bases are ignored.
            * Only reads with `MAPQ` > `min MAPQ` and bases with `BQ` > `min BQ` are considered

            Considering only bases usable for variant calling, _i.e._ excluding:

            1. Clipped bases
            2. Bases in duplicate reads
            3. Reads with `MAPQ` < `min MAPQ` (default `20`)
            4. Bases with `BQ` < `min BQ` (default `10`)
            5. Reads with `MAPQ` = `0` (multimappers)
            6. Overlapping mates are double-counted
            """,
            plot = make_own_plots(),
        )
        return cov_data.keys()


"""
The cov_patterns list contains the regexes used to match the metrics. The regexes must be properly sorted.
All regexes are constructed to be as general and simple as possible to speed up the matching of metrics.

The "make_auto_configs" dictionary is closely related to the cov_patterns.
The PCT_RGX and UNI_RGX are also used in the "make_user_configs".
"""
PCT_RGX = "^PCT of (?P<region>.+?) with coverage (?P<entity>.+)"
AVG_RGX = "^Average (?P<entity>.+?) coverage over (?P<region>.+)"
ALN_RGX = "^Aligned (?P<entity>.+)"
UNI_RGX = "^Uniformity of coverage (?P<entity>.+?) over (?P<region>.+)"
MED_RGX = "^Median (?P<entity>.+?) coverage over (?P<region>.+)"
RAT_RGX = "(?P<entity>.+?) ratio over (?P<region>.+)"
ANY_RGX = ".+"

cov_patterns = [PCT_RGX, AVG_RGX, ALN_RGX, UNI_RGX, MED_RGX, RAT_RGX, ANY_RGX]

# Initialise the logger.
log = logging.getLogger(__name__)

"""
The log_data is the common holder for errors and warnings,
which are stored in groups to make the output log file more readable.
errors:
- invalid_file_names, which do not conform to the:
    <output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv
- duplicate_file_names, which theoretically can be found in subdirectories.
- duplicate_clean_sample_names is a special case for catching bugs associated with cleaning.
  Details are in the DragenCoverageMetrics.
- invalid_file_lines, which do not conform to COVERAGE SUMMARY,,<metric>,<value1>  or
                                              COVERAGE SUMMARY,,<metric>,<value1>,<value2>
warnings:
- unknown_metrics are not presented in the standard.
- unusual_values are those except for int/float/NA.
"""
log_data = {
    "invalid_file_names": defaultdict(list),
    "duplicate_file_names": defaultdict(list),
    "duplicate_clean_sample_names": {},
    "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
    "unknown_metrics": defaultdict(lambda: defaultdict(list)),
    "unusual_values": defaultdict(lambda: defaultdict(dict))
}

# Used to simplify the code and speed up execution.
cov_headers_support = defaultdict(dict)

def parse(file_handler):
    """Parser for coverage metrics csv files.
    Input:  file_handler, which gives access to the necessary information - file name/content/root
    Output: 0 as failure or 1 as success in extracting of metrics."""

    """ cov_data, cov_headers, log_data are updated within the function body.

    File name is checked below, because there is no guarantee that it would be as in the official standard:
    <output prefix>.<coverage region prefix>_coverage_metrics.csv

    Accepted structure of files:
    <output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv

    Some possible file names:

    T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv
    T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv
    T_SRR7890936_50pc.target_bed_coverage_metrics.csv
    SAMPLE_MONO_0_dragen.qc-coverage-region-1_coverage_metrics.csv
    SAMPLE_MONO_0_dragen.qc-coverage-region-2_coverage_metrics.csv
    SAMPLE_MONO_1_dragen.wgs_coverage_metrics.csv
    """
    root_name = file_handler["root"]
    file_name = file_handler["fn"]

    file_match = re.search(r"^([^.]+)\.(.+?)_coverage_metrics(.*).csv", file_name)
    if not file_match:
        log_data["invalid_file_names"][root_name].append(file_name)
        return 0

    sample, subregion = file_match.group(1), file_match.group(2)

    # subregion is just concatenated with an arbitrary group, because it is non-standard.
    # subregion is used instead of sample, because 2D plots can be created only if sample names are the same.
    if file_match.group(3): subregion += file_match.group(3)

    for region in cov_data[sample]:
        if subregion in cov_data[sample][region]:
            log_data["duplicate_file_names"][root_name].append(file_name)
            return 0

    success = 0
    region = ""
    data = {}
    for line in file_handler["f"].splitlines():

        # Check the general structure. Maybe a new section will be added in the future.
        line_match = re.search("^COVERAGE SUMMARY,,([^,]+),([^,]+),?([^,]*)", line)

        # If line is fine then extract the necessary fields.
        if line_match:
            metric, value1, value2 = line_match.group(1), line_match.group(2), line_match.group(3)
            success = 1

        # Otherwise check if line is empty. If not then report it and go to the next line.
        else:
            if not re.search("^\s*$", line): log_data["invalid_file_lines"][root_name][file_name].append(line)
            continue

        """ Check the extracted values. These are int/float and in some cases NA for now. All other values will be reported.
        The type conversion is only performed to int and float. Other types of data are left unchanged as strings.
        It has some consequences for the "modify" configuration. Details are in the description for the "make_auto_configs".
        No additional try blocks are required to support other simple types of values (eg Inf, None). Just adapt the regexes.
        """
        try:
            value1 = int(value1)
        except ValueError:
            try:
                value1 = float(value1)
            except ValueError:
                if not re.search("^NA$", value1, re.IGNORECASE):
                    log_data["unusual_values"][root_name][file_name][metric + V2] = value1

        if value2:
            try:
                value2 = float(value2)
            except ValueError:
                try:
                    value2 = int(value2)
                except ValueError:
                    if not re.search("^NA$", value2, re.IGNORECASE):
                        log_data["unusual_values"][root_name][file_name][metric + V2] = value2

        """
        Unfortunately there is no guarantee for consistency in the structure of metrics.
        So each metric must be modified. The details are in the next big comment.
        """
        metric_id = make_metric_id(metric)

        if metric_id not in cov_headers:
            for metric_pattern in cov_patterns:
                metric_match = re.search(metric_pattern, metric_id, re.IGNORECASE)
                if metric_match:
                    output = make_auto_configs[metric_pattern](metric_match)

                    if "region" in output:
                        region = output["region"]
                        cov_headers_support[metric_id]["region"] = region

                    if "warning" in output: cov_headers_support[metric_id]["warning"] = 1

                    configs = get_std_configs(SINGLE_HEADER)

                    if auto_configs_enabled:
                        auto_configs = output["configs"]
                        configs.update(auto_configs)

                    if user_configs_enabled:
                        user_configs = make_user_configs(metric_id, region)
                        configs.update(user_configs)

                    cov_headers[metric_id] = configs

                    if value2:
                        configs = get_std_configs(SINGLE_HEADER)

                        if auto_configs_enabled and "configs2" in output:
                            auto_configs = output["configs2"]
                            configs.update(auto_configs)

                        if user_configs_enabled:
                            user_configs = make_user_configs(metric_id + V2, region)
                            configs.update(user_configs)

                        cov_headers[metric_id + V2] = configs

                    break

        if not region and "region" in cov_headers_support[metric_id]:
            region = cov_headers_support[metric_id]["region"]

        if "warning" in cov_headers_support[metric_id]:
            log_data["unknown_metrics"][root_name][file_name].append(metric)

        data[metric_id] = value1
        if value2: data[metric_id + V2] = value2

    if success:
        if not region: region = "Unknown region"
        cov_data[sample][region][subregion] = data

    return success



'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The DRAGEN v3.7 standard does not provide any info about consistency of coverage metrics.
In practice two metrics may be inconsistent:
    PCT of region with coverage [ix, inf)
    PCT of region with coverage [ix, jx)

Presumably some part of the DRAGEN tries to improve the appearance of content in csv files.
The two mentioned metrics may have a different number of spaces inserted within them.
Amount and insertion places seem to depend on the highest number of digits in ix and jx.
So, for example:
    if max ix = jx = 100x   [100x: inf), [ 50x:100x)   then [  0x: inf),  [ 10x: inf)
                                                        and [  0x:  1x),  [ 10x: 15x)
but if max ix = jx = 1500x  [1500x: inf),[1000x:1500x) then [   0x: inf), [  10x: inf)
                                                        and [   0x:   1x),[  10x:  15x)

Using such versions of the same metric as keys would result in two columns instead of one.
There is a possibility that an input folder contains a mixed data with different versions.

Assumption: each metric may be inconsistent due to unknown reasons.
Solution:   modify each metric.

To not exaggerate the problem only the most obvious possible issues are solved.
If some new metric-specific peculiarities arise, they can be added to the make_metric_id.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''

def make_metric_id(metric):
    """Tries to fix consistency issues."""

    metric_id = re.sub("\s+", " ", metric).strip() # Single backslashes are escaped if presented.
    metric_id = metric_id.lower()

    pct_case = re.search(PCT_RGX, metric_id, re.IGNORECASE)
    if pct_case:
        metric_id = "pct of " + pct_case["region"] + " with coverage " + pct_case["entity"].replace(" ", "")
        return metric_id

    return metric_id


def make_user_configs(metric, region):
    """Prepares the user-defined configurations."""

    configs = {}

    # First check if empty. If not set general region-specific configs.
    if region:
        metric = metric.replace(region, "region")
        region = get_REGION(region)
        if region in REGIONS and REGIONS[region]:
            configs.update(get_std_configs(REGIONS[region]))

    pct_case = re.search(PCT_RGX, metric, re.IGNORECASE)
    if pct_case or re.search(UNI_RGX, metric, re.IGNORECASE):
        metric = re.sub("\d+", "0", metric)

    # Now set region-specific parameters for a given metric.
    if metric in METRICS:
        metric = METRICS[metric]
        configs.update(get_std_configs(metric))

        # If region exists and not empty.
        if region and region in metric:
            configs.update(get_std_configs(metric[region]))

    # Try to set (ix:inf)/(ix:jx)-specific settings.
    if pct_case and "extra" in metric and metric["extra"]:

        ix_jx = re.search(r"\[(\d+x):(inf|\d+x)\)", pct_case["entity"])
        if ix_jx:
            _group = (ix_jx.group(1), ix_jx.group(2))

            if _group in metric["extra"]:
                metric = metric["extra"][_group]
                configs.update(get_std_configs(metric))
                if region in metric: configs.update(get_std_configs(metric[region]))

    return configs

def get_std_configs(configs):
    """Copies the standard real/virtual configurations"""
    std_configs = {}

    for key, val in configs.items():
        if key in SINGLE_HEADER or key in EXTRA_HEADER:
            if val is not None: std_configs[key] = val

    return std_configs

def make_general_stats():
    """Prepare data and headers for the general stats."""

    gen_data = defaultdict(dict)
    gen_headers = defaultdict(dict)

    for sample in cov_data:
        for region in cov_data[sample]:
            for subregion in cov_data[sample][region]:
                for metric in cov_data[sample][region][subregion]:
                    new_metric_id = metric + subregion
                    gen_data[sample][new_metric_id] = cov_data[sample][region][subregion][metric]
                    gen_headers[new_metric_id]={
                        key: val for key, val in cov_headers[metric].items() if key not in EXTRA_HEADER
                    }
                    subregion_cleaned = clean_subregion(subregion, mode = "Column header")
                    if "title" in gen_headers[new_metric_id]:
                        gen_headers[new_metric_id]["title"]=(
                            re.sub("&nbsp;", "", gen_headers[new_metric_id]["title"]) + " " + subregion_cleaned
                        )
    return gen_data, gen_headers


def make_own_plots():
    """Create plots for the non-general section."""
    """
    This function prepares a plot for each subregion, that is <coverage region prefix>

    The "table_config" feature will be implemented later.

    Data may be transfered directly if needed:
    "<script>DRAGEN_data={var}</script>".format(var = json.dumps(var))
    """
    plots = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    for sample in cov_data:
        for region in cov_data[sample]:
            reg = cov_data[sample][region]
            for subregion in reg:
                data = defaultdict(dict)
                headers = defaultdict(dict)
                for metric in reg[subregion]:
                    data[sample][metric] = reg[subregion][metric]
                    headers[metric] = get_std_configs(cov_headers[metric])
                    if "hidden2" in headers[metric]:
                        headers[metric]["hidden"] = headers[metric]["hidden2"]
                        del headers[metric]["hidden2"]

                if sample not in plots[region][subregion]["data"]:
                    plots[region][subregion]["data"][sample] = data[sample]

                for head in headers:
                    if head not in plots[region][subregion]["headers"]:
                        plots[region][subregion]["headers"][head] = headers[head]


    html = ""
    for region in plots:
        if re.search("genome", region, re.IGNORECASE):
            for subregion in plots[region]:
                html = ("<hr><p style='text-align:center;'><b>" +
                        subregion.upper() + "</b></p>" +
                        table.plot(plots[region][subregion]["data"],
                                   plots[region][subregion]["headers"],
                                   {"namespace": NAMESPACE, "table_title": subregion})
                        + html)
            continue

        for subregion in plots[region]:
            html += ("<hr><p style='text-align:center;'><b>" +
                     clean_subregion(subregion, mode = "Table header") + "</b></p>" +
                     table.plot(plots[region][subregion]["data"],
                                plots[region][subregion]["headers"],
                               {"namespace": NAMESPACE, "table_title": subregion}
                     )
            )
    return html



##-----------------------------------------------------------------------------##
##                              Various functions                              ##
##-----------------------------------------------------------------------------##

def check_cleaned_names(sample_names):
    """Check cleaned sample names for duplicates."""

    temp = defaultdict(list)
    for orig_s_name in sample_names: temp[sample_names[orig_s_name]].append(orig_s_name)

    for clean_s_name in temp:
        if len(temp[clean_s_name]) > 1:
            log_data["duplicate_clean_sample_names"][clean_s_name] = temp[clean_s_name]
        else:
            orig_s_name = temp[clean_s_name][0]
            if clean_s_name != orig_s_name: cov_data[clean_s_name] = cov_data.pop(orig_s_name)


def make_log_report():
    """The only purpose of this function is to create a readable and informative log output."""

    if log_data["invalid_file_names"]:
        log_message = ("\n\nThe file names must conform to the following structure:\n" +
                       "<output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv\n\n" +
                       "The following files are not valid:\n")

        for root in log_data["invalid_file_names"]:
            log_message += "  " + root + ":\n"
            for file in log_data["invalid_file_names"][root]:
                log_message += "    " + file + "\n"

        log.error(log_message + "\n")

    if log_data["duplicate_file_names"]:
        log_message = "\n\nThe following duplicates were found:\n"

        for root in log_data["duplicate_file_names"]:
            log_message += "  " + root + ":\n"
            for file in log_data["duplicate_file_names"][root]:
                log_message += "    " + file + "\n"

        log.error(log_message + "\n")

    if log_data["duplicate_clean_sample_names"]:
        log_message = ("\n\nThe following sample names are left as they are " +
                       "because their cleaned versions are duplicates.\n")

        for clean_s_name in log_data["duplicate_clean_sample_names"]:
            log_message += "  " + clean_s_name + " is formed from:\n"
            for orig_s_name in log_data["duplicate_clean_sample_names"][clean_s_name]:
                log_message += "    " + orig_s_name + "\n"

        log.error(log_message + "\n")

    if log_data["invalid_file_lines"]:
        log_message = ("\n\nThe lines in files must be:\n" +
                       "COVERAGE SUMMARY,,<metric>,<value1> or " +
                       "COVERAGE SUMMARY,,<metric>,<value1>,<value2>\n\n"  +
                       "The following files contain invalid lines:\n" )

        for root in log_data["invalid_file_lines"]:
            log_message += "  " + root + ":\n"
            for file in log_data["invalid_file_lines"][root]:
                log_message += "    " + file + ":\n"
                for line in log_data["invalid_file_lines"][root][file]:
                    log_message += "      " + line + "\n"

        log.error(log_message + "\n")

    if log_data["unknown_metrics"]:
        log_message = "\n\nThe following files contain unknown metrics:\n"

        for root in log_data["unknown_metrics"]:
            log_message += "  " + root + ":\n"
            for file in log_data["unknown_metrics"][root]:
                log_message += "    " + file + ":\n"
                for metric in log_data["unknown_metrics"][root][file]:
                    log_message += "      " + metric + "\n"

        log.warning(log_message + "\n")

    if log_data["unusual_values"]:
        log_message = ("\n\nAll metrics' values except int, float and NA are non-standard.\n"  +
                       "The following files contain non-standard values:\n" )

        for root in log_data["unusual_values"]:
            log_message += "  " + root + ":\n"
            for file in log_data["unusual_values"][root]:
                log_message += "    " + file + ":\n"
                for metric in log_data["unusual_values"][root][file]:
                    log_message += ("      " + metric + " = " +
                                    log_data["unusual_values"][root][file][metric] + "\n")

        log.warning(log_message + "\n")


def clean_subregion(subregion, mode = ""):
    """Change <coverage region prefix> to make the output table better looking.
    Input:  subregion, that is the <coverage region prefix> in the official standard.
    Output: a cleaned version of the subregion."""

    subregion = subregion.strip()

    if re.search("QC", subregion, re.IGNORECASE):
        subregion = re.sub("qc", "QC", subregion, flags = re.I)
        if mode == "Table header":
            subregion = re.sub(r"-|_|\.", " ", subregion)
        elif mode == "Column header":
            subregion = re.sub(r"-|_|\.", "", subregion)
            subregion = re.sub("coverage|region", "", subregion, flags = re.I)
        return subregion

    subregion = re.sub(r"-|_|\.", " ", subregion)
    return subregion


# This cleaner will be exported into the ./dragen/utils.py later.
# Turned off for compatibility with other dragen modules.
def clean_sample(sample):
    """Cleans sample names to make the output tables more attractive."""
    return sample

    sample = re.sub("dragen|sample", "", sample, flags = re.I)
    sample = re.sub("-|_", " ", sample)
    sample = sample.strip()

    return sample


"""
The following functions are mainly used to create configurations for columns of html tables.

Input:
- match_object, which stores the result of matching a metric.

Output:
A dictionary with keys:
- obligatory configs  is a dictionary with real/logical configurations from the SINGLE_HEADER.
- obligatory configs2 if the second value is presented or in case when metric can not be recognized.
- optional   region   if presented in a metric and may be extracted easily. Returning it is encouraged.
- optional   warning  may be returned if a metric could not be recognized. The associated value is irrelevant.

Please note: modify lambda/func must check for string input.
"""

##---------------------------------------------------##
## This little segment was taken from ./dragen/utils.py
## as an adoptation for the previous code version.
## It is left unchanged.
read_format = "{:,.1f}"
read_count_multiplier = config.read_count_multiplier
if read_count_multiplier == 1:
    read_format = "{:,.0f}"

base_format = "{:,.1f}&nbsp;"
base_count_multiplier = config.base_count_multiplier
if base_count_multiplier == 1:
    base_format = "{:,.0f}"
elif base_count_multiplier == 0.000000001:
    base_format = "{:,.2f}"
##---------------------------------------------------##

def get_REGION(region):
    """Matches input region with well known regions."""
    for RG in REGIONS:
        if re.search(RG, region, re.IGNORECASE): return RG
    return region

def get_Aligned_configs(match_object):
    """The following metrics are supported:
    Aligned bases, Aligned reads
    Aligned bases in region, Aligned reads in region
    """
    configs = {}
    entity = match_object["entity"].lower()

    metric_match = re.search("^(bases|reads)$", entity)
    if metric_match:
        configs.update({
            "min": 0, "format": base_format if entity == "bases" else read_format,
            "description": "Total number of aligned " + entity + ".",
            "title": (config.base_count_prefix if entity == "bases" else config.read_count_prefix) + " Aln " + entity
        })
        if entity == "reads": configs["modify"] = lambda x: x if isinstance(x, str) else x * read_count_multiplier
        else:                 configs["modify"] = lambda x: x if isinstance(x, str) else x * base_count_multiplier

        return {"configs": configs}

    configs2 = {}
    metric_match = re.search("^(?P<entity>bases|reads) in (?P<region>.+)", match_object["entity"])
    if metric_match:
        entity = metric_match["entity"].lower()
        region = get_REGION(metric_match["region"])
        configs.update({
            "min": 0, "format": base_format if entity == "bases" else read_format,
            "title": (config.base_count_prefix if entity == "bases" else config.read_count_prefix) + " Aln " + entity,
            "description":(
                ("Number of uniquely mapped bases to " + region + ".") if entity == "bases" else
                ("Number of uniquely mapped reads to " + region + ". When " + region + " is the target BED, " +
                 "this metric is equivalent to and replaces Capture Specificity based on target region."))
        })
        if entity == "reads": configs["modify"] = lambda x: x if isinstance(x, str) else x * read_count_multiplier
        else:                 configs["modify"] = lambda x: x if isinstance(x, str) else x * base_count_multiplier

        configs2.update({"min": 0, "max": 100, "suffix": " %",
                        "format": base_format if entity == "bases" else read_format,
                        "description": "Percentage of aligned " + entity + " in " + region + ".",
                        "title": "Aln " + entity + " on trg"})

        return {"configs": configs, "configs2": configs2, "region": metric_match["region"]}

    # Else unrecognized.
    configs["title"] = "Aln " + entity
    configs2["title"] = "Aln " + entity + " value2"
    return {"configs": configs, "configs2": configs2, "warning": 1}

def get_Average_configs(match_object):
    """The following metrics are supported:
    Average alignment coverage over region
    Average chromosome X coverage over region, Average chromosome Y coverage over region
    Average mitochondrial coverage over region
    Average autosomal coverage over region
    """
    configs = {}
    entity = match_object["entity"]
    region = get_REGION(match_object["region"])


    entity_match = re.search("^(alignment|chr x|chr y|mitochondrial|autosomal)$", entity)
    if entity_match:
        configs.update({"suffix": " x", "min": 0, "hidden": True, "title": "Depth"})

        if entity == "alignment":
            configs["description"] = ("Number of uniquely mapped bases to " + region +
                                      " divided by the number of sites in " + region + ".")


        elif entity == "autosomal":
            configs["description"] = (
                "Total number of bases that aligned to the autosomal loci in " + region +
                " divided by the total number of loci in the autosomal loci in " + region +
                ". If there is no autosome in the reference genome, or the " + region +
                " does not intersect autosomes, this metric shows as NA."
            )
            configs["title"] = "Avg aut cov"

        elif entity == "mitochondrial":
            configs["description"] = (
                "Total number of bases that aligned to the intersection of the mitochondrial chromosome with " + region +
                " divided by the total number of loci in the intersection of the mitochondrial chromosome with " + region +
                ". If there is no mitochondrial chromosome in the reference genome or the " + region +
                " does not intersect mitochondrial chromosome, this metric shows as NA."
            )
            configs["title"] = "Avg mit cov"

        else:
            entity = re.search("(x|y)", entity).group(1).capitalize()
            configs["title"] = "Avg " + entity + " cov"
            entity = "chromosome " + entity
            configs["description"] = (
                "Total number of bases that aligned to the intersection of " + entity +
                " with " + region + " divided by the total number of loci in the intersection of " +
                entity + " with " + region + ". If there is no " + entity + " in the reference genome or the " +
                region + " does not intersect " + entity + ", this metric shows as NA.")
        return {"configs": configs, "region": match_object["region"]}

    # Else unrecognized.
    entity = match_object["entity"]
    configs["title"] = "Avg " + entity
    configs2 = {"title": "Avg " + entity + " value2"}
    return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

def get_Uniformity_configs(match_object):
    """The following metrics are supported:
    Uniformity of coverage (PCT > float*mean) over region
    """
    configs = {}

    entity = match_object["entity"]
    region = get_REGION(match_object["region"])


    entity_match = re.search(r"\(PCT\s*>\s*(\d+\.\d+)\*mean\)", entity, re.IGNORECASE)
    if entity_match:
        multiplier = entity_match.group(1)
        percent = str(float(multiplier)*100) + "%"
        configs.update({"suffix": " %", "min": 0, "max": 100, "title": ">" + multiplier + "*mean",
                        "description": "Percentage of sites with coverage greater than " +
                        percent + " of the mean coverage in " + region + ".",})

        return {"configs": configs, "region": match_object["region"]}

    # Else unrecognized.
    configs2 = {}
    configs["title"] = "Uniformity " + entity
    configs2["title"] = "Uniformity " + entity + " value2"
    return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

def get_PCT_configs(match_object):
    """The following metrics are supported:
    PCT of region with coverage [ix, inf)
    PCT of region with coverage [ix, jx)
    """
    configs = {}
    entity = match_object["entity"]
    region = get_REGION(match_object["region"])

    entity_match = re.search(r"\[(\d+)x:(inf|\d+)x?\)", entity)
    if entity_match:
        left_ix  = entity_match.group(1)
        right_jx = entity_match.group(2)
        configs.update({
            "suffix": " %", "min": 0, "max": 100, "hidden": True, "title": entity,
            "description":(
                "Percentage of sites in " + region + " with no coverage." if (left_ix == "0" and right_jx == "1") else
               ("Percentage of sites in " + region + " with at least " + left_ix + "x coverage.")
                if re.search("inf", right_jx, re.IGNORECASE) else
               ("Percentage of sites in " + region + " with at least " + left_ix + "x but less than " + right_jx + "x coverage.")),
        })
        return {"configs": configs, "region": match_object["region"]}

    # Else unrecognized.
    configs2 = {}
    configs["title"] ="PCT " + entity
    configs2["title"] = "PCT " + entity + " value2"
    return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

def get_Median_configs(match_object):
    """The following metrics are supported:
    Median autosomal coverage over region
    """
    configs = {}
    entity = match_object["entity"]
    region = get_REGION(match_object["region"])

    metric_match = re.search("^(autosomal)$", entity)
    if metric_match:
        configs.update({"suffix": " x", "min": 0, "hidden": True, "title": "Med aut cov",
                        "description": ("Median alignment coverage over the autosomal loci in " + region +
                                        ". If there is no autosome in the reference genome or the " + region +
                                        " does not intersect autosomes, this metric shows as NA."),
        })
        return {"configs": configs, "region": match_object["region"]}

    # Else unrecognized.
    configs2 = {}
    configs["title"] = "Med " + entity
    configs2["title"] = "Med " + entity + " value2"
    return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

def get_ratio_configs(match_object):
    """The following metrics are supported:
    Mean/Median autosomal coverage ratio over region
    XAvgCov/YAvgCov ratio over region
    XAvgCov/AutosomalAvgCov ratio over region
    YAvgCov/AutosomalAvgCov ratio over region
    """
    configs = {}
    entity = match_object["entity"]
    region = get_REGION(match_object["region"])

    entity_match = re.search("^Mean/Median autosomal coverage$", entity, re.IGNORECASE)
    if entity_match:
        configs.update({"suffix": " x", "title": "Mean/med aut cov", "min": 0, "hidden": True,
                        "description": "Mean autosomal coverage in " + region +
                        " divided by the median autosomal coverage in " + region +
                        ". If there is no autosome in the reference genome or the " + region +
                        " does not intersect autosomes, this metric shows as NA."
        })
        return {"configs": configs, "region": match_object["region"]}

    entity_match = re.search("^(XAvgCov|YAvgCov|AutosomalAvgCov)/(XAvgCov|YAvgCov|AutosomalAvgCov)$", entity, re.IGNORECASE)
    if entity_match:
        devident = re.search("(X|Y|Aut)", entity_match.group(1), re.IGNORECASE).group(1).capitalize()
        devisor  = re.search("(X|Y|Aut)", entity_match.group(2), re.IGNORECASE).group(1).capitalize()
        configs.update({"min": 0, "hidden": True, "title": devident + "AvgCov/" + devisor + "AvgCov"})

        chrom = ("chromosome " + devident) if (len(devident) == 1) else "autosomal chromosome"
        chrom2  = ("chromosome " + devisor)  if (len(devisor)  == 1) else "autosomal chromosome"
        configs["description"] = ("Average " + chrom + " alignment coverage in " + region +
                                  " divided by the average " + chrom2 + " alignment coverage in  " + region)

        return {"configs": configs, "region": match_object["region"]}

    # Else unrecognized.
    configs2 = {}
    configs["title"] = entity
    configs2["title"] = entity + " value2"
    return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}


def get_any_configs(match_object):
    """This function handles all other metrics, whose general structure could not be recognized."""
    configs = {"title": match_object.string, "scale": "Purples",}
    configs2 = {"title":  match_object.string + " value2", "scale": "Reds",}
    return {"configs": configs, "configs2": configs2, "warning": 1}


'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The make_auto_configs is a dictionary with keys being the regexes defined in the "cov_patterns".
Each regex is associated with function, which automatically produces configurations from the "SINGLE_HEADER".
It always catches a metric, but it is not capable to produce good configs for all possible variations of metrics.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''
make_auto_configs = {
    ALN_RGX: get_Aligned_configs,
    AVG_RGX: get_Average_configs,
    PCT_RGX: get_PCT_configs,
    UNI_RGX: get_Uniformity_configs,
    MED_RGX: get_Median_configs,
    RAT_RGX: get_ratio_configs,
    ANY_RGX: get_any_configs
}
