import logging
import os
import re
from typing import Dict, List, Optional, Union

from multiqc.utils import config

# Initialise the logger
log = logging.getLogger(__name__)


def read_histogram(self, program_key, headers, formats, picard_tool, sentieon_algo=None):
    """
    Reads a Picard HISTOGRAM file.

    Args:
        self: the Picard QC module
        program_key: the key used to find the program (ex. picard/quality_by_cycle)
        headers: the list of expected headers for the histogram
        formats: the list of methods to apply to re-format each field (on a given row)
        picard_tool: the name of the Picard tool to be found in the header, e.g. MeanQualityByCycle
        sentieon_algo: the name of the Sentieon algorithm to be found in the header, e.g. MeanQualityByCycle
    """
    all_data = dict()
    assert len(formats) == len(headers)
    sample_data = None

    # Go through logs and find Metrics
    for f in self.find_log_files(program_key, filehandles=True):
        s_name = f["s_name"]
        for line in f["f"]:
            maybe_s_name = self.extract_sample_name(line, f, picard_tool=picard_tool, sentieon_algo=sentieon_algo)
            if maybe_s_name:
                s_name = maybe_s_name
                sample_data = None

            if self.is_line_right_before_table(line, sentieon_algo=sentieon_algo):
                # check the header
                line = f["f"].readline()
                if line.strip().split("\t") == headers:
                    sample_data = dict()
                else:
                    sample_data = None

            elif sample_data is not None:
                fields = line.strip().split("\t")
                if len(fields) == len(headers):
                    for i in range(len(fields)):
                        fields[i] = formats[i](fields[i])
                    sample_data[fields[0]] = dict(zip(headers, fields))

        # append the data
        if sample_data:
            if s_name in all_data:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
            all_data[s_name] = sample_data

            self.add_data_source(f, s_name, section="Histogram")
            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

    data = self.ignore_samples(all_data)

    # Write data to file
    self.write_data_file(data, f"{self.anchor}_histogram")

    return data


def is_line_right_before_table(
    line: str,
    picard_class: Union[None, str, List[str]] = None,
    sentieon_algo: Optional[str] = None,
) -> bool:
    """
    Picard logs from different samples can be concatenated together, so the module
    needs to know a marker to find where new sample information starts.

    The command line Picard tools themselves use the "## METRICS CLASS" header
    line for that purpose; however, the Picard classes often used by other
    tools and platforms - e.g. Sentieon and Parabricks  - while adding their own
    headers, so we need to handle them as well.
    """
    if isinstance(picard_class, list):
        picard_classes = picard_class
    elif picard_class is None:
        picard_classes = []
    else:
        picard_classes = [picard_class]
    return (
        (line.startswith("## METRICS CLASS") or line.startswith("## HISTOGRAM"))
        and (not picard_classes or any(c in line for c in picard_classes))
        or sentieon_algo
        and line.startswith("#SentieonCommandLine:")
        and f" --algo {sentieon_algo}" in line
    )


def extract_sample_name(
    mod,
    line: str,
    f: Dict,
    picard_tool: str,
    sentieon_algo: Optional[str] = None,
) -> Optional[str]:
    """
    Historically, MultiQC supported Picard tools QC outputs merged together into
    one file from different samples and tools. In order to handle this correctly,
    MultiQC needs to take the sample name elsewhere rather than the file name,
    so it tries to parse the command line recorded in the output header.

    This approach can be disabled with `config.picard_config.s_name_filenames`
    """
    if getattr(config, "picard_config", {}).get("s_name_filenames", False):
        return None
    picard_command = (
        line.startswith("# ")
        and picard_tool in line
        and (
            " input=" in line.lower()
            or " --input=" in line.lower()
            or "input" in line.lower().split()
            or "--input" in line.lower().split()
        )
    )
    sentieon_command = (
        sentieon_algo
        and line.startswith("#SentieonCommandLine:")
        and f" --algo {sentieon_algo}" in line
        and " -i " in line
    )
    # Pull sample name from the input file name, recorded in the command line:
    fn_search = None
    if picard_command:
        fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", line, flags=re.IGNORECASE)
    elif sentieon_command:
        fn_search = re.search(r" -i\s+(\[?\S+\]?)", line, flags=re.IGNORECASE)
    if fn_search:
        f_name = os.path.basename(fn_search.group(1).strip("[]"))
        s_name = mod.clean_s_name(f_name, f)
        return s_name
    return None


def multiply_hundred(val):
    try:
        val = float(val) * 100
    except ValueError:
        pass
    return val
