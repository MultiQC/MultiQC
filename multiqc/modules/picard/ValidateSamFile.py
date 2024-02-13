""" MultiQC submodule to parse output from Picard ValidateSamFile """

import logging

from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)

# Possible warnings and descriptions
WARNING_DESCRIPTIONS = {
    "ADJACENT_INDEL_IN_CIGAR": "CIGAR string contains an insertion (I) followed by " "deletion (D), or vice versa",
    "BAM_FILE_MISSING_TERMINATOR_BLOCK": "BAM appears to be healthy, but is an older "
    "file so doesn't have terminator block",
    "E2_BASE_EQUALS_PRIMARY_BASE": "Secondary base calls should not be the same as "
    "primary, unless one or the other is N",
    "INVALID_DATE_STRING": "Date string is not ISO-8601",
    "INVALID_QUALITY_FORMAT": "Quality encodings out of range; appear to be Solexa or "
    "Illumina when Phred expected. Avoid exception being "
    "thrown as a result of no qualities being read.",
    "MISSING_TAG_NM": "The NM tag (nucleotide differences) is missing",
    "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND": "Pair flag set but not marked as " "first or second of pair",
    "RECORD_MISSING_READ_GROUP": "A SAMRecord is found with no read group id",
}

# All possible errors and descriptions
ERROR_DESCRIPTIONS = {
    "CIGAR_MAPS_OFF_REFERENCE": "Bases corresponding to M operator in CIGAR extend " "beyond reference",
    "DUPLICATE_PROGRAM_GROUP_ID": "Same program group id appears more than once",
    "DUPLICATE_READ_GROUP_ID": "Same read group id appears more than once",
    "EMPTY_READ": "Indicates that a read corresponding to the first strand has a "
    "length of zero and/or lacks flow signal intensities (FZ)",
    "HEADER_RECORD_MISSING_REQUIRED_TAG": "Header tag missing in header line",
    "HEADER_TAG_MULTIPLY_DEFINED": "Header tag appears more than once in header line " "" "" "with different value",
    "INVALID_ALIGNMENT_START": "Alignment start position is incorrect",
    "INVALID_CIGAR": "CIGAR string error for either read or mate",
    "INVALID_FLAG_FIRST_OF_PAIR": "First of pair flag set for unpaired read",
    "INVALID_FLAG_MATE_NEG_STRAND": "Mate negative strand flag set for unpaired read",
    "INVALID_FLAG_MATE_UNMAPPED": "Mate unmapped flag is incorrectly set",
    "INVALID_FLAG_NOT_PRIM_ALIGNMENT": "Not primary alignment flag set for unmapped " "read",
    "INVALID_FLAG_PROPER_PAIR": "Proper pair flag set for unpaired read",
    "INVALID_FLAG_READ_UNMAPPED": "Mapped read flat not set for mapped read",
    "INVALID_FLAG_SECOND_OF_PAIR": "Second of pair flag set for unpaired read",
    "INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT": "Supplementary alignment flag set for " "unmapped read",
    "INVALID_INDEX_FILE_POINTER": "Invalid virtualFilePointer in index",
    "INVALID_INDEXING_BIN": "Indexing bin set on SAMRecord does not agree with " "computed value",
    "INVALID_INSERT_SIZE": "Inferred insert size is out of range",
    "INVALID_MAPPING_QUALITY": "Mapping quality set for unmapped read or is >= 256",
    "INVALID_MATE_REF_INDEX": "Mate reference index (MRNM) set for unpaired read",
    "INVALID_PLATFORM_VALUE": "The read group has an invalid value set for its PL " "field",
    "INVALID_PREDICTED_MEDIAN_INSERT_SIZE": "PI tag value is not numeric",
    "INVALID_REFERENCE_INDEX": "Reference index not found in sequence dictionary",
    "INVALID_TAG_NM": "The NM tag (nucleotide differences) is incorrect",
    "INVALID_VERSION_NUMBER": "Does not match any of the acceptable versions",
    "MATE_CIGAR_STRING_INVALID_PRESENCE": "A cigar string for a read whose mate is " "NOT mapped",
    "MATE_FIELD_MISMATCH": "Read alignment fields do not match its mate",
    "MATE_NOT_FOUND": "Read is marked as paired, but its pair was not found",
    "MATES_ARE_SAME_END": "Both mates of a pair are marked either as first or second " "" "" "mates",
    "MISMATCH_FLAG_MATE_NEG_STRAND": "Mate negative strand flag does not match read " "strand flag",
    "MISMATCH_FLAG_MATE_UNMAPPED": "Mate unmapped flag does not match read unmapped " "flag of mate",
    "MISMATCH_MATE_ALIGNMENT_START": "Mate alignment does not match alignment start " "of mate",
    "MISMATCH_MATE_CIGAR_STRING": "The mate cigar tag does not match its mate's cigar " "" "" "string",
    "MISMATCH_MATE_REF_INDEX": "Mate reference index (MRNM) does not match reference " "" "" "index of mate",
    "MISMATCH_READ_LENGTH_AND_E2_LENGTH": "Lengths of secondary base calls tag " "values" " and read should match",
    "MISMATCH_READ_LENGTH_AND_QUALS_LENGTH": "Length of sequence string and length of"
    " base quality string do not match",
    "MISMATCH_READ_LENGTH_AND_U2_LENGTH": "Secondary base quals tag values should " "match read length",
    "MISSING_HEADER": "The SAM/BAM file is missing the header",
    "MISSING_PLATFORM_VALUE": "The read group is missing its PL (platform unit) field",
    "MISSING_READ_GROUP": "The header is missing read group information",
    "MISSING_SEQUENCE_DICTIONARY": "There is no sequence dictionary in the header",
    "MISSING_VERSION_NUMBER": "Header has no version number",
    "POORLY_FORMATTED_HEADER_TAG": "Header tag does not have colon",
    "READ_GROUP_NOT_FOUND": "A read group ID on a SAMRecord is not found in the header",
    "RECORD_OUT_OF_ORDER": "The record is out of order",
    "TAG_VALUE_TOO_LARGE": "Unsigned integer tag value is deprecated in BAM. Template " "" "" "length",
    "TRUNCATED_FILE": "BAM file does not have terminator block",
    "UNRECOGNIZED_HEADER_TYPE": "Header record is not one of the standard types",
}


def _default_data_entry():
    return {"WARNING_count": 0, "ERROR_count": 0, "file_validation_status": "pass"}


def parse_reports(module):
    """
    Find Picard ValidateSamFile reports and parse their data based on wether we
    think it's a VERBOSE or SUMMARY report
    """

    # Get data
    data = _parse_reports_by_type(module)

    if data:
        #  Filter to strip out ignored sample names (REQUIRED)
        data = module.ignore_samples(data)

        # Populate the general stats table
        module.general_stats_addcols(data, _get_general_stats_headers(), namespace="ValidateSamFile")

        # Add any found data to the report
        _add_section_to_report(module, data)

        # Write parsed data to a file
        module.write_data_file(data, "multiqc_picard_validatesamfile")

    return len(data)


def _parse_reports_by_type(module):
    """Returns a data dictionary

    Goes through logs and parses them based on 'No errors found', VERBOSE or SUMMARY
    type.
    """

    data_by_sample = dict()

    for f in module.find_log_files("picard/sam_file_validation", filehandles=True):
        module.add_data_source(f, "ValidateSamFile")

        s_name = f["s_name"]

        if s_name in data_by_sample:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

        fh = f["f"]
        first_line = fh.readline().rstrip()
        fh.seek(0)  # Rewind reading of the file

        if "No errors found" in first_line:
            sample_data = _parse_no_error_report()
        elif first_line.startswith("ERROR") or first_line.startswith("WARNING"):
            sample_data = _parse_verbose_report(fh)
        else:
            sample_data = _parse_summary_report(fh)

        data_by_sample[s_name] = sample_data

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        module.add_software_version(None, s_name)

    return data_by_sample


def _parse_no_error_report():
    return _default_data_entry()


def _parse_verbose_report(filehandle):
    sample_data = _default_data_entry()
    errors = warnings = 0
    for line in filehandle:
        if line.startswith("WARNING:"):
            warnings += 1
        elif line.startswith("ERROR:"):
            errors += 1

    sample_data["WARNING_count"] = warnings
    sample_data["ERROR_count"] = errors

    if errors:
        sample_data["file_validation_status"] = "fail"
    elif warnings:
        sample_data["file_validation_status"] = "warn"

    return sample_data


def _parse_summary_report(filehandle):
    sample_data = _default_data_entry()
    for problem_type, name, count in _histogram_data(filehandle):
        sample_data[name] = count
        sample_data[problem_type + "_count"] += count

    if sample_data["ERROR_count"]:
        sample_data["file_validation_status"] = "fail"
    elif sample_data["WARNING_count"]:
        sample_data["file_validation_status"] = "warn"

    return sample_data


def _histogram_data(iterator):
    """Yields only the row contents that contain the histogram entries"""
    histogram_started = False
    header_passed = False
    for line in iterator:
        if "## HISTOGRAM" in line:
            histogram_started = True
        elif histogram_started:
            if header_passed:
                values = line.rstrip().split("\t")
                if len(values) == 1:
                    continue
                try:
                    problem_type, name = values[0].split(":")
                except ValueError:
                    log.warning(
                        "Line did not look like normal picard 'ERROR:NAME' format, " "" "" "ignoring: {}".format(
                            values[0]
                        )
                    )
                    continue
                yield problem_type, name, int(values[1])
            elif line.startswith("Error Type"):
                header_passed = True


def _add_section_to_report(module, data):
    """
    Adds found data to the report via several HTML generators
    """

    # Count samples that have errors and/or warnings
    pass_count = error_count = only_warning_count = 0
    for sample_data in data.values():
        if sample_data["file_validation_status"] == "fail":
            error_count += 1
        elif sample_data["file_validation_status"] == "warn":
            only_warning_count += 1
        else:
            pass_count += 1

    # Add overview note
    plot_html = []
    note_html = _generate_overview_note(
        pass_count=pass_count, only_warning_count=only_warning_count, error_count=error_count, total_count=len(data)
    )
    plot_html.append(note_html)

    # Add the detailed table, but only if we have something to show
    if error_count or only_warning_count:
        table_html = _generate_detailed_table(data)
        plot_html.append(table_html)

    # Finally, add the html to the report as a section
    module.add_section(
        name="SAM/BAM File Validation",
        anchor="picard_validatesamfile",
        description=(
            "This tool reports on the validity of a SAM or BAM " "file relative to " "the SAM-format specification."
        ),
        helptext="""
            A detailed table is only shown if errors or warnings are found. Details 
            about the errors and warnings are only shown if a `SUMMARY` report was 
            parsed.

            For more information on the warnings, errors and possible fixes please 
            read [this broadinstitute article](
            https://software.broadinstitute.org/gatk/documentation/article.php?id
            =7571).""",
        plot="\n".join(plot_html),
    )


def _get_general_stats_headers():
    """
    Returns a header dict for the general stats columns
    """
    headers = {
        "file_validation_status": {
            "title": "Validation",
            "description": "ValidateSamFile (File Validation)",
        },
        "WARNING_count": {
            "title": "# Warnings",
            "description": "ValidateSamFile (number of warnings)",
            "scale": "Oranges",
            "shared_key": "ValidateSamEntries",
            "colour": "255,237,160",
            "format": "{:.0f}",
            "hidden": True,
        },
        "ERROR_count": {
            "title": "# Errors",
            "description": "ValidateSamFile (number of errors)",
            "scale": "Reds",
            "shared_key": "ValidateSamEntries",
            "colour": "252,146,114",
            "format": "{:.0f}",
            "hidden": True,
        },
    }

    return headers


def _generate_overview_note(pass_count, only_warning_count, error_count, total_count):
    """Generates and returns the HTML note that provides a summary of validation
    status."""

    note_html = ['<div class="progress">']
    pbars = [
        [float(error_count), "danger", "had errors"],
        [float(only_warning_count), "warning", "had warnings"],
        [float(pass_count), "success", "passed"],
    ]
    for b in pbars:
        if b[0]:
            note_html.append(
                '<div class="progress-bar progress-bar-{pbcol}" style="width: {pct}%" '
                'data-toggle="tooltip" title="{count} {sample} {txt}">{'
                "count}</div>".format(
                    pbcol=b[1],
                    count=int(b[0]),
                    pct=(b[0] / float(total_count)) * 100.0,
                    txt=b[2],
                    sample="samples" if b[0] > 1 else "sample",
                )
            )
    note_html.append("</div>")

    return "\n".join(note_html)


def _generate_detailed_table(data):
    """
    Generates and returns the HTML table that overviews the details found.
    """
    headers = _get_general_stats_headers()

    # Only add headers for errors/warnings we have found
    for problems in data.values():
        for problem in problems:
            if problem not in headers and problem in WARNING_DESCRIPTIONS:
                headers["WARNING_count"]["hidden"] = False
                headers[problem] = {
                    "description": WARNING_DESCRIPTIONS[problem],
                    "namespace": "WARNING",
                    "scale": headers["WARNING_count"]["scale"],
                    "format": "{:.0f}",
                    "shared_key": "warnings",
                    "hidden": True,  # Hide by default; to unclutter things.
                }
            if problem not in headers and problem in ERROR_DESCRIPTIONS:
                headers["ERROR_count"]["hidden"] = False
                headers[problem] = {
                    "description": ERROR_DESCRIPTIONS[problem],
                    "namespace": "ERROR",
                    "scale": headers["ERROR_count"]["scale"],
                    "format": "{:.0f}",
                    "shared_key": "errors",
                    "hidden": True,  # Hide by default; to unclutter things.
                }

    table_config = {
        "table_title": "Picard: SAM/BAM File Validation",
        "id": "picard_validatesamfile_table",
    }

    return table.plot(data=data, headers=headers, pconfig=table_config)
