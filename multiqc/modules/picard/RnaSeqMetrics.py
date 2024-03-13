""" MultiQC submodule to parse output from Picard RnaSeqMetrics """

import logging

from multiqc.modules.picard import util
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard RnaSeqMetrics reports and parse their data"""

    data_by_sample = dict()
    histogram_by_sample = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/rnaseqmetrics", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]
        in_hist = False

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="RnaSeqMetrics",
            )
            if maybe_s_name:
                s_name = maybe_s_name

            if s_name is None:
                continue

            # Catch the histogram values
            if in_hist:
                try:
                    sections = line.split("\t")
                    pos = int(sections[0])
                    coverage = float(sections[1])
                    histogram_by_sample[s_name][pos] = coverage
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False

            if util.is_line_right_before_table(line, picard_class="RnaSeqMetrics"):
                keys = f["f"].readline().strip("\n").split("\t")
                vals = f["f"].readline().strip("\n").split("\t")
                if len(vals) != len(keys):
                    continue

                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

                module.add_data_source(f, s_name, section="RnaSeqMetrics")
                data_by_sample[s_name] = dict()
                histogram_by_sample[s_name] = dict()

                for k, v in zip(keys, vals):
                    if not v:
                        v = "NA"
                    else:
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                        else:
                            # Multiply percentages by 100
                            if k.startswith("PCT_"):
                                v = v * 100.0
                    data_by_sample[s_name][k] = v
                # Calculate some extra numbers
                if "PF_BASES" in keys and "PF_ALIGNED_BASES" in keys:
                    data_by_sample[s_name]["PF_NOT_ALIGNED_BASES"] = (
                        data_by_sample[s_name]["PF_BASES"] - data_by_sample[s_name]["PF_ALIGNED_BASES"]
                    )

            elif line.startswith("## HISTOGRAM"):
                keys = f["f"].readline().strip("\n").split("\t")
                assert len(keys) >= 2, (keys, f)
                in_hist = True
                histogram_by_sample[s_name] = dict()

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    histogram_by_sample = module.ignore_samples(histogram_by_sample)
    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, "multiqc_picard_RnaSeqMetrics")

    # Add to general stats table
    headers = {
        "PCT_RIBOSOMAL_BASES": {
            "title": "rRNA",
            "description": "Percent of aligned bases overlapping ribosomal RNA regions",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "Reds",
        },
        "PCT_MRNA_BASES": {
            "title": "mRNA",
            "description": "Percent of aligned bases overlapping UTRs and coding regions of mRNA transcripts",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "Greens",
        },
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="RnaSeqMetrics")

    # Bar plot of bases assignment
    bg_cats = {
        "CODING_BASES": {"name": "Coding"},
        "UTR_BASES": {"name": "UTR"},
        "INTRONIC_BASES": {"name": "Intronic"},
        "INTERGENIC_BASES": {"name": "Intergenic"},
        "RIBOSOMAL_BASES": {"name": "Ribosomal"},
        "PF_NOT_ALIGNED_BASES": {"name": "PF not aligned"},
    }

    # Warn user if any samples are missing 'RIBOSOMAL_BASES' data; ie picard was run without an rRNA interval file.
    warn_rrna = ""
    rrna_missing = []
    for s_name, metrics in data_by_sample.items():
        if metrics["RIBOSOMAL_BASES"] == "NA":
            rrna_missing.append(s_name)
    if rrna_missing:
        if len(rrna_missing) < 5:
            missing_samples = f"for samples <code>{'</code>, <code>'.join(rrna_missing)}</code>"
        else:
            missing_samples = f"<strong>{len(rrna_missing)} samples</strong>"
        warn_rrna = f"""
        <div class="alert alert-warning">
          <span class="glyphicon glyphicon-warning-sign"></span>
          Picard was run without an rRNA annotation file {missing_samples}, therefore the ribosomal assignment is not available. To correct, rerun with the <code>RIBOSOMAL_INTERVALS</code> parameter, as documented <a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics" target="_blank">here</a>.
        </div>
        """

    pconfig = {
        "id": "picard_rnaseqmetrics_assignment_plot",
        "title": "Picard: RnaSeqMetrics Base Assignments",
        "ylab": "Number of bases",
    }
    module.add_section(
        name="RnaSeqMetrics Assignment",
        anchor="picard-rna-assignment",
        description="Number of bases in primary alignments that align to regions in the reference genome." + warn_rrna,
        plot=bargraph.plot(data_by_sample, bg_cats, pconfig),
    )

    # Bar plot of strand mapping
    bg_cats = dict()
    bg_cats["CORRECT_STRAND_READS"] = {"name": "Correct"}
    bg_cats["INCORRECT_STRAND_READS"] = {"name": "Incorrect", "color": "#8e123c"}

    pdata = dict()
    for s_name, d in data_by_sample.items():
        if d["CORRECT_STRAND_READS"] > 0 and d["INCORRECT_STRAND_READS"] > 0:
            pdata[s_name] = d
    if len(pdata) > 0:
        pconfig = {
            "id": "picard_rnaseqmetrics_strand_plot",
            "title": "Picard: RnaSeqMetrics Strand Mapping",
            "ylab": "Number of reads",
            "hide_zero_cats": False,
        }
        module.add_section(
            name="RnaSeqMetrics Strand Mapping",
            anchor="picard-rna-strand",
            description="Number of aligned reads that map to the correct strand.",
            plot=bargraph.plot(data_by_sample, bg_cats, pconfig),
        )

    # Section with histogram plot
    if len(histogram_by_sample) > 0:
        # Plot the data and add section
        pconfig = {
            "smooth_points": 500,
            "smooth_points_sumcounts": [True, False],
            "id": "picard_rna_coverage",
            "title": "Picard: Normalized Gene Coverage",
            "ylab": "Coverage",
            "xlab": "Percent through gene",
            "xDecimals": False,
            "tt_label": "<b>{point.x}%</b>: {point.y:.0f}",
            "ymin": 0,
        }
        module.add_section(
            name="Gene Coverage",
            anchor="picard-rna-coverage",
            plot=linegraph.plot(histogram_by_sample, pconfig),
        )

    # Return the number of detected samples to the parent module
    return len(data_by_sample)
