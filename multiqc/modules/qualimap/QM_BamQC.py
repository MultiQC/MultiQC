""" MultiQC Submodule to parse output from Qualimap BamQC """


import logging
import math
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Qualimap BamQC reports and parse their data"""

    # General stats - genome_results.txt
    self.qualimap_bamqc_genome_results = dict()
    for f in self.find_log_files("qualimap/bamqc/genome_results"):
        parse_genome_results(self, f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None, f["s_name"])

    self.qualimap_bamqc_genome_results = self.ignore_samples(self.qualimap_bamqc_genome_results)
    if len(self.qualimap_bamqc_genome_results) > 0:
        self.write_data_file(self.qualimap_bamqc_genome_results, "multiqc_qualimap_bamqc_genome_results")

    # Coverage - coverage_histogram.txt
    self.qualimap_bamqc_coverage_hist = dict()
    for f in self.find_log_files("qualimap/bamqc/coverage", filehandles=True):
        parse_coverage(self, f)
    self.qualimap_bamqc_coverage_hist = self.ignore_samples(self.qualimap_bamqc_coverage_hist)

    # Insert size - insert_size_histogram.txt
    self.qualimap_bamqc_insert_size_hist = dict()
    for f in self.find_log_files("qualimap/bamqc/insert_size", filehandles=True):
        parse_insert_size(self, f)
    self.qualimap_bamqc_insert_size_hist = self.ignore_samples(self.qualimap_bamqc_insert_size_hist)

    # GC distribution - mapped_reads_gc-content_distribution.txt
    self.qualimap_bamqc_gc_content_dist = dict()
    self.qualimap_bamqc_gc_by_species = dict()  # {'HUMAN': data_dict, 'MOUSE': data_dict}
    for f in self.find_log_files("qualimap/bamqc/gc_dist", filehandles=True):
        parse_gc_dist(self, f)
    self.qualimap_bamqc_gc_content_dist = self.ignore_samples(self.qualimap_bamqc_gc_content_dist)
    self.qualimap_bamqc_gc_by_species = self.ignore_samples(self.qualimap_bamqc_gc_by_species)

    num_parsed = max(
        len(self.qualimap_bamqc_genome_results),
        len(self.qualimap_bamqc_coverage_hist),
        len(self.qualimap_bamqc_insert_size_hist),
        len(self.qualimap_bamqc_gc_content_dist),
    )
    # Go no further if nothing found
    if num_parsed == 0:
        return 0

    try:
        covs = config.qualimap_config["general_stats_coverage"]
        assert type(covs) == list
        assert len(covs) > 0
        covs = [str(i) for i in covs]
        log.debug("Custom Qualimap thresholds: {}".format(", ".join([i for i in covs])))
    except (AttributeError, TypeError, AssertionError):
        covs = [1, 5, 10, 30, 50]
        covs = [str(i) for i in covs]
        log.debug("Using default Qualimap thresholds: {}".format(", ".join([i for i in covs])))
    self.covs = covs

    # Make the plots for the report
    report_sections(self)

    # Set up the general stats table
    general_stats_headers(self)

    # Return the number of reports we found
    return num_parsed


def parse_genome_results(self, f):
    """Parse the contents of the Qualimap BamQC genome_results.txt file"""
    regexes = {
        "Input": {
            "bam_file": r"bam file = (.+)",
        },
        "Globals": {
            "total_reads": r"number of reads = ([\d,]+)",
            "mapped_reads": r"number of mapped reads = ([\d,]+)",
            "mapped_bases": r"number of mapped bases = ([\d,]+)",
            "sequenced_bases": r"number of sequenced bases = ([\d,]+)",
        },
        "Insert size": {
            "mean_insert_size": r"mean insert size = ([\d,\.]+)",
            "median_insert_size": r"median insert size = ([\d,\.]+)",
        },
        "Mapping quality": {
            "mean_mapping_quality": r"mean mapping quality = ([\d,\.]+)",
        },
        "Mismatches and indels": {
            "general_error_rate": r"general error rate = ([\d,\.]+)",
        },
        "Coverage": {
            "mean_coverage": r"mean coverageData = ([\d,\.]+)",
        },
        "Globals inside": {
            "regions_size": r"regions size = ([\d,\.]+)",
            "regions_mapped_reads": r"number of mapped reads = ([\d,]+)",  # WARNING: Same as in Globals
        },
    }
    d = dict()
    section = None
    for line in f["f"].splitlines():
        if line.startswith(">>>>>>>"):
            section = line[8:]
        elif section:
            for k, r in regexes.get(section, {}).items():
                r_search = re.search(r, line)
                if r_search:
                    if "\d" in r:
                        try:
                            d[k] = float(r_search.group(1).replace(",", ""))
                        except ValueError:
                            d[k] = r_search.group(1)
                    else:
                        d[k] = r_search.group(1)

    # Check we have an input filename
    if "bam_file" not in d:
        log.debug("Couldn't find an input filename in genome_results file {}".format(f["fn"]))
        return None

    # Get a nice sample name
    s_name = self.clean_s_name(d["bam_file"], f)

    # Add to general stats table & calculate a nice % aligned
    try:
        self.general_stats_data[s_name]["total_reads"] = d["total_reads"]
        self.general_stats_data[s_name]["mapped_reads"] = d["mapped_reads"]
        self.general_stats_data[s_name]["general_error_rate"] = d["general_error_rate"] * 100.0
        self.general_stats_data[s_name]["mean_coverage"] = d["mean_coverage"]
        self.general_stats_data[s_name]["regions_size"] = d["regions_size"]
        self.general_stats_data[s_name]["regions_mapped_reads"] = d["regions_mapped_reads"]
        try:
            d["percentage_aligned"] = (d["mapped_reads"] / d["total_reads"]) * 100.0
            self.general_stats_data[s_name]["percentage_aligned"] = d["percentage_aligned"]
            d["percentage_aligned_on_target"] = (d["regions_mapped_reads"] / d["mapped_reads"]) * 100.0
            self.general_stats_data[s_name]["percentage_aligned_on_target"] = d["percentage_aligned_on_target"]
        except ZeroDivisionError:
            pass
    except KeyError:
        pass

    # Save results
    if s_name in self.qualimap_bamqc_genome_results:
        log.debug("Duplicate genome results sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_genome_results[s_name] = d
    self.add_data_source(f, s_name=s_name, section="genome_results")


def parse_coverage(self, f):
    """Parse the contents of the Qualimap BamQC Coverage Histogram file"""
    # Get the sample name from the parent parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/coverage_histogram.txt
    s_name = self.get_s_name(f)

    d = dict()
    for l in f["f"]:
        if l.startswith("#"):
            continue
        coverage, count = l.split(None, 1)
        coverage = int(round(float(coverage)))
        count = float(count)
        d[coverage] = count

    if len(d) == 0:
        log.debug("Couldn't parse contents of coverage histogram file {}".format(f["fn"]))
        return None

    # Find median without importing anything to do it for us
    num_counts = sum(d.values())
    cum_counts = 0
    total_cov = 0
    median_coverage = None
    for thiscov, thiscount in d.items():
        cum_counts += thiscount
        total_cov += thiscov * thiscount
        if cum_counts >= num_counts / 2:
            median_coverage = thiscov
            break
    self.general_stats_data[s_name]["median_coverage"] = median_coverage
    # Save results
    if s_name in self.qualimap_bamqc_coverage_hist:
        log.debug("Duplicate coverage histogram sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_coverage_hist[s_name] = d
    self.add_data_source(f, s_name=s_name, section="coverage_histogram")


def parse_insert_size(self, f):
    """Parse the contents of the Qualimap BamQC Insert Size Histogram file"""
    # Get the sample name from the parent parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/insert_size_histogram.txt
    s_name = self.get_s_name(f)

    d = dict()
    zero_insertsize = 0
    for l in f["f"]:
        if l.startswith("#"):
            continue
        insertsize, count = l.split(None, 1)
        insertsize = int(round(float(insertsize)))
        count = float(count) / 1000000
        if insertsize == 0:
            zero_insertsize = count
        else:
            d[insertsize] = count

    # Find median without importing anything to do it for us
    num_counts = sum(d.values())
    cum_counts = 0
    median_insert_size = None
    for thisins, thiscount in d.items():
        cum_counts += thiscount
        if cum_counts >= num_counts / 2:
            median_insert_size = thisins
            break
    # Add the median insert size to the general stats table
    self.general_stats_data[s_name]["median_insert_size"] = median_insert_size

    # Save results
    if s_name in self.qualimap_bamqc_insert_size_hist:
        log.debug("Duplicate insert size histogram sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_insert_size_hist[s_name] = d
    self.add_data_source(f, s_name=s_name, section="insert_size_histogram")


def parse_gc_dist(self, f):
    """Parse the contents of the Qualimap BamQC Mapped Reads GC content distribution file"""
    # Get the sample name from the parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt
    s_name = self.get_s_name(f)

    d = dict()
    reference_species = None
    reference_d = dict()
    avg_gc = 0
    for l in f["f"]:
        if l.startswith("#"):
            sections = l.strip("\n").split("\t", 3)
            if len(sections) > 2:
                reference_species = sections[2]
            continue
        sections = l.strip("\n").split("\t", 3)
        gc = int(round(float(sections[0])))
        content = float(sections[1])
        avg_gc += gc * content
        d[gc] = content
        if len(sections) > 2:
            reference_content = float(sections[2])
            reference_d[gc] = reference_content

    # Add average GC to the general stats table
    self.general_stats_data[s_name]["avg_gc"] = avg_gc

    # Save results
    if s_name in self.qualimap_bamqc_gc_content_dist:
        log.debug("Duplicate Mapped Reads GC content distribution sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_gc_content_dist[s_name] = d
    if reference_species and reference_species not in self.qualimap_bamqc_gc_by_species:
        self.qualimap_bamqc_gc_by_species[reference_species] = reference_d
    self.add_data_source(f, s_name=s_name, section="mapped_gc_distribution")


coverage_histogram_helptext = """
For a set of DNA or RNA reads mapped to a reference sequence, such as a genome
or transcriptome, the depth of coverage at a given base position is the number
of high-quality reads that map to the reference at that position
(<a href="https://doi.org/10.1038/nrg3642" target="_blank">Sims et al. 2014</a>).

Bases of a reference sequence (y-axis) are groupped by their depth of coverage
(*0&#215;, 1&#215;, &#8230;, N&#215;*) (x-axis). This plot shows
the frequency of coverage depths relative to the reference sequence for each
read dataset, which provides an indirect measure of the level and variation of
coverage depth in the corresponding sequenced sample.

If reads are randomly distributed across the reference sequence, this plot
should resemble a Poisson distribution (<a href="https://doi.org/10.1016/0888-7543(88)90007-9"
target="_blank">Lander & Waterman 1988</a>), with a peak indicating approximate
depth of coverage, and more uniform coverage depth being reflected in a narrower
spread. The optimal level of coverage depth depends on the aims of the
experiment, though it should at minimum be sufficiently high to adequately
address the biological question; greater uniformity of coverage is generally
desirable, because it increases breadth of coverage for a given depth of
coverage, allowing equivalent results to be achieved at a lower sequencing depth
(<a href="https://doi.org/10.1002/gepi.20575" target="_blank">Sampson
et al. 2011</a>; <a href="https://doi.org/10.1038/nrg3642" target="_blank">Sims
et al. 2014</a>). However, it is difficult to achieve uniform coverage
depth in practice, due to biases introduced during sample preparation
(<a href="https://doi.org/10.1016/j.yexcr.2014.01.008" target="_blank">van
Dijk et al. 2014</a>), sequencing (<a href="https://doi.org/10.1186/gb-2013-14-5-r51"
target="_blank">Ross et al. 2013</a>) and read mapping
(<a href="https://doi.org/10.1038/nrg3642" target="_blank">Sims et al. 2014</a>).

This plot may include a small peak for regions of the reference sequence with
zero depth of coverage. Such regions may be absent from the given sample (due
to a deletion or structural rearrangement), present in the sample but not
successfully sequenced (due to bias in sequencing or preparation), or sequenced
but not successfully mapped to the reference (due to the choice of mapping
algorithm, the presence of repeat sequences, or mismatches caused by variants
or sequencing errors). Related factors cause most datasets to contain some
unmapped reads (<a href="https://doi.org/10.1038/nrg3642" target="_blank">Sims
et al. 2014</a>)."""

genome_fraction_helptext = """
For a set of DNA or RNA reads mapped to a reference sequence, such as a genome
or transcriptome, the depth of coverage at a given base position is the number
of high-quality reads that map to the reference at that position, while the
breadth of coverage is the fraction of the reference sequence to which reads
have been mapped with at least a given depth of coverage
(<a href="https://doi.org/10.1038/nrg3642" target="_blank">Sims et al. 2014</a>).

Defining coverage breadth in terms of coverage depth is useful, because
sequencing experiments typically require a specific minimum depth of coverage
over the region of interest (<a href="https://doi.org/10.1038/nrg3642"
target="_blank">Sims et al. 2014</a>), so the extent of the reference sequence
that is amenable to analysis is constrained to lie within regions that have
sufficient depth. With inadequate sequencing breadth, it can be difficult to
distinguish the absence of a biological feature (such as a gene) from a lack
of data (<a href="https://doi.org/10.1101/gr.7050807" target="_blank">Green 2007</a>).

For increasing coverage depths (*1&#215;, 2&#215;, &#8230;, N&#215;*),
coverage breadth is calculated as the percentage of the reference
sequence that is covered by at least that number of reads, then plots
coverage breadth (y-axis) against coverage depth (x-axis). This plot
shows the relationship between sequencing depth and breadth for each read
dataset, which can be used to gauge, for example, the likely effect of a
minimum depth filter on the fraction of a genome available for analysis."""


def report_sections(self):
    """Add results from Qualimap BamQC parsing to the report"""
    # Append to self.sections list

    if len(self.qualimap_bamqc_coverage_hist) > 0:
        # Chew back on histogram to prevent long flat tail
        # (find a sensible max x - lose 1% of longest tail)
        max_x = 0
        total_bases_by_sample = dict()
        for s_name, d in self.qualimap_bamqc_coverage_hist.items():
            total_bases_by_sample[s_name] = sum(d.values())
            cumulative = 0
            for count in sorted(d.keys(), reverse=True):
                cumulative += d[count]
                if cumulative / total_bases_by_sample[s_name] > 0.01:
                    max_x = max(max_x, count)
                    break

        rates_within_threshs = dict()
        for s_name, hist in self.qualimap_bamqc_coverage_hist.items():
            total = total_bases_by_sample[s_name]
            # Make a range of depths that isn't stupidly huge for high coverage expts
            depth_range = list(range(0, max_x + 1, math.ceil(float(max_x) / 400.0) if max_x > 0 else 1))
            # Check that we have our specified coverages in the list
            for c in self.covs:
                if int(c) not in depth_range:
                    depth_range.append(int(c))
            # Calculate the coverage rates for this range of coverages
            rates_within_threshs[s_name] = _calculate_bases_within_thresholds(hist, total, depth_range)
            # Add requested coverage levels to the General Statistics table
            for c in self.covs:
                if int(c) in rates_within_threshs[s_name]:
                    self.general_stats_data[s_name]["{}_x_pc".format(c)] = rates_within_threshs[s_name][int(c)]
                else:
                    self.general_stats_data[s_name]["{}_x_pc".format(c)] = 0

        # Section 1 - BamQC Coverage Histogram
        self.add_section(
            name="Coverage histogram",
            anchor="qualimap-coverage-histogram",
            description="Distribution of the number of locations in the reference genome with a given depth of coverage.",
            helptext=coverage_histogram_helptext,
            plot=linegraph.plot(
                self.qualimap_bamqc_coverage_hist,
                {
                    "id": "qualimap_coverage_histogram",
                    "title": "Qualimap BamQC: Coverage histogram",
                    "ylab": "Genome bin counts",
                    "xlab": "Coverage (X)",
                    "ymin": 0,
                    "xmin": 0,
                    "xmax": max_x,
                    "xDecimals": False,
                    "tt_label": "<b>{point.x}X</b>: {point.y}",
                },
            ),
        )
        # Section 2 - BamQC cumulative coverage genome fraction
        self.add_section(
            name="Cumulative genome coverage",
            anchor="qualimap-cumulative-genome-fraction-coverage",
            description="Percentage of the reference genome with at least the given depth of coverage.",
            helptext=genome_fraction_helptext,
            plot=linegraph.plot(
                rates_within_threshs,
                {
                    "id": "qualimap_genome_fraction",
                    "title": "Qualimap BamQC: Genome fraction covered by at least X reads",
                    "ylab": "Fraction of reference (%)",
                    "xlab": "Coverage (X)",
                    "ymax": 100,
                    "ymin": 0,
                    "xmin": 0,
                    "xmax": max_x,
                    "xDecimals": False,
                    "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
                },
            ),
        )

    # Section 3 - Insert size histogram
    if len(self.qualimap_bamqc_insert_size_hist) > 0:
        insert_size_helptext = """
        To overcome limitations in the length of DNA or RNA sequencing reads,
        many sequencing instruments can produce two or more shorter reads from
        one longer fragment in which the relative position of reads is
        approximately known, such as paired-end or mate-pair reads
        (<a href="https://doi.org/10.1146/annurev-anchem-062012-092628"
        target="_blank">Mardis 2013</a>). Such techniques can extend the reach
        of sequencing technology, allowing for more accurate placement of reads
        (<a href="https://doi.org/10.1146/annurev-genom-090413-025358"
        target="_blank">Reinert et al. 2015</a>) and better resolution of repeat
        regions (<a href="https://doi.org/10.1146/annurev-genom-090413-025358"
        target="_blank">Reinert et al. 2015</a>), as well as detection of
        structural variation (<a href="https://doi.org/10.1038/nrg2958"
        target="_blank">Alkan et al. 2011</a>) and chimeric transcripts
        (<a href="https://doi.org/10.1073/pnas.0904720106"
        target="_blank">Maher et al. 2009</a>).

        All these methods assume that the approximate size of an insert is known.
        (Insert size can be defined as the length in bases of a sequenced DNA or
        RNA fragment, excluding technical sequences such as adapters, which are
        typically removed before alignment.) This plot allows for that assumption
        to be assessed. With the set of mapped fragments for a given sample, QualiMap
        groups the fragments by insert size, then plots the frequency of mapped
        fragments (y-axis) over a range of insert sizes (x-axis). In an ideal case,
        the distribution of fragment sizes for a sequencing library would culminate
        in a single peak indicating average insert size, with a narrow spread
        indicating highly consistent fragment lengths.

        QualiMap calculates insert sizes as follows: for each fragment in which
        every read mapped successfully to the same reference sequence, it
        extracts the insert size from the `TLEN` field of the leftmost read
        (see the <a href="http://qualimap.bioinfo.cipf.es/doc_html/index.html"
        target="_blank">Qualimap 2 documentation</a>), where the `TLEN` (or
        'observed Template LENgth') field contains 'the number of bases from the
        leftmost mapped base to the rightmost mapped base'
        (<a href="https://samtools.github.io/hts-specs/" target="_blank">SAM
        format specification</a>). Note that because it is defined in terms of
        alignment to a reference sequence, the value of the `TLEN` field may
        differ from the insert size due to factors such as alignment clipping,
        alignment errors, or structural variation or splicing in a gap between
        reads from the same fragment."""
        self.add_section(
            name="Insert size histogram",
            anchor="qualimap-insert-size-histogram",
            description="Distribution of estimated insert sizes of mapped reads.",
            helptext=insert_size_helptext,
            plot=linegraph.plot(
                self.qualimap_bamqc_insert_size_hist,
                {
                    "id": "qualimap_insert_size",
                    "title": "Qualimap BamQC: Insert size histogram",
                    "ylab": "Fraction of reads",
                    "xlab": "Insert Size (bp)",
                    "ymin": 0,
                    "xmin": 0,
                    "tt_label": "<b>{point.x} bp</b>: {point.y}",
                },
            ),
        )

    # Section 4 - GC-content distribution
    if len(self.qualimap_bamqc_gc_content_dist) > 0:
        gc_content_helptext = """
        GC bias is the difference between the guanine-cytosine content
        (GC-content) of a set of sequencing reads and the GC-content of the DNA
        or RNA in the original sample. It is a well-known issue with sequencing
        systems, and may be introduced by PCR amplification, among other factors
        (<a href="https://doi.org/10.1093/nar/gks001" target="_blank">Benjamini
        & Speed 2012</a>; <a href="https://doi.org/10.1186/gb-2013-14-5-r51"
        target="_blank">Ross et al. 2013</a>).

        QualiMap calculates the GC-content of individual mapped reads, then
        groups those reads by their GC-content (*1%, 2%, &#8230;, 100%*), and
        plots the frequency of mapped reads (y-axis) at each level of GC-content
        (x-axis). This plot shows the GC-content distribution of mapped reads
        for each read dataset, which should ideally resemble that of the
        original sample. It can be useful to display the GC-content distribution
        of an appropriate reference sequence for comparison, and QualiMap has an
        option to do this (see the <a href="http://qualimap.bioinfo.cipf.es/doc_html/index.html"
        target="_blank">Qualimap 2 documentation</a>)."""
        extra_series = []
        for i, (species_name, species_data) in enumerate(sorted(self.qualimap_bamqc_gc_by_species.items())):
            extra_series.append(
                {
                    "name": species_name,
                    "data": list(species_data.items()),
                    "dashStyle": "Dash",
                    "lineWidth": 1,
                    "color": ["#000000", "#E89191"][i % 2],
                }
            )
        if len(self.qualimap_bamqc_gc_content_dist) == 1:
            desc = "The solid line represents the distribution of GC content of mapped reads for the sample."
        else:
            desc = "Each solid line represents the distribution of GC content of mapped reads for a given sample."
        lg_config = {
            "id": "qualimap_gc_content",
            "title": "Qualimap BamQC: GC content distribution",
            "ylab": "Fraction of reads",
            "xlab": "GC content (%)",
            "ymin": 0,
            "xmin": 0,
            "xmax": 100,
            "tt_label": "<b>{point.x}%</b>: {point.y:.3f}",
        }
        if len(extra_series) == 1:
            desc += " The dotted line represents a pre-calculated GC distribution for the reference genome."
            lg_config["extra_series"] = extra_series
        elif len(extra_series) > 1:
            desc += " Each dotted line represents a pre-calculated GC distribution for a specific reference genome."
            lg_config["extra_series"] = extra_series

        self.add_section(
            name="GC content distribution",
            anchor="qualimap-gc-distribution",
            description=desc,
            helptext=gc_content_helptext,
            plot=linegraph.plot(self.qualimap_bamqc_gc_content_dist, lg_config),
        )


def general_stats_headers(self):
    try:
        hidecovs = config.qualimap_config["general_stats_coverage_hidden"]
        assert type(hidecovs) == list
        log.debug("Hiding Qualimap thresholds: {}".format(", ".join([i for i in hidecovs])))
    except (AttributeError, TypeError, KeyError, AssertionError):
        hidecovs = [1, 5, 10, 50]
    hidecovs = [str(i) for i in hidecovs]

    self.general_stats_headers["avg_gc"] = {
        "title": "% GC",
        "description": "Mean GC content",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "PuRd",
        "format": "{:,.0f}",
    }
    self.general_stats_headers["median_insert_size"] = {
        "title": "Ins. size",
        "description": "Median insert size",
        "min": 0,
        "scale": "PuOr",
        "format": "{:,.0f}",
    }
    for c in self.covs:
        self.general_stats_headers["{}_x_pc".format(c)] = {
            "title": "&ge; {}X".format(c),
            "description": "Fraction of genome with at least {}X coverage".format(c),
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
            "hidden": c in hidecovs,
        }
    self.general_stats_headers["median_coverage"] = {
        "title": "Median cov",
        "description": "Median coverage",
        "min": 0,
        "suffix": "X",
        "scale": "BuPu",
    }
    self.general_stats_headers["mean_coverage"] = {
        "title": "Mean cov",
        "description": "Mean coverage",
        "min": 0,
        "suffix": "X",
        "scale": "BuPu",
    }
    self.general_stats_headers["percentage_aligned_on_target"] = {
        "title": "% On target",
        "description": "% mapped reads on target region",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "YlGn",
    }
    self.general_stats_headers["regions_size"] = {
        "title": "{} Region size".format(config.read_count_prefix),
        "description": "Size of target region",
        "suffix": " bp",
        "scale": "PuBuGn",
        "hidden": True,
    }
    self.general_stats_headers["regions_mapped_reads"] = {
        "title": "{} On target".format(config.read_count_prefix),
        "description": "Number of mapped reads on target region ({})".format(config.read_count_desc),
        "scale": "RdYlGn",
        "shared_key": "read_count",
        "hidden": True,
    }
    self.general_stats_headers["general_error_rate"] = {
        "title": "Error rate",
        "description": "Alignment error rate. Total edit distance (SAM NM field) over the number of mapped bases",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "OrRd",
        "format": "{0:.2f}",
        "hidden": True,
    }
    self.general_stats_headers["percentage_aligned"] = {
        "title": "% Aligned",
        "description": "% mapped reads",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "YlGn",
    }
    self.general_stats_headers["mapped_reads"] = {
        "title": "{} Aligned".format(config.read_count_prefix),
        "description": "Number of mapped reads ({})".format(config.read_count_desc),
        "scale": "RdYlGn",
        "shared_key": "read_count",
        "hidden": True,
    }
    self.general_stats_headers["total_reads"] = {
        "title": "{} Total reads".format(config.read_count_prefix),
        "description": "Number of reads ({})".format(config.read_count_desc),
        "scale": "Blues",
        "shared_key": "read_count",
        "hidden": True,
    }


def _calculate_bases_within_thresholds(bases_by_depth, total_size, depth_thresholds):
    bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
    rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)

    dt = sorted(depth_thresholds, reverse=True)
    c = 0
    for depth in sorted(bases_by_depth.keys(), reverse=True):
        while depth < dt[c]:
            c += 1
            bases_within_threshs[dt[c]] = bases_within_threshs[dt[c - 1]]
        if depth >= dt[c]:
            bases_within_threshs[dt[c]] += bases_by_depth[depth]
    while c + 1 < len(dt):
        c += 1
        bases_within_threshs[dt[c]] = total_size
    for t in dt:
        bs = bases_within_threshs[t]
        if total_size > 0:
            rate = 100.0 * bases_within_threshs[t] / total_size
            assert rate <= 100, (
                "Error: rate is > 1: rate = " + str(rate) + ", bases = " + str(bs) + ", size = " + str(total_size)
            )
            rates_within_threshs[t] = rate
    return rates_within_threshs
