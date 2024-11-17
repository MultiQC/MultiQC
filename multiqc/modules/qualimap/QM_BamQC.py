import logging
import os
from typing import Union, Dict, Tuple, List

import math
import re

from multiqc import config, BaseMultiqcModule
from multiqc.modules.qualimap import parse_numerals, get_s_name
from multiqc.plots import linegraph
from multiqc.utils.util_functions import update_dict

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule):
    """Find Qualimap BamQC reports and parse their data"""

    all_general_stats: Dict[str, Dict] = dict()

    # General stats - genome_results.txt
    genome_results, general_stats = parse_genome_results(module)
    if genome_results:
        module.write_data_file(genome_results, "multiqc_qualimap_bamqc_genome_results")
    update_dict(all_general_stats, general_stats)

    # Coverage - coverage_histogram.txt
    coverage_hist, general_stats = parse_coverage(module)
    update_dict(all_general_stats, general_stats)

    # Insert size - insert_size_histogram.txt
    insert_size_hist, general_stats = parse_insert_size(module)
    update_dict(all_general_stats, general_stats)

    # GC distribution - mapped_reads_gc-content_distribution.txt
    gc_content_dist, gc_by_species, general_stats = parse_gc_dist(module)
    update_dict(all_general_stats, general_stats)

    num_parsed = max(
        len(genome_results),
        len(coverage_hist),
        len(insert_size_hist),
        len(gc_content_dist),
    )
    # Go no further if nothing found
    if num_parsed == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    threshs, hidden_threshs = config.get_cov_thresholds("qualimap_config")

    # Make the plots for the report
    general_stats = report_sections(
        module,
        threshs,
        coverage_hist,
        insert_size_hist,
        gc_content_dist,
        gc_by_species,
    )
    update_dict(all_general_stats, general_stats)

    # Set up the general stats table
    module.general_stats_addcols(
        all_general_stats,
        general_stats_headers(threshs, hidden_threshs),
        namespace="BamQC",
    )

    # Return the number of reports we found
    return num_parsed


def parse_genome_results(module: BaseMultiqcModule) -> Tuple[Dict, Dict]:
    """Parse the contents of the Qualimap BamQC genome_results.txt file"""
    metrics_by_sample: Dict[str, Dict[str, Union[int, float, str]]] = dict()
    general_stats: Dict = dict()

    for f in module.find_log_files("qualimap/bamqc/genome_results"):
        int_metrics = {
            "number of reads": "total_reads",
            "number of mapped reads": "mapped_reads",
            "number of mapped reads inside": "regions_mapped_reads",
            "number of mapped bases": "mapped_bases",
            "number of sequenced bases": "sequenced_bases",
            "regions size": "regions_size",
        }
        float_metrics = {
            "mean insert size": "mean_insert_size",
            "median insert size": "median_insert_size",
            "mean mapping quality": "mean_mapping_quality",
            "mean coverageData": "mean_coverage",
        }
        rate_metrics = {  # useful to determine decimal separator
            "general error rate": "general_error_rate",
            "GC percentage": "gc_percentage",
            "duplication rate": "duplication_rate",
            "homopolymer indels": "homopolymer_indels",
        }

        value_regex = re.compile(r"\s+[\d,\.\xa0]+\s+")
        preparsed_d = dict()
        # Keeping track of the section because "number of mapped reads" occurs both
        # under "Globals" and "Globals inside"
        section = None
        for line in f["f"].splitlines():
            if line.startswith(">>>>>>>"):
                section = line[8:]
            elif not section:
                continue
            if "=" in line:
                key, val = line.split("=", 1)
                m = re.search(value_regex, val)
                if m:
                    val = m.group(0)
                    val = re.sub(r"\xa0", "", val)
                key = key.strip()
                if "Globals inside" in section and key == "number of mapped reads":
                    key = "number of mapped reads inside"
                preparsed_d[key] = val.strip()

        # Check we have an input filename
        if "bam file" not in preparsed_d:
            log.debug(f"Couldn't find an input filename in genome_results file {f['fn']}")
            continue

        s_name = module.clean_s_name(preparsed_d["bam file"], f)
        if s_name in metrics_by_sample:
            log.debug(f"Duplicate genome results sample name found! Overwriting: {s_name}")

        module.add_data_source(f, s_name=s_name, section="genome_results")

        d: Dict[str, Union[int, float, str]] = parse_numerals(
            preparsed_d,
            float_metrics=float_metrics,
            int_metrics=int_metrics,
            rate_metrics=rate_metrics,
            fpath=os.path.join(f["root"], f["fn"]),
        )

        if "general_error_rate" in d and isinstance(d["general_error_rate"], float):
            d["general_error_rate"] = d["general_error_rate"] * 100.0

        if (
            "mapped_reads" in d
            and "total_reads" in d
            and isinstance(d["total_reads"], int)
            and d["total_reads"] > 0
            and isinstance(d["mapped_reads"], int)
        ):
            d["percentage_aligned"] = (d["mapped_reads"] / d["total_reads"]) * 100.0

        if (
            "regions_mapped_reads" in d
            and "mapped_reads" in d
            and isinstance(d["mapped_reads"], int)
            and d["mapped_reads"] > 0
            and isinstance(d["regions_mapped_reads"], int)
        ):
            d["percentage_aligned_on_target"] = (d["regions_mapped_reads"] / d["mapped_reads"]) * 100.0

        # Add to general stats table & calculate a nice % aligned
        general_stats[s_name] = {
            k: v
            for k, v in d.items()
            if k
            in [
                "total_reads",
                "mapped_reads",
                "general_error_rate",
                "mean_coverage",
                "regions_size",
                "regions_mapped_reads",
                "percentage_aligned",
                "percentage_aligned_on_target",
            ]
        }

        metrics_by_sample[s_name] = d

    metrics_by_sample = module.ignore_samples(metrics_by_sample)
    general_stats = module.ignore_samples(general_stats)
    return metrics_by_sample, general_stats


def parse_coverage(module: BaseMultiqcModule) -> Tuple[Dict, Dict]:
    """Parse the contents of the Qualimap BamQC Coverage Histogram file"""
    coverage_hist: Dict = dict()
    general_stats: Dict = dict()

    for f in module.find_log_files("qualimap/bamqc/coverage", filehandles=True):
        # Get the sample name from the parent directory
        # Typical path: <sample name>/raw_data_qualimapReport/coverage_histogram.txt
        s_name = get_s_name(module, f)
        if s_name in coverage_hist:
            log.debug(f"Duplicate coverage histogram sample name found! Overwriting: {s_name}")

        module.add_data_source(f, s_name=s_name, section="coverage_histogram")

        d = dict()
        for line in f["f"]:
            if line.startswith("#"):
                continue
            coverage, count = line.split(None, 1)
            coverage = int(round(float(coverage.replace(",", "."))))
            count = float(count)
            d[coverage] = count

        if len(d) == 0:
            log.debug(f"Couldn't parse contents of coverage histogram file {f['fn']}")
            continue

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

        coverage_hist[s_name] = d
        general_stats[s_name] = {"median_coverage": median_coverage}

    coverage_hist = module.ignore_samples(coverage_hist)
    general_stats = module.ignore_samples(general_stats)
    return coverage_hist, general_stats


def parse_insert_size(module) -> Tuple[Dict, Dict]:
    """Parse the contents of the Qualimap BamQC Insert Size Histogram file"""
    data_by_sample: Dict = dict()
    general_stats: Dict = dict()

    for f in module.find_log_files("qualimap/bamqc/insert_size", filehandles=True):
        # Get the sample name from the parent directory
        # Typical path: <sample name>/raw_data_qualimapReport/insert_size_histogram.txt
        s_name = get_s_name(module, f)
        if s_name in data_by_sample:
            log.debug(f"Duplicate insert size histogram sample name found! Overwriting: {s_name}")

        module.add_data_source(f, s_name=s_name, section="insert_size_histogram")

        d = dict()
        for line in f["f"]:
            if line.startswith("#"):
                continue
            insertsize, count = line.split(None, 1)
            insertsize = int(round(float(insertsize)))
            count = float(count) / 1000000
            if insertsize != 0:
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
        general_stats[s_name] = {"median_insert_size": median_insert_size}

        # Save results
        data_by_sample[s_name] = d

    data_by_sample = module.ignore_samples(data_by_sample)
    general_stats = module.ignore_samples(general_stats)
    return data_by_sample, general_stats


def parse_gc_dist(module) -> Tuple[Dict, Dict, Dict]:
    """Parse the contents of the Qualimap BamQC Mapped Reads GC content distribution file"""
    gc_content_dist: Dict = dict()
    gc_by_species: Dict = dict()
    general_stats: Dict = dict()

    for f in module.find_log_files("qualimap/bamqc/gc_dist", filehandles=True):
        # Get the sample name from the parent directory
        # Typical path: <sample name>/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt
        s_name = get_s_name(module, f)
        if s_name in gc_content_dist:
            log.debug(f"Duplicate Mapped Reads GC content distribution sample name found! Overwriting: {s_name}")

        module.add_data_source(f, s_name=s_name, section="mapped_gc_distribution")

        d = dict()
        reference_species = None
        reference_d = dict()
        avg_gc = 0.0
        for line in f["f"]:
            if line.startswith("#"):
                sections = line.strip("\n").split("\t", 3)
                if len(sections) > 2:
                    reference_species = sections[2]
                continue
            sections = line.strip("\n").split("\t", 3)
            gc = int(round(float(sections[0])))
            content = float(sections[1])
            avg_gc += gc * content
            d[gc] = content
            if len(sections) > 2:
                reference_content = float(sections[2])
                reference_d[gc] = reference_content

        # Add average GC to the general stats table
        general_stats[s_name] = {"avg_gc": avg_gc}

        # Save results
        gc_content_dist[s_name] = d
        if reference_species and reference_species not in gc_by_species:
            gc_by_species[reference_species] = reference_d

    gc_content_dist = module.ignore_samples(gc_content_dist)
    gc_by_species = module.ignore_samples(gc_by_species)
    general_stats = module.ignore_samples(general_stats)
    return gc_content_dist, gc_by_species, general_stats


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


def report_sections(
    module: BaseMultiqcModule,
    threshs: List[int],
    coverage_hist: Dict,
    insert_size_hist: Dict,
    gc_content_dist: Dict,
    gc_by_species: Dict,
) -> Dict:
    """Add results from Qualimap BamQC parsing to the report"""
    # Append to self.sections list
    general_stats: Dict = {s_name: {} for s_name in coverage_hist}

    if len(coverage_hist) > 0:
        # Chew back on histogram to prevent long flat tail
        # (find a sensible max x - lose 1% of longest tail)
        max_x = 20
        total_bases_by_sample = dict()
        for s_name, d in coverage_hist.items():
            total_bases_by_sample[s_name] = sum(d.values())
            cumulative = 0
            for count in sorted(d.keys(), reverse=True):
                cumulative += d[count]
                if cumulative / total_bases_by_sample[s_name] > 0.01:
                    max_x = max(max_x, count)
                    break

        rates_within_threshs = dict()
        for s_name, hist in coverage_hist.items():
            total = total_bases_by_sample[s_name]
            # Make a range of depths that isn't stupidly huge for high coverage expts
            depth_range = list(range(0, max_x + 1, math.ceil(float(max_x) / 400.0) if max_x > 0 else 1))
            # Check that we have our specified coverages in the list
            for c in threshs:
                if int(c) not in depth_range:
                    depth_range.append(int(c))
            # Calculate the coverage rates for this range of coverages
            rates_within_threshs[s_name] = _calculate_bases_within_thresholds(hist, total, depth_range)
            # Add requested coverage levels to the General Statistics table
            for c in threshs:
                if int(c) in rates_within_threshs[s_name]:
                    general_stats[s_name][f"{c}_x_pc"] = rates_within_threshs[s_name][int(c)]
                else:
                    general_stats[s_name][f"{c}_x_pc"] = 0

        # Section 1 - BamQC Coverage Histogram
        module.add_section(
            name="Coverage histogram",
            anchor="qualimap-coverage-histogram",
            description="Distribution of the number of locations in the reference genome with a given depth of coverage.",
            helptext=coverage_histogram_helptext,
            plot=linegraph.plot(
                coverage_hist,
                {
                    "id": "qualimap_coverage_histogram",
                    "title": "QualiMap: BamQC: Coverage histogram",
                    "ylab": "Genome bin counts",
                    "xlab": "Coverage (X)",
                    "ymin": 0,
                    "xmin": 0,
                    "xmax": max_x,
                    "x_decimals": False,
                    "tt_label": "<b>{point.x}X</b>: {point.y}",
                },
            ),
        )
        # Section 2 - BamQC cumulative coverage genome fraction
        module.add_section(
            name="Cumulative genome coverage",
            anchor="qualimap-cumulative-genome-fraction-coverage",
            description="Percentage of the reference genome with at least the given depth of coverage.",
            helptext=genome_fraction_helptext,
            plot=linegraph.plot(
                rates_within_threshs,
                {
                    "id": "qualimap_genome_fraction",
                    "title": "QualiMap: BamQC: Genome fraction covered by at least X reads",
                    "ylab": "Fraction of reference (%)",
                    "xlab": "Coverage (X)",
                    "ymax": 100,
                    "ymin": 0,
                    "xmin": 0,
                    "xmax": max_x,
                    "x_decimals": False,
                    "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
                },
            ),
        )

    # Section 3 - Insert size histogram
    if len(insert_size_hist) > 0:
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
        module.add_section(
            name="Insert size histogram",
            anchor="qualimap-insert-size-histogram",
            description="Distribution of estimated insert sizes of mapped reads.",
            helptext=insert_size_helptext,
            plot=linegraph.plot(
                insert_size_hist,
                {
                    "id": "qualimap_insert_size",
                    "title": "QualiMap: BamQC: Insert size histogram",
                    "ylab": "Fraction of reads",
                    "xlab": "Insert Size (bp)",
                    "ymin": 0,
                    "xmin": 0,
                    "tt_label": "<b>{point.x} bp</b>: {point.y}",
                },
            ),
        )

    # Section 4 - GC-content distribution
    if len(gc_content_dist) > 0:
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
        for i, (species_name, species_data) in enumerate(sorted(gc_by_species.items())):
            extra_series.append(
                {
                    "name": species_name,
                    "pairs": list(species_data.items()),
                    "width": 1,
                    "dash": "dash",
                    "color": ["#000000", "#E89191"][i % 2],
                }
            )
        if len(gc_content_dist) == 1:
            desc = "The solid line represents the distribution of GC content of mapped reads for the sample."
        else:
            desc = "Each solid line represents the distribution of GC content of mapped reads for a given sample."
        lg_config = {
            "id": "qualimap_gc_content",
            "title": "QualiMap: BamQC: GC content distribution",
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

        module.add_section(
            name="GC content distribution",
            anchor="qualimap-gc-distribution",
            description=desc,
            helptext=gc_content_helptext,
            plot=linegraph.plot(gc_content_dist, lg_config),
        )

    return general_stats


def general_stats_headers(threshs, hidden_threshs) -> Dict:
    headers = dict()
    headers["avg_gc"] = {
        "title": "% GC",
        "description": "Mean GC content",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "PuRd",
        "format": "{:,.0f}",
    }
    headers["median_insert_size"] = {
        "title": "Ins. size",
        "description": "Median insert size",
        "min": 0,
        "scale": "PuOr",
        "format": "{:,.0f}",
    }
    for c in threshs:
        headers[f"{c}_x_pc"] = {
            "title": f"≥ {c}X",
            "description": f"Fraction of genome with at least {c}X coverage",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
            "hidden": c in hidden_threshs,
        }
    headers["median_coverage"] = {
        "title": "Median cov",
        "description": "Median coverage",
        "min": 0,
        "suffix": "X",
        "scale": "BuPu",
    }
    headers["mean_coverage"] = {
        "title": "Mean cov",
        "description": "Mean coverage",
        "min": 0,
        "suffix": "X",
        "scale": "BuPu",
    }
    headers["percentage_aligned_on_target"] = {
        "title": "% On target",
        "description": "% mapped reads on target region",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "YlGn",
    }
    headers["regions_size"] = {
        "title": f"{config.read_count_prefix} Region size",
        "description": "Size of target region",
        "suffix": " bp",
        "scale": "PuBuGn",
        "hidden": True,
    }
    headers["regions_mapped_reads"] = {
        "title": f"{config.read_count_prefix} On target",
        "description": f"Number of mapped reads on target region ({config.read_count_desc})",
        "scale": "RdYlGn",
        "shared_key": "read_count",
        "hidden": True,
    }
    headers["general_error_rate"] = {
        "title": "Error rate",
        "description": "Alignment error rate. Total edit distance (SAM NM field) over the number of mapped bases",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "OrRd",
        "format": "{0:.2f}",
        "hidden": True,
    }
    headers["percentage_aligned"] = {
        "title": "% Aligned",
        "description": "% mapped reads",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "scale": "YlGn",
    }
    headers["mapped_reads"] = {
        "title": f"{config.read_count_prefix} Aligned",
        "description": f"Number of mapped reads ({config.read_count_desc})",
        "scale": "RdYlGn",
        "shared_key": "read_count",
        "hidden": True,
    }
    headers["total_reads"] = {
        "title": f"{config.read_count_prefix} Total reads",
        "description": f"Number of reads ({config.read_count_desc})",
        "scale": "Blues",
        "shared_key": "read_count",
        "hidden": True,
    }
    return headers


def _calculate_bases_within_thresholds(bases_by_depth, total_size, depth_thresholds):
    bases_within_threshs = {depth: 0 for depth in depth_thresholds}
    rates_within_threshs = {depth: None for depth in depth_thresholds}

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
