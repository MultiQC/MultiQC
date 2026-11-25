import logging
import math

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

# Reasons that leafcutter_ds gives for not testing a cluster, with the keys and labels used here.
# The strings come from differential_splicing() in leafcutter/R/differential_splicing.R and from
# leafcutter/differential_splicing/differential_splicing.py in leafcutter-ds.
SKIP_STATUSES: dict[str, tuple[str, str]] = {
    "Too many introns in cluster": ("too_many_introns", "Not tested: too many introns"),
    "<=1 junction in cluster": ("one_junction", "Not tested: one junction"),
    "<=1 sample with coverage>0": ("one_sample_covered", "Not tested: one sample with reads"),
    "<=1 sample with coverage>min_coverage": ("one_sample_min_coverage", "Not tested: one sample above min coverage"),
    "<2 introns used in >=min_samples_per_intron samples": (
        "too_few_introns_used",
        "Not tested: too few introns in enough samples",
    ),
    "Not enough valid samples": ("not_enough_samples", "Not tested: not enough valid samples"),
}

# The --FDR default of leafviz/prepare_results.R
DEFAULT_PADJ_CUTOFF = 0.05


def parse_leafcutter_cluster_significance(module: BaseMultiqcModule) -> int:
    """Parse `*_cluster_significance.txt` files from `leafcutter_ds`, one comparison per file."""
    padj_cutoff = float(getattr(config, "leafcutter", {}).get("padj_cutoff", DEFAULT_PADJ_CUTOFF))

    data: dict[str, dict[str, int]] = {}
    for f in module.find_log_files("leafcutter/cluster_significance", filehandles=True):
        try:
            parsed = parse_cluster_significance(f["f"], padj_cutoff)
        except ValueError as e:
            log.warning(f"Skipping {f['fn']}: {e}")
            continue
        s_name = f["s_name"]
        if s_name in data:
            log.debug(f"Duplicate comparison name found in {f['fn']}! Overwriting: {s_name}")
        data[s_name] = parsed
        module.add_data_source(f, s_name=s_name, section="cluster_significance")

    data = module.ignore_samples(data)
    if len(data) == 0:
        return 0

    module.add_software_version(None)

    keys = {
        "significant": {"name": f"Significant (adjusted p < {padj_cutoff:g})", "color": "#7b4ea3"},
        "not_significant": {"name": "Not significant", "color": "#437bb1"},
        "no_pvalue": {"name": "Tested, no adjusted p-value", "color": "#a6c4e2"},
    }
    not_tested_colours = ["#e0e0e0", "#cccccc", "#b1b1b1", "#969696", "#7a7a7a", "#5e5e5e", "#424242", "#262626"]
    not_tested = [*SKIP_STATUSES.values(), ("timeout", "Not tested: timeout"), ("error", "Not tested: error")]
    for (key, label), colour in zip(not_tested, not_tested_colours):
        keys[key] = {"name": label, "color": colour}

    module.add_section(
        name="Cluster Test Outcomes",
        anchor="leafcutter_cluster_test_outcomes",
        description=f"""
            Outcome of the differential splicing test for each intron cluster, one bar per
            comparison. Significant: adjusted p < {padj_cutoff:g}.
            """,
        helptext=f"""
            `leafcutter_ds` only tests clusters that pass its coverage filters. The reasons for
            not testing a cluster, as LeafCutter reports them in the `status` column, are:

            * **Too many introns**: `Too many introns in cluster` (more than `--max_cluster_size`)
            * **One junction**: `<=1 junction in cluster`
            * **One sample with reads**: `<=1 sample with coverage>0`
            * **One sample above min coverage**: `<=1 sample with coverage>min_coverage`
            * **Too few introns in enough samples**: `<2 introns used in >=min_samples_per_intron samples`
            * **Not enough valid samples**: fewer than `--min_samples_per_group` samples in each
              group reach `--min_coverage` reads
            * **Timeout**: the model fit took longer than `--timeout` (R version only)
            * **Error**: the model fit failed, and LeafCutter reported the error message

            A comparison with few tested clusters usually has too few samples per group or too
            little depth for the filters, rather than little differential splicing.

            The p-values are adjusted with the Benjamini-Hochberg method. The default threshold
            is the FDR threshold of LeafCutter's `prepare_results.R`. Set it in a MultiQC config
            file (default shown):

            ```yaml
            leafcutter:
              padj_cutoff: {DEFAULT_PADJ_CUTOFF:g}
            ```
            """,
        plot=bargraph.plot(
            data,
            keys,
            {
                "id": "leafcutter_cluster_test_outcomes_plot",
                "title": "LeafCutter: Cluster Test Outcomes",
                "ylab": "Clusters",
                "cpswitch_counts_label": "Number of clusters",
            },
        ),
    )

    module.write_data_file(data, "multiqc_leafcutter_cluster_significance")

    return len(data)


def parse_cluster_significance(f, padj_cutoff: float) -> dict[str, int]:
    """
    Count the clusters in a `leafcutter_ds` cluster significance table by test outcome.

    The R version writes `status loglr df p cluster p.adjust` without an exon file and moves
    `cluster` to the front with a trailing `genes` column when given one. leafcutter-ds always
    writes `cluster status loglr df p p.adjust genes`, sometimes with an `annotations` column.
    The columns are therefore read by name.
    """
    header = f.readline().rstrip("\r\n").split("\t")
    try:
        status_col = header.index("status")
        padj_col = header.index("p.adjust")
    except ValueError:
        raise ValueError("missing 'status' or 'p.adjust' column in the header")

    counts = {key: 0 for key in ["significant", "not_significant", "no_pvalue", "timeout", "error"]}
    for key, _ in SKIP_STATUSES.values():
        counts[key] = 0

    for line_number, line in enumerate(f, start=2):
        fields = line.rstrip("\r\n").split("\t")
        if fields == [""]:
            continue
        if len(fields) != len(header):
            raise ValueError(f"expected {len(header)} fields on line {line_number}, found {len(fields)}")
        status = fields[status_col]
        if status == "Success":
            padj = float(fields[padj_col]) if fields[padj_col] != "NA" else math.nan
            if math.isnan(padj):
                counts["no_pvalue"] += 1
            elif padj < padj_cutoff:
                counts["significant"] += 1
            else:
                counts["not_significant"] += 1
        elif status in SKIP_STATUSES:
            counts[SKIP_STATUSES[status][0]] += 1
        elif status == "timeout":
            counts["timeout"] += 1
        else:
            counts["error"] += 1

    counts["clusters"] = sum(counts.values())
    counts["tested"] = counts["significant"] + counts["not_significant"] + counts["no_pvalue"]
    return counts
