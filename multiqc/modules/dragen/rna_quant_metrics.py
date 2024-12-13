import logging

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class DragenRnaQuantMetrics(BaseMultiqcModule):
    def add_rna_metrics(self):
        data_by_sample = dict()
        has_old_style = False
        has_new_style = False

        for f in self.find_log_files("dragen/rna_quant_metrics"):
            data = parse_metrics_file(f)
            s_name = f["s_name"]

            # Detect style and normalize data if needed
            if "Forward transcript fragments" in data:
                has_new_style = True
                data["Transcript fragments"] = data.get("Forward transcript fragments", 0) + data.get(
                    "Reverse transcript fragments", 0
                )
            elif "Transcript fragments" in data:
                has_old_style = True

            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="rna_quant_metrics")
            data_by_sample[s_name] = data

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        # Only plot data, don't want to write this to a file
        # (can do so with --export-plots already)
        # self.write_data_file(data_by_sample, "dragen_quant_metrics")

        # Warn if mixed styles are detected
        if has_old_style and has_new_style:
            log.warning(
                "Mixed old-style (Transcript fragments) and new-style (Forward/Reverse transcript fragments) "
                "metrics detected. Forward and Reverse fragments will be summed up to match old-style format, "
                "and the old style will be used for all samples."
            )

        if has_old_style:
            # Always use old-style categories for consistency when mixing
            categories = ["Transcript fragments"]
        else:
            categories = ["Reverse transcript fragments", "Forward transcript fragments"]

        categories += [
            "Strand mismatched fragments",
            "Ambiguous strand fragments",
            "Unknown transcript fragments",
            "Intron fragments",
            "Intergenic fragments",
        ]

        self.add_section(
            name="RNA Quantification Metrics",
            anchor="dragen-rna-quant-metrics",
            description="RNA quantification metrics for DRAGEN.",
            plot=bargraph.plot(
                [
                    {
                        sample: {metric: stat for metric, stat in data.items() if metric in categories}
                        for sample, data in data_by_sample.items()
                    }
                ],
                cats=categories,
                pconfig={
                    "id": "dragen_rna_quant_metrics",
                    "title": "Dragen: RNA Quant Metrics",
                    "ylab": "Fragments",
                    "cpswitch_counts_label": "Fragments",
                    "data_labels": [
                        {
                            "name": "Transcript / Intronic / Intergenic Fragments",
                            "ylab": "Fragments",
                            "cpswitch_counts_label": "Fragments",
                        },
                        {
                            "name": "Fragment Orientations",
                            "ylab": "Fragments",
                            "cpswitch_counts_label": "Fragments",
                        },
                    ],
                },
            ),
        )

        return data_by_sample.keys()


def parse_metrics_file(f):
    """
    sample.time_metrics.csv

    RNA QUANTIFICATION STATISTICS,,Library orientation,ISR
    RNA QUANTIFICATION STATISTICS,,Total Genes,58735
    RNA QUANTIFICATION STATISTICS,,Total Transcripts,206601
    RNA QUANTIFICATION STATISTICS,,Coding Genes,28815
    RNA QUANTIFICATION STATISTICS,,Median transcript CV coverage,0.61
    RNA QUANTIFICATION STATISTICS,,Median 5' coverage bias,0.1868
    RNA QUANTIFICATION STATISTICS,,Median 3' coverage bias,0.1228
    RNA QUANTIFICATION STATISTICS,,Number of genes with coverage > 1x,14014,23.86
    RNA QUANTIFICATION STATISTICS,,Number of genes with coverage > 10x,10084,17.17
    RNA QUANTIFICATION STATISTICS,,Number of genes with coverage > 30x,7769,13.23
    RNA QUANTIFICATION STATISTICS,,Number of genes with coverage > 100x,4603,7.84
    RNA QUANTIFICATION STATISTICS,,Forward transcript fragments,0,0.00
    RNA QUANTIFICATION STATISTICS,,Reverse transcript fragments,27628845,87.03
    RNA QUANTIFICATION STATISTICS,,Strand mismatched fragments,908889,2.86
    RNA QUANTIFICATION STATISTICS,,Ambiguous strand fragments,0,0.00
    RNA QUANTIFICATION STATISTICS,,Unknown transcript fragments,1613938,5.08
    RNA QUANTIFICATION STATISTICS,,Intron fragments,550923,1.74
    RNA QUANTIFICATION STATISTICS,,Intergenic fragments,1045169,3.29
    RNA QUANTIFICATION STATISTICS,,Fold coverage of all exons,89.60
    RNA QUANTIFICATION STATISTICS,,Fold coverage of introns,0.22
    RNA QUANTIFICATION STATISTICS,,Fold coverage of intergenic regions,0.14
    RNA QUANTIFICATION STATISTICS,,Fold coverage of coding exons,111.21
    """
    data = {}
    for line in f["f"].splitlines():
        tokens = line.split(",")
        if len(tokens) == 4:
            analysis, _, metric, stat = tokens
            percentage = None
        elif len(tokens) == 5:
            analysis, _, metric, stat, percentage = tokens
        else:
            raise ValueError(f"Unexpected number of tokens in line {line}")

        try:
            stat = float(stat)
        except ValueError:
            pass
        data[metric] = stat

    return data
