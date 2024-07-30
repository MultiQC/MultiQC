import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .analyze_saturation_mutagenesis import AnalyzeSaturationMutagenesisMixin
from .base_recalibrator import BaseRecalibratorMixin
from .varianteval import VariantEvalMixin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, AnalyzeSaturationMutagenesisMixin, BaseRecalibratorMixin, VariantEvalMixin):
    """
    Supported tools:

    - `AnalyzeSaturationMutagenesis`
    - `BaseRecalibrator`
    - `VariantEval`

    #### AnalyzeSaturationMutagenesis

    [AnalyzeSaturationMutagenesis](https://gatk.broadinstitute.org/hc/en-us/articles/4404604903451-AnalyzeSaturationMutagenesis-BETA-)
    is a (beta!) tool for counting variants in saturation mutagenesis experiments. It accepts mapped reads and a reference sequence and outputs
    a number of files for further analysis.

    #### BaseRecalibrator

    [BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)
    is a tool for detecting systematic errors in read base quality scores of aligned high-throughput
    sequencing reads. It outputs a base quality score recalibration table that can be used in
    conjunction with the
    [PrintReads](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php)
    tool to recalibrate base quality scores.

    #### VariantEval

    [VariantEval](https://software.broadinstitute.org/gatk/gatkdocs/current/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php)
    is a general-purpose tool for variant evaluation. It gives information about percentage of
    variants in dbSNP, genotype concordance, Ti/Tv ratios and a lot more.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="GATK",
            anchor="gatk",
            target="GATK",
            href="https://www.broadinstitute.org/gatk/",
            info="Wide variety of tools with a primary focus on variant discovery and genotyping.",
            doi=["10.1101/201178", "10.1002/0471250953.bi1110s43", "10.1038/ng.806", "10.1101/gr.107524.110"],
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = {}
        self.general_stats_data = {}

        # Call submodule functions
        n_reports_found = 0
        n_reports_found += self.parse_gatk_analyze_saturation_mutagenesis()
        n_reports_found += self.parse_gatk_base_recalibrator()
        n_reports_found += self.parse_gatk_varianteval()

        # Exit if we didn't find anything
        if n_reports_found == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    def parse_report(self, lines, table_names):
        """Parse a GATK report https://software.broadinstitute.org/gatk/documentation/article.php?id=1244

        Only GATTable entries are parsed.  Tables are returned as a dict of tables.
        Each table is a dict of arrays, where names correspond to column names, and arrays
        correspond to column values.

        Args:
            lines (file handle): an iterable over the lines of a GATK report.
            table_names (dict): a dict with keys that are GATK report table names
                (e.g. "#:GATKTable:Quantized:Quality quantization map"), and values that are the
                keys in the returned dict.

        Returns:
            {
                table_1:
                    {
                        col_1: [ val_1, val_2, ... ]
                        col_2: [ val_1, val_2, ... ]
                        ...
                    }
                table_2:
                    ...
            }
        """

        report = dict()
        lines = (line for line in lines)
        for line in lines:
            line = line.rstrip()
            if line in table_names.keys():
                report[table_names[line]] = self.parse_gatk_report_table(lines)
        return report

    @staticmethod
    def parse_gatk_report_table(lines):
        headers = next(lines).rstrip().split()
        table = {h: [] for h in headers}
        for line in lines:
            line = line.rstrip()

            # testing to see if we have reached the end of a table in a GATKReport
            if line == "":
                break

            for index, value in enumerate(line.split()):
                table[headers[index]].append(value)
        return table
