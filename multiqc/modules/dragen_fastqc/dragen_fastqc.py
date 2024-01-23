import logging
import os

from ..base_module import ModuleNoSamplesFound
from .base_metrics import DragenBaseMetrics
from .content_metrics import DragenContentMetrics
from .gc_metrics import DragenFastqcGcMetrics
from .read_metrics import DragenReadMetrics
from .util import parse_fastqc_metrics_file

log = logging.getLogger(__name__)


class MultiqcModule(DragenBaseMetrics, DragenReadMetrics, DragenFastqcGcMetrics, DragenContentMetrics):
    """DRAGEN provides a number of different pipelines and outputs, including base calling, DNA and RNA alignment,
    post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
    However, it can be treated as a fast aligner with additional features on top, as users will unlikely use any
    features without enabling DRAGEN mapping. So we will treat this module as an alignment tool module and
    place it accordingly in the module_order list, in docs, etc.

    The QC metrics DRAGEN generates resemble those of samtools-stats, qualimap, mosdepth, bcftools-stats and alike.
    Whenver possible, the visual output is made similar to those modules.

    Note that this MultiQC module supports some of DRAGEN output but not all. Contributions are welcome!

    The code is structured in a way so every mix-in parses one type of QC file that DRAGEN generates
    (e.g. *.mapping_metrics.csv, *.wgs_fine_hist_normal.csv, etc). If a corresponding file is found, a mix-in adds
    a section into the report.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="DRAGEN-FastQC",
            anchor="dragen-fastqc",
            href="https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html",
            info=(
                " is a Bio-IT Platform that provides ultra-rapid secondary analysis of sequencing data"
                " using field-programmable gate array technology (FPGA)."
            ),
            # Can't find a DOI // doi=
        )

        self.css = {
            "assets/css/multiqc_fastqc.css": os.path.join(
                os.path.dirname(__file__), "..", "fastqc", "assets", "css", "multiqc_fastqc.css"
            )
        }
        self.js = {
            "assets/js/multiqc_dragen_fastqc.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "multiqc_dragen_fastqc.js"
            )
        }
        self.intro += '<script type="application/json" class="fastqc_passfails">["dragen_fastqc", {"per_base_sequence_content": {"TEST": "pass"}}]</script>'

        data_by_sample = {}
        for f in self.find_log_files("dragen_fastqc"):
            s_name, data_by_mate = parse_fastqc_metrics_file(f)

            # Clean sample name and key
            new_s_name = self.clean_s_name(s_name, f)
            data_by_mate[new_s_name] = data_by_mate.pop(s_name)
            s_name = new_s_name
            del new_s_name

            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.add_data_source(f, s_name)
            data_by_sample.update(data_by_mate)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names:
        self.dragen_fastqc_data = self.ignore_samples(data_by_sample)

        # TODO: Split this up and write the interesting bits to files
        # Currently is 13760 lines in my test, which is too much.
        # self.write_data_file(self.dragen_fastqc_data, "dragen_fastqc")

        samples_found = set()

        # "POSITIONAL QUALITY" and "POSITIONAL BASE MEAN QUALITY" metrics
        samples_found |= self.add_base_metrics()

        # "READ MEAN QUALITY" and "READ LENGTHS" metrics
        samples_found |= self.add_read_metrics()

        # "GC CONTENT" and "GC CONTENT QUALITY" metrics
        samples_found |= self.add_gc_metrics()

        # "POSITIONAL BASE CONTENT" metrics
        samples_found |= self.add_content_metrics()

        if len(samples_found) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(samples_found)} reports")
