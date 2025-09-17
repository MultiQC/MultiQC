import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .stats import parse_seqfu_stats

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Supported commands:

    - `stats`:

    ### seqfu stats

    #### Input files

    `seqfu stats` can generated reports in multiple formats, see https://telatin.github.io/seqfu2/tools/stats.html. Only TSVs with headers (default `seqfu stats` output) are currently detected and parsed by MultiQC.

    :::note
    `seqfu stats` has a `--multiqc` option that generates a `_mqc.txt` file can be used with MuliQC as custom content. This is different from this module which enables additional features.
    :::

    #### Configuration

    Sample names are automatically extracted from the "File" columns by default. If you only have one sample per file and prefer to use the filename as the sample name instead, you can set the global `use_filename_as_sample_name` option to `true` or list `seqfu` under it.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Seqfu",
            anchor="seqfu",
            target="seqfu",
            href="https://telatin.github.io/seqfu2",
            info="Manipulate FASTA/FASTQ files.",
            doi="10.3390/bioengineering8050059",
        )

        n = dict()

        # Call submodule functions
        n["stats"] = parse_seqfu_stats(self)
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} seqfu stats reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
