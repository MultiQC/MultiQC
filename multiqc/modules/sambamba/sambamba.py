import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from multiqc.modules.sambamba.markdup import parse_sambamba_markdup

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Supported commands:

    - `markdup`

    #### markdup

    This module parses key phrases in the output log files to find duplicate +
    unique reads and then calculates duplicate rate per sample. It will work for both
    single and paired-end data. The absolute number of reads by type are displayed
    in a stacked bar plot, and duplicate rates are in the general statistics table.

    Duplicate rates are calculated as follows:

    #### Paired end

    > `duplicate_rate = duplicateReads / (sortedEndPairs * 2 + singleEnds - singleUnmatchedPairs) * 100`

    #### Single end

    > `duplicate_rate = duplicateReads / singleEnds * 100`

    If Sambamba Markdup is invoked using Snakemake, the following bare-bones
    rule should work fine:

    ```
    rule markdup:
      input:
        "data/align/{sample}.bam"
      output:
        "data/markdup/{sample}.markdup.bam"
      log:
        "data/logs/{sample}.log"
      shell:
        "sambamba markdup {input} {output} > {log} 2>&1"
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Sambamba",
            anchor="sambamba",
            href="https://lomereiter.github.io/sambamba/",
            info="Toolkit for interacting with BAM/CRAM files.",
            extra="It is functionally similar to Samtools, but the source code is written in the D Language. "
            "It allows for faster performance while still being easy to use.",
            doi="10.1093/bioinformatics/btv098",
        )

        n = dict()

        # Call submodule functions
        n["markdup"] = parse_sambamba_markdup(self)
        if n["markdup"] > 0:
            log.info(f"Found {n['markdup']} sambamba markdup reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
