import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .diffsplice import parse_suppa_diffsplice
from .psi import parse_suppa_psi

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    SUPPA2 computes PSI (percent spliced in) values for local alternative splicing events
    and transcript isoforms from transcript quantification, and tests them for
    differential splicing between conditions.

    Supported commands:

    - `psiPerEvent` and `psiPerIsoform`: `.psi` files
    - `diffSplice`: `.dpsi` files

    #### psiPerEvent and psiPerIsoform

    The module reads the PSI matrix written by `psiPerEvent` (`<output>.psi`) and
    `psiPerIsoform` (`<output>_isoform.psi`). Sample names come from the header, so one
    file can hold many samples. Files of local events and of isoforms are told apart by
    their ids (`<gene>;SE:...` against `<gene>;<transcript>`) and shown as separate tabs.
    For each sample it shows how many events have a PSI value, and the distribution of
    PSI values.

    The per-condition `.psi` files that are passed to `diffSplice` have the same format.
    Pass MultiQC either the full matrix or the per-condition files, not both, or the
    samples will be reported twice.

    #### diffSplice

    The module reads the `<output>.dpsi` file written by `diffSplice` and counts, for each
    comparison, the events tested and the events with a significant change. By default
    an event is significant with p < 0.05 and |dPSI| > 0.05, the defaults of SUPPA
    `clusterEvents`.

    Comparisons are named after the file, without the `.dpsi` extension. Files with more
    than one comparison (three or more conditions) add the comparison label from the
    header, `<condition 1>-<condition 2>`.

    Genome-wide `.psi` files with many samples can exceed the default MultiQC file size
    limit of 50 MB. Raise `log_filesize_limit` in the MultiQC config if they are skipped.

    #### Configuration

    Set these options under `suppa` in a MultiQC config file:

    - `pvalue_cutoff` (default `0.05`): a tested event is significant when its `diffSplice`
      p-value is below this value...
    - `dpsi_cutoff` (default `0.05`): ...and its absolute dPSI is above this value.

    ```yaml
    suppa:
      pvalue_cutoff: 0.05
      dpsi_cutoff: 0.05
    ```
    """

    def __init__(self):
        super().__init__(
            name="SUPPA2",
            anchor="suppa",
            href="https://github.com/comprna/SUPPA",
            info="Fast and accurate quantification and differential analysis of alternative splicing from transcript abundances.",
            doi=["10.1186/s13059-018-1417-1", "10.1261/rna.051557.115"],
            license="MIT License",
            license_url="https://github.com/comprna/SUPPA/blob/master/LICENSE.txt",
        )

        n = {
            "psi": parse_suppa_psi(self),
            "diffsplice": parse_suppa_diffsplice(self),
        }
        for subtool, count in n.items():
            if count > 0:
                log.info(f"Found {count} {subtool} reports")

        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
