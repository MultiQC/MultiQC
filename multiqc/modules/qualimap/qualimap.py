import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module supports the Qualimap commands `BamQC` and `RNASeq`.
    Note that Qualimap must be run with the `-outdir` option as well as
    `-outformat HTML` (which is on by default). MultiQC uses files
    found within the `raw_data_qualimapReport` folder (as well as `genome_results.txt`).

    Qualimap adds lots of columns to the General Statistics table. To avoid making the table
    too wide and bloated, some of these are hidden by default (`Error Rate`, `M Aligned`, `M Total reads`).
    You can override these defaults in your MultiQC config file - for example, to show
    `Error Rate` by default and hide `Ins. size` by default, add the following:

    ```yaml
    table_columns_visible:
      QualiMap:
        general_error_rate: True
        median_insert_size: False
    ```

    See the [relevant section of the documentation](http://multiqc.info/docs/#hiding-columns) for more detail.

    In addition to this, it's possible to customise which coverage thresholds calculated
    by the Qualimap BamQC module _(default: 1, 5, 10, 30, 50)_ and which of these are hidden in the
    General Statistics tablewhen the report loads _(default: all hidden except 30X)_.

    To do this, add something like the following to your MultiQC config file:

    ```yaml
    qualimap_config:
      general_stats_coverage:
        - 10
        - 20
        - 40
        - 200
        - 30000
      general_stats_coverage_hidden:
        - 10
        - 20
        - 200
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="QualiMap",
            anchor="qualimap",
            href="http://qualimap.bioinfo.cipf.es/",
            info="Quality control of alignment data and its derivatives like feature counts.",
            doi=["10.1093/bioinformatics/btv566", "10.1093/bioinformatics/bts503"],
        )

        # Initialise the submodules
        from . import QM_BamQC, QM_RNASeq

        # Set up class objects to hold parsed data()
        n = dict()

        # Call submodule functions
        n["BamQC"] = QM_BamQC.parse_reports(self)
        if n["BamQC"] > 0:
            log.info(f"Found {n['BamQC']} BamQC reports")

        n["RNASeq"] = QM_RNASeq.parse_reports(self)
        if n["RNASeq"] > 0:
            log.info(f"Found {n['RNASeq']} RNASeq reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
