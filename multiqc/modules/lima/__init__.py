from __future__ import absolute_import

import click

from .lima import MultiqcModule

lima_barcodes = click.option(
        '--lima-barcodes',
        type=click.Path(exists=True),
        required=False,
        help=
        """
        Tab-separated file containing the mapping from barcodes to sample names
        in the format:\n
        samplename first_barcode second_barcode

        This is used to determine the sample names in the lima output. If this
        is not specified, the barcode names will be used as sample names for
        the lima results.
        """
)
