import logging

from multiqc.base_module import BaseMultiqcModule

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="CheckM2",
            anchor="checkm2",
            href="https://github.com/chklovski/CheckM2",
            info="A rapid, scalable and accurate tool for assessing microbial genome quality using machine learning.",
            doi=["10.1038/s41592-023-01940-w"],
        )
        print('Hello world')