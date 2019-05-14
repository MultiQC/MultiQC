from __future__ import absolute_import
from .seqyclean import MultiqcModule

from multiqc.modules.base_module import BaseMultiqcModule

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='My Module', anchor='mymod',
        href="http://www.awesome_bioinfo.com/my_module",
        info="is an example analysis module used for writing documentation.")
