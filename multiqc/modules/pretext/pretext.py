#!/usr/bin/env python

""" MultiQC module to parse output from BUSCO """

from __future__ import print_function
from collections import OrderedDict
import logging
import re, base64
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """BUSCO module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="PRETEXT",
            anchor="pretext",
            href="https://github.com/wtsi-hpag/PretextView",
            info="is a tool for visualizing HiC contact maps",
        )
        self.pretext_data = dict()
        for f in self.find_log_files("pretext"):
            self.pretext_data[f["s_name"].replace("_pretext","")]=(base64.b64encode(f["f"].read()).decode("utf-8"))
        image_format = "png"
        log.info("Found {} pretext image".format(len(self.pretext_data)))
        for image_name,image_string in self.pretext_data.items():
            img_html = (
            '<div class="mqc-custom-content-image"><img src="data:image/{};base64,{}" /></div>'.format(image_format, image_string ))
            self.add_section (
                name=image_name,
                anchor="pretext",
                content = img_html
            )
