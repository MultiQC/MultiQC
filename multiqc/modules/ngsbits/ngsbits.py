#!/usr/bin/env python

""" MultiQC module to parse output from ngs-bits """

import logging
import xml.etree.ElementTree
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule


# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ngs-bits",
            anchor="ngsbits",
            href="https://github.com/imgag/ngs-bits",
            info="is a collection of short-read sequencing tools.",
            doi="10.1093/bioinformatics/btx032",
        )

        # Get the list of submodules (can be customised)
        ngsbits_sections = getattr(config, "ngsbits_sections", [])
        if len(ngsbits_sections) == 0:
            ngsbits_sections = [
                "readqc",
                "mappingqc",
            ]

        # Call submodule functions
        n = dict()
        for sm in ngsbits_sections:
            try:
                # Import the submodule and call parse_reports()
                #   Function returns number of parsed logs
                module = __import__(f"multiqc.modules.ngsbits.{sm}", fromlist=[""])
                n[sm] = getattr(module, "parse_reports")(self)
                if n[sm] > 0:
                    log.info(f"Found {n[sm]} {sm} reports")
            except (ImportError, AttributeError):
                log.warning(f"Could not find ngs-bits section '{sm}'")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

    @staticmethod
    def parse_qcml_by(qcml_contents, tag):
        """Parse a qcML file and return key-value pairs from the quality parameter entries."""
        root = xml.etree.ElementTree.fromstring(qcml_contents)
        values = dict()
        params = dict()

        for qp in root.findall(".//{http://www.prime-xs.eu/ms/qcml}%s" % tag):
            # skip n/a values
            if qp.attrib["value"].startswith("n/a"):
                continue

            # replace 'percentage' with '%'
            qp_name = re.sub(r" percentage$", " %", qp.attrib["name"])

            try:
                values[qp_name] = float(qp.attrib["value"])
            except ValueError:
                values[qp_name] = qp.attrib["value"]

            # add description and accession number of the parameter to the header
            params[qp_name] = (qp.attrib["description"], qp.attrib["accession"])
        return values, params
