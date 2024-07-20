from typing import Dict, Union, Tuple

import xml.etree.ElementTree
import re


def parse_qcml_by(qcml_contents, tag) -> Tuple[Dict[str, Union[float, str]], Dict[str, Tuple[str, str]]]:
    """Parse a qcML file and return key-value pairs from the quality parameter entries."""
    root = xml.etree.ElementTree.fromstring(qcml_contents)
    values: Dict[str, Union[float, str]] = dict()
    params: Dict[str, Tuple[str, str]] = dict()

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
