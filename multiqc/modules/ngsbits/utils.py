import re
import xml.etree.ElementTree


def parse_qcml_by(qcml_contents: str, tag: str) -> tuple[dict[str, float | str], dict[str, tuple[str, str]]]:
    """Parse a qcML file and return key-value pairs from the quality parameter entries."""
    root = xml.etree.ElementTree.fromstring(qcml_contents)
    values: dict[str, float | str] = {}
    params: dict[str, tuple[str, str]] = {}

    for qp in root.findall(f".//{{http://www.prime-xs.eu/ms/qcml}}{tag}"):
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
