import logging
import os
from typing import Dict, Union

from multiqc import BaseMultiqcModule
from multiqc.modules.qualimap.qualimap import MultiqcModule

__all__ = ["MultiqcModule"]

log = logging.getLogger(__name__)


def parse_numerals(
    preparsed_d: Dict[str, str],
    float_metrics: Dict[str, str],
    int_metrics: Dict[str, str],
    rate_metrics: Dict[str, str],
    fpath: str,
) -> Dict[str, Union[int, float, str]]:
    """
    Take pre-parsed Qualimap report (keys to string values), and properly parse
    numeral values, taking regional formats into account.
    """
    # Determine if decimal separator is dot or comma
    decimalcomma = None
    for k in rate_metrics:
        if k in preparsed_d:
            val = preparsed_d[k]
            if "," in val and "." in val:
                log.error(
                    f"Couldn't determine decimal separator for file {fpath}, as both . and , are "
                    f"found in a rational value: {val}"
                )
                return {}
            if "," in val:
                if decimalcomma is False:
                    log.error(
                        f"Couldn't determine decimal separator for file {fpath}, as differently formatted"
                        f"rational values are found"
                    )
                    return {}
                decimalcomma = True
            if "." in val:
                if decimalcomma is True:
                    log.error(
                        f"Couldn't determine decimal separator for file {fpath}, as differently formatted"
                        f"rational values are found"
                    )
                    return {}
                decimalcomma = False
    if decimalcomma is None:
        # All expected float numbers are integer, so attempt to instead determine the
        # thousands separator from large int values.
        for k in int_metrics:
            if k in preparsed_d:
                val = preparsed_d[k]
                if "," in val and "." in val:
                    log.error(
                        f"Couldn't determine decimal separator for file {fpath}, as both . and , are "
                        f"found in a rational value: {val}"
                    )
                    return {}
                if "," in val:
                    if decimalcomma is True:
                        log.error(
                            f"Couldn't determine decimal separator for file {fpath}, as differently formatted"
                            f"rational values are found"
                        )
                        return {}
                    decimalcomma = False
                if "." in val:
                    if decimalcomma is False:
                        log.error(
                            f"Couldn't determine decimal separator for file {fpath}, as differently formatted"
                            f"rational values are found"
                        )
                        return {}
                    decimalcomma = True
    if decimalcomma is None:
        log.debug(f"Couldn't determine decimal separator for file {fpath}")

    d: Dict = {}
    for k, v in preparsed_d.items():
        v = v.strip("X").strip("%")
        if k in float_metrics or k in rate_metrics or k in int_metrics:
            if decimalcomma is True:
                v = v.replace(".", "").replace(",", ".")
            v = v.replace(",", "")
            if k in int_metrics:
                d[int_metrics[k]] = int(v)
            elif k in float_metrics:
                d[float_metrics[k]] = float(v)
            elif k in rate_metrics:
                d[rate_metrics[k]] = float(v)

    return d


def get_s_name(module: BaseMultiqcModule, f):
    s_name = os.path.basename(os.path.dirname(f["root"]))
    s_name = module.clean_s_name(s_name, f)
    if s_name.endswith(".qc"):
        s_name = s_name[:-3]
    return s_name
