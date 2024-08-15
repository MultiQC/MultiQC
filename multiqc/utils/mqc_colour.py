"""
Helper functions to manipulate colours and colour scales
"""

import functools
import hashlib

# Default logger will be replaced by caller
import logging
import re
from typing import Tuple

import numpy as np
import spectra  # type: ignore

from multiqc import config, report

logger = logging.getLogger(__name__)


@functools.lru_cache(128)  # 34 unique colourmaps found using multiqc-test-data
def cached_spectra_colour_scale(colours: Tuple[str]):
    """Caches spectra color scale calls as these are expensive"""
    return spectra.scale(list(colours))


class mqc_colour_scale(object):
    """Class to hold a colour scheme."""

    # ColorBrewer colours, taken from Chroma.js source code
    # https://github.com/gka/chroma.js
    ###    Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The
    ###    Pennsylvania State University.
    ###
    ###    Licensed under the Apache License, Version 2.0 (the "License");
    ###    you may not use this file except in compliance with the License.
    ###    You may obtain a copy of the License at
    ###    http://www.apache.org/licenses/LICENSE-2.0
    ###
    ###    Unless required by applicable law or agreed to in writing, software distributed
    ###    under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
    ###    CONDITIONS OF ANY KIND, either express or implied. See the License for the
    ###    specific language governing permissions and limitations under the License.
    COLORBREWER_SCALES = {
        # sequential
        "OrRd": ["#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"],
        "PuBu": ["#fff7fb", "#ece7f2", "#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#023858"],
        "BuPu": ["#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b"],
        "Oranges": [
            "#fff5eb",
            "#fee6ce",
            "#fdd0a2",
            "#fdae6b",
            "#fd8d3c",
            "#f16913",
            "#d94801",
            "#a63603",
            "#7f2704",
        ],
        "BuGn": ["#f7fcfd", "#e5f5f9", "#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#006d2c", "#00441b"],
        "YlOrBr": [
            "#ffffe5",
            "#fff7bc",
            "#fee391",
            "#fec44f",
            "#fe9929",
            "#ec7014",
            "#cc4c02",
            "#993404",
            "#662506",
        ],
        "YlGn": ["#ffffe5", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443", "#006837", "#004529"],
        "Reds": ["#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"],
        "RdPu": ["#fff7f3", "#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a"],
        "Greens": [
            "#f7fcf5",
            "#e5f5e0",
            "#c7e9c0",
            "#a1d99b",
            "#74c476",
            "#41ab5d",
            "#238b45",
            "#006d2c",
            "#00441b",
        ],
        "YlGnBu": [
            "#ffffd9",
            "#edf8b1",
            "#c7e9b4",
            "#7fcdbb",
            "#41b6c4",
            "#1d91c0",
            "#225ea8",
            "#253494",
            "#081d58",
        ],
        "Purples": [
            "#fcfbfd",
            "#efedf5",
            "#dadaeb",
            "#bcbddc",
            "#9e9ac8",
            "#807dba",
            "#6a51a3",
            "#54278f",
            "#3f007d",
        ],
        "GnBu": ["#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081"],
        "Greys": [
            "#ffffff",
            "#f0f0f0",
            "#d9d9d9",
            "#bdbdbd",
            "#969696",
            "#737373",
            "#525252",
            "#252525",
            "#000000",
        ],
        "YlOrRd": [
            "#ffffcc",
            "#ffeda0",
            "#fed976",
            "#feb24c",
            "#fd8d3c",
            "#fc4e2a",
            "#e31a1c",
            "#bd0026",
            "#800026",
        ],
        "PuRd": ["#f7f4f9", "#e7e1ef", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256", "#980043", "#67001f"],
        "Blues": [
            "#f7fbff",
            "#deebf7",
            "#c6dbef",
            "#9ecae1",
            "#6baed6",
            "#4292c6",
            "#2171b5",
            "#08519c",
            "#08306b",
        ],
        "PuBuGn": [
            "#fff7fb",
            "#ece2f0",
            "#d0d1e6",
            "#a6bddb",
            "#67a9cf",
            "#3690c0",
            "#02818a",
            "#016c59",
            "#014636",
        ],
        # diverging
        "Spectral": [
            "#9e0142",
            "#d53e4f",
            "#f46d43",
            "#fdae61",
            "#fee08b",
            "#ffffbf",
            "#e6f598",
            "#abdda4",
            "#66c2a5",
            "#3288bd",
            "#5e4fa2",
        ],
        "RdYlGn": [
            "#a50026",
            "#d73027",
            "#f46d43",
            "#fdae61",
            "#fee08b",
            "#ffffbf",
            "#d9ef8b",
            "#a6d96a",
            "#66bd63",
            "#1a9850",
            "#006837",
        ],
        "RdBu": [
            "#67001f",
            "#b2182b",
            "#d6604d",
            "#f4a582",
            "#fddbc7",
            "#f7f7f7",
            "#d1e5f0",
            "#92c5de",
            "#4393c3",
            "#2166ac",
            "#053061",
        ],
        "PiYG": [
            "#8e0152",
            "#c51b7d",
            "#de77ae",
            "#f1b6da",
            "#fde0ef",
            "#f7f7f7",
            "#e6f5d0",
            "#b8e186",
            "#7fbc41",
            "#4d9221",
            "#276419",
        ],
        "PRGn": [
            "#40004b",
            "#762a83",
            "#9970ab",
            "#c2a5cf",
            "#e7d4e8",
            "#f7f7f7",
            "#d9f0d3",
            "#a6dba0",
            "#5aae61",
            "#1b7837",
            "#00441b",
        ],
        "RdYlBu": [
            "#a50026",
            "#d73027",
            "#f46d43",
            "#fdae61",
            "#fee090",
            "#ffffbf",
            "#e0f3f8",
            "#abd9e9",
            "#74add1",
            "#4575b4",
            "#313695",
        ],
        "BrBG": [
            "#543005",
            "#8c510a",
            "#bf812d",
            "#dfc27d",
            "#f6e8c3",
            "#f5f5f5",
            "#c7eae5",
            "#80cdc1",
            "#35978f",
            "#01665e",
            "#003c30",
        ],
        "RdGy": [
            "#67001f",
            "#b2182b",
            "#d6604d",
            "#f4a582",
            "#fddbc7",
            "#ffffff",
            "#e0e0e0",
            "#bababa",
            "#878787",
            "#4d4d4d",
            "#1a1a1a",
        ],
        "PuOr": [
            "#7f3b08",
            "#b35806",
            "#e08214",
            "#fdb863",
            "#fee0b6",
            "#f7f7f7",
            "#d8daeb",
            "#b2abd2",
            "#8073ac",
            "#542788",
            "#2d004b",
        ],
        # qualitative
        "Set2": ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"],
        "Accent": ["#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", "#666666"],
        "Set1": ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"],
        "Set3": [
            "#8dd3c7",
            "#ffffb3",
            "#bebada",
            "#fb8072",
            "#80b1d3",
            "#fdb462",
            "#b3de69",
            "#fccde5",
            "#d9d9d9",
            "#bc80bd",
            "#ccebc5",
            "#ffed6f",
        ],
        "Dark2": ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"],
        "Paired": [
            "#a6cee3",
            "#1f78b4",
            "#b2df8a",
            "#33a02c",
            "#fb9a99",
            "#e31a1c",
            "#fdbf6f",
            "#ff7f00",
            "#cab2d6",
            "#6a3d9a",
            "#ffff99",
            "#b15928",
        ],
        "Pastel2": ["#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9", "#fff2ae", "#f1e2cc", "#cccccc"],
        "Pastel1": [
            "#fbb4ae",
            "#b3cde3",
            "#ccebc5",
            "#decbe4",
            "#fed9a6",
            "#ffffcc",
            "#e5d8bd",
            "#fddaec",
            "#f2f2f2",
        ],
        # Originally from Highcharts
        "plot_defaults": [
            "#7cb5ec",
            "#434348",
            "#90ed7d",
            "#f7a35c",
            "#8085e9",
            "#f15c80",
            "#e4d354",
            "#2b908f",
            "#f45b5b",
            "#91e8e1",
        ],
    }

    def __init__(self, name="GnBu", minval=0, maxval=100, id=None):
        """Initialise class with a colour scale"""

        self.name = name
        self.id = id
        self.colours = self.get_colours(name)

        # Sanity checks
        minval = re.sub(r"[^0-9\.\-e]", "", str(minval))
        maxval = re.sub(r"[^0-9\.\-e]", "", str(maxval))
        if minval == "":
            minval = 0
        if maxval == "":
            maxval = 100
        if float(minval) == float(maxval):
            self.minval = float(minval)
            self.maxval = float(minval) + 1.0
        elif float(minval) > float(maxval):
            self.minval = float(maxval)
            self.maxval = float(minval)
        else:
            self.minval = float(minval)
            self.maxval = float(maxval)

    def get_colour(self, val, colformat="hex", lighten=0.3, source=None):
        """Given a value, return a colour within the colour scale"""

        if val is None:
            return ""

        # Ported from the original JavaScript for continuity
        # Seems to work better than adjusting brightness / saturation / luminosity
        def rgb_converter(x):
            return max(0, min(1, 1 + ((x - 1) * lighten)))

        try:
            if self.name in mqc_colour_scale.qualitative_scales and isinstance(val, float):
                if config.strict:
                    sequential_scales = [
                        s for s in mqc_colour_scale.COLORBREWER_SCALES if s not in mqc_colour_scale.qualitative_scales
                    ]
                    source = f"{source}: " if source else ""
                    errmsg = (
                        f"{source}A categorical scale '{self.name}' is used for float values ({val}). "
                        f"Consider using one of the sequential scales instead: {', '.join(sequential_scales)}"
                    )
                    logger.error(errmsg)
                    report.lint_errors.append(errmsg)
            elif self.name in mqc_colour_scale.qualitative_scales:
                if not isinstance(val, int):
                    # When we have non-numeric values (e.g. Male/Female, Yes/No, chromosome names, etc.), and a qualitative
                    # scale (Set1, Set3, etc.), we don't want to attempt to parse numbers, otherwise we might end up with all
                    # values assigned with the same color. But instead we will get a hash from a string to hope to assign
                    # a unique color for each possible enumeration value.
                    val = deterministic_hash(val)
                thecolour = spectra.html(self.colours[val % len(self.colours)])
                thecolour = spectra.rgb(*[rgb_converter(v) for v in thecolour.rgb])
                return thecolour.hexcode

            # When there is only 1 color in scale, spectra.scale() will crash with DivisionByZero
            elif len(self.colours) == 1:
                thecolour = spectra.html(self.colours[0])
                thecolour = spectra.rgb(*[rgb_converter(v) for v in thecolour.rgb])
                return thecolour.hexcode

            else:
                # Sanity checks
                val_stripped = re.sub(r"[^0-9\.\-e]", "", str(val))
                val_float: float
                if val_stripped == "":
                    val_float = self.minval
                else:
                    try:
                        val_float = float(val_stripped)
                    except ValueError:
                        # No color formatting for non-numeric values
                        return ""
                    val_float = max(val_float, self.minval)
                    val_float = min(val_float, self.maxval)

                domain_nums = list(np.linspace(self.minval, self.maxval, len(self.colours)))
                my_spectra_scale = cached_spectra_colour_scale(tuple(self.colours))
                my_scale = my_spectra_scale.domain(domain_nums)

                # Lighten colours
                thecolour = spectra.rgb(*[rgb_converter(v) for v in my_scale(val_float).rgb])

                return thecolour.hexcode

        except Exception as e:
            # Shouldn't crash all of MultiQC just for colours
            logger.warning(f"{self.id + ': ' if self.id else ''}Error getting colour: {e}")
            return ""

    def get_colours(self, name="GnBu"):
        """Function to get a colour scale by name
        Input: Name of colour scale (suffix with -rev for reversed)
               Defaults to 'GnBu' if scale not found.
        Returns: List of hex colours
        """
        if name.startswith("#"):
            return [name]

        if name in mqc_colour_scale.html_colors:
            return [name]

        # Detect reverse colour scales
        reverse = False
        if str(name).endswith("-rev"):
            reverse = True
            name = name[:-4]

        # Default colour scale
        if name not in mqc_colour_scale.COLORBREWER_SCALES:
            errmsg = f"{self.id + ': ' if self.id else ''}Colour scale {name} not found - defaulting to GnBu"
            if config.strict:
                logger.error(errmsg)
                report.lint_errors.append(errmsg)
            else:
                logger.debug(errmsg)
            name = "GnBu"

        # Return colours
        if reverse:
            return list(reversed(mqc_colour_scale.COLORBREWER_SCALES[name]))
        else:
            return mqc_colour_scale.COLORBREWER_SCALES[name]

    qualitative_scales = ["Set2", "Accent", "Set1", "Set3", "Dark2", "Paired", "Pastel2", "Pastel1", "plot_defaults"]

    html_colors = {
        "black": "#000000",
        "silver": "#C0C0C0",
        "gray": "#808080",
        "grey": "#808080",
        "white": "#FFFFFF",
        "maroon": "#800000",
        "red": "#FF0000",
        "purple": "#800080",
        "fuchsia": "#FF00FF",
        "green": "#008000",
        "lime": "#00FF00",
        "olive": "#808000",
        "yellow": "#FFFF00",
        "navy": "#000080",
        "blue": "#0000FF",
        "teal": "#008080",
        "aqua": "#00FFFF",
        "darkblue": "#00008B",
        "mediumblue": "#0000CD",
        "darkgreen": "#006400",
        "darkcyan": "#008B8B",
        "deepskyblue": "#00BFFF",
        "darkturquoise": "#00CED1",
        "mediumspringgreen": "#00FA9A",
        "springgreen": "#00FF7F",
        "cyan": "#00FFFF",
        "midnightblue": "#191970",
        "dodgerblue": "#1E90FF",
        "lightseagreen": "#20B2AA",
        "forestgreen": "#228B22",
        "seagreen": "#2E8B57",
        "darkslategray": "#2F4F4F",
        "darkslategrey": "#2F4F4F",
        "limegreen": "#32CD32",
        "mediumseagreen": "#3CB371",
        "turquoise": "#40E0D0",
        "royalblue": "#4169E1",
        "steelblue": "#4682B4",
        "darkslateblue": "#483D8B",
        "mediumturquoise": "#48D1CC",
        "indigo": "#4B0082",
        "darkolivegreen": "#556B2F",
        "cadetblue": "#5F9EA0",
        "cornflowerblue": "#6495ED",
        "rebeccapurple": "#663399",
        "mediumaquamarine": "#66CDAA",
        "dimgray": "#696969",
        "dimgrey": "#696969",
        "slateblue": "#6A5ACD",
        "olivedrab": "#6B8E23",
        "slategray": "#708090",
        "slategrey": "#708090",
        "lightslategray": "#778899",
        "lightslategrey": "#778899",
        "mediumslateblue": "#7B68EE",
        "lawngreen": "#7CFC00",
        "chartreuse": "#7FFF00",
        "aquamarine": "#7FFFD4",
        "skyblue": "#87CEEB",
        "lightskyblue": "#87CEFA",
        "blueviolet": "#8A2BE2",
        "darkred": "#8B0000",
        "darkmagenta": "#8B008B",
        "saddlebrown": "#8B4513",
        "darkseagreen": "#8FBC8F",
        "lightgreen": "#90EE90",
        "mediumpurple": "#9370DB",
        "darkviolet": "#9400D3",
        "palegreen": "#98FB98",
        "darkorchid": "#9932CC",
        "yellowgreen": "#9ACD32",
        "sienna": "#A0522D",
        "brown": "#A52A2A",
        "darkgray": "#A9A9A9",
        "darkgrey": "#A9A9A9",
        "lightblue": "#ADD8E6",
        "greenyellow": "#ADFF2F",
        "paleturquoise": "#AFEEEE",
        "lightsteelblue": "#B0C4DE",
        "powderblue": "#B0E0E6",
        "firebrick": "#B22222",
        "darkgoldenrod": "#B8860B",
        "mediumorchid": "#BA55D3",
        "rosybrown": "#BC8F8F",
        "darkkhaki": "#BDB76B",
        "mediumvioletred": "#C71585",
        "indianred": "#CD5C5C",
        "peru": "#CD853F",
        "chocolate": "#D2691E",
        "tan": "#D2B48C",
        "lightgray": "#D3D3D3",
        "lightgrey": "#D3D3D3",
        "thistle": "#D8BFD8",
        "orchid": "#DA70D6",
        "goldenrod": "#DAA520",
        "palevioletred": "#DB7093",
        "crimson": "#DC143C",
        "gainsboro": "#DCDCDC",
        "plum": "#DDA0DD",
        "burlywood": "#DEB887",
        "lightcyan": "#E0FFFF",
        "lavender": "#E6E6FA",
        "darksalmon": "#E9967A",
        "violet": "#EE82EE",
        "palegoldenrod": "#EEE8AA",
        "lightcoral": "#F08080",
        "khaki": "#F0E68C",
        "aliceblue": "#F0F8FF",
        "honeydew": "#F0FFF0",
        "azure": "#F0FFFF",
        "sandybrown": "#F4A460",
        "wheat": "#F5DEB3",
        "beige": "#F5F5DC",
        "whitesmoke": "#F5F5F5",
        "mintcream": "#F5FFFA",
        "ghostwhite": "#F8F8FF",
        "salmon": "#FA8072",
        "antiquewhite": "#FAEBD7",
        "linen": "#FAF0E6",
        "lightgoldenrodyellow": "#FAFAD2",
        "oldlace": "#FDF5E6",
        "magenta": "#FF00FF",
        "deeppink": "#FF1493",
        "orangered": "#FF4500",
        "tomato": "#FF6347",
        "hotpink": "#FF69B4",
        "coral": "#FF7F50",
        "darkorange": "#FF8C00",
        "lightsalmon": "#FFA07A",
        "orange": "#FFA500",
        "lightpink": "#FFB6C1",
        "pink": "#FFC0CB",
        "gold": "#FFD700",
        "peachpuff": "#FFDAB9",
        "navajowhite": "#FFDEAD",
        "moccasin": "#FFE4B5",
        "bisque": "#FFE4C4",
        "mistyrose": "#FFE4E1",
        "blanchedalmond": "#FFEBCD",
        "papayawhip": "#FFEFD5",
        "lavenderblush": "#FFF0F5",
        "seashell": "#FFF5EE",
        "cornsilk": "#FFF8DC",
        "lemonchiffon": "#FFFACD",
        "floralwhite": "#FFFAF0",
        "snow": "#FFFAFA",
        "lightyellow": "#FFFFE0",
        "ivory": "#FFFFF0",
    }


def deterministic_hash(x):
    """
    Deterministic hash function for strings. This is useful for assigning a unique color
    to each possible value of a categorical variable.

    Unlike the built-in hash function, this function always returns the same value for
    the same input string (see https://docs.python.org/3/using/cmdline.html#cmdoption-R)
    """
    return int(hashlib.sha1(x.encode("utf-8")).hexdigest(), 16) % (10**8)
