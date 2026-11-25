import re

EVENT_TYPES: dict[str, str] = {
    "SE": "Skipped exon",
    "A5": "Alternative 5' splice site",
    "A3": "Alternative 3' splice site",
    "MX": "Mutually exclusive exons",
    "RI": "Retained intron",
    "AF": "Alternative first exon",
    "AL": "Alternative last exon",
}
ISOFORM = "Isoform"

EVENT_TYPE_COLOURS: dict[str, str] = {
    "Skipped exon": "#437bb1",
    "Alternative 5' splice site": "#e8853a",
    "Alternative 3' splice site": "#3a9a5b",
    "Mutually exclusive exons": "#8e5ea2",
    "Retained intron": "#c0392b",
    "Alternative first exon": "#2a9d8f",
    "Alternative last exon": "#9c6b30",
}

# Local event ids are `<gene>;<type>:<coordinates>`, isoform ids are `<gene>;<transcript>`
LOCAL_EVENT_ID_RE = re.compile(r"^[^;\t]+;(" + "|".join(EVENT_TYPES) + r"):\S+$")


def event_type(event_id: str) -> str:
    """Return the event type label for a SUPPA event id, or `ISOFORM` for a transcript id."""
    match = LOCAL_EVENT_ID_RE.match(event_id)
    if match:
        return EVENT_TYPES[match.group(1)]
    return ISOFORM
