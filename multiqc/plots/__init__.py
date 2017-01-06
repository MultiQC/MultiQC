""" HighCharts plots need a globally unique ID embedded in the HTML.
    - id(pconfig) seems the obvious solution but this might be recycled if
    a new pconfig object is generated at the same memory address.
    - random letters work ok but it seems messy to put that code in every chart
    module.
    - UUID seems overkill
    - asking the report for a UID would work but I want to sever the dependency
    between plots and report classes
"""

_uid = 0

def get_uid():
    global _uid
    _uid += 1
    return '{:04d}'.format(_uid)
