import sys

# NewType introduced in 3.10.
# Previous version had typing_extensions.NewType, but it's not pickable.
if sys.version_info[1] < 10:
    from typing import NewType
else:
    from typing_extensions import NewType

AnchorT = NewType("AnchorT", str)
ModuleIdT = NewType("ModuleIdT", str)
SectionIdT = NewType("SectionIdT", str)

ColumnKeyT = NewType("ColumnKeyT", str)
SampleNameT = NewType("SampleNameT", str)
SampleGroupT = NewType("SampleGroupT", str)
