import sys


if sys.version_info[1] >= 10:
    # NewType introduced in 3.10.
    # Previous version had typing_extensions.NewType, but it's not pickable.
    from typing import NewType

    AnchorT = NewType("AnchorT", str)
    ModuleIdT = NewType("ModuleIdT", str)
    SectionIdT = NewType("SectionIdT", str)

else:
    AnchorT = str
    ModuleIdT = str
    SectionIdT = str
