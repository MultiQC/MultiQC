import sys


if sys.version_info[1] >= 10:
    # NewType introduced in 3.10.
    # Previous version had typing_extensions.NewType, but it's not pickable.
    from typing import NewType

    AnchorT = NewType("AnchorT", str)
    ModuleIdT = NewType("ModuleIdT", str)
    SectionIdT = NewType("SectionIdT", str)

    ColumnKeyT = NewType("ColumnKeyT", str)
    SampleNameT = NewType("SampleNameT", str)
    SampleGroupT = NewType("SampleGroupT", str)

else:

    class AnchorT(str):
        pass

    class ModuleIdT(str):
        pass

    class SectionIdT(str):
        pass

    class ColumnKeyT(str):
        pass

    class SampleNameT(str):
        pass

    class SampleGroupT(str):
        pass
