import io
from typing import NewType, Optional, TypedDict, Union

Anchor = NewType("Anchor", str)
ModuleId = NewType("ModuleId", str)
SectionId = NewType("SectionId", str)

ColumnKey = NewType("ColumnKey", str)
SampleName = NewType("SampleName", str)
SampleGroup = NewType("SampleGroup", str)


class FileDict(TypedDict):
    fn: str
    root: str


class LoadedFileDict(FileDict):
    sp_key: str
    s_name: str
    f: Optional[Union[str, io.IOBase]]
