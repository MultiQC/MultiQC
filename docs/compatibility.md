# Updating for compatibility

When releasing new versions of MultiQC we aim to maintain compatibility so that your existing
modules and plugins will keep working. However, in some cases we have to make changes that
require code to be modified. This section summarises the changes by MultiQC release.

## v0.10

Modules that referenced ```multiqc.BaseMultiqcModule``` instead need to reference
```multiqc.modules.base_module.BaseMultiqcModule```.

Modules that ```import plots``` need to import the specific plots needed. For example:

```
import plots
plots.beeswarm.plot(...)
```
becomes...
```
from plots import beeswarm
beeswarm.plot(...)
```

