# Updating for compatibility

When releasing new versions of MultiQC we aim to maintain compatibility so that your existing
modules and plugins will keep working. However, in some cases we have to make changes that
require code to be modified. This section summarises the changes by MultiQC release.

## v1.0 Updates
MultiQC v1.0 brings a few changes that could break external plugins. These are:

* Module import refactoring to allow a new testing environment
  * This should allow better, more modular, unit testing. This should equate to more
    reliable and maintainable code
  * All modules need to change some of their import statements. This includes plugin
    modules outside of the core MultiQC package.
  * Many thanks to [@tbooth](https://github.com/tbooth) at
    [@EdinburghGenomics](https://github.com/EdinburghGenomics) for his patient work with this.

There are two things that you probably need to change in your plugin modules to
make them work with the updated version of MultiQC, both to do with imports.
Instead of this style of importing modules:
```python
from multiqc import config, BaseMultiqcModule, plots
```

You now need this:
```python
from multiqc import config
from multiqc.plots import bargraph   # Load specific plot types here
from multiqc.modules.base_module import BaseMultiqcModule
```

Modules that directly reference `multiqc.BaseMultiqcModule` instead need to reference
`multiqc.modules.base_module.BaseMultiqcModule`.

Secondly, modules that use `import plots` now need to import the specific plots needed.
You will also need to update any plotting functions, removing the `plot.` prefix.

For example, change this:
```python
import plots
return plots.bargraph.plot(data, keys, pconfig)
```
to this:
```python
from plots import bargraph
return bargraph.plot(data, keys, pconfig)
```

These changes have been made to simplify the module imports within MultiQC,
allowing specific parts of the codebase to be imported into a Python script
on their own. This enables small, atomic, clean unit testing.

If you have any questions, please open an issue.
