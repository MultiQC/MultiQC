# Updating for compatibility

When releasing new versions of MultiQC we aim to maintain compatibility so that your existing
modules and plugins will keep working. However, in some cases we have to make changes that
require code to be modified. This section summarises the changes by MultiQC release.

## v1.0 Updates
MultiQC v1.0 brings a few changes in the way that MultiQC modules and plugins are written. Most are backwards-compatible, but there are a couple that could break external plugins.

#### Module imports

New MultiQC module imports have been refactored to make them less inter-dependent and fragile. This has a bunch of advantages, notably allowing better, more modular, unit testing (and hopefully more reliable and maintainable code).

All MultiQC modules and plugins will need to change some of their import statements.

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

> Many thanks to [@tbooth](https://github.com/tbooth) at [@EdinburghGenomics](https://github.com/EdinburghGenomics) for his patient work with this.

#### Searching for files
The core `find_log_files` function has been rewritten and now works a little differently. Instead of searching all analysis files each time it's called (by every module), all files are searched once at the start of the MultiQC execution. This makes MultiQC run much faster.

To use the new syntax, add your search pattern to `config.sp` using the new `before_config` plugin hook:

`setup.py`:
```python
# [..]
  'multiqc.hooks.v1': [
    'before_config = myplugin.mymodule:load_config'
  ]
```

`mymodule.py`:
```python
from multiqc.utils import config
def load_config():
    my_search_patterns = {
        'my_plugin/my_mod': {'fn': '*_somefile.txt'},
        'my_plugin/my_other_mod': {'fn': '*other_file.txt'},
    }
    config.update_dict(config.sp, my_search_patterns)
```

This will add in your search patterns to the default MultiQC config, before user config files are loaded (allowing people to overwrite your defaults as with other modules).

Now, you can find your files much as before, using the string specified above:

```python
for f in self.find_log_files('my_plugin/my_mod'):
  # do something
```

The old syntax (supplying a `dict` instead of a string to the function without any previous config setup) will still work, but you will get a depreciation notice. This functionality may be removed in the future.

#### Adding report sections
Until now, report sections were added by creating a list called `self.sections` and adding to it. If you only had a single section, the routine was to instead append to the `self.intro` string.

These methods have been depreciated in favour of a new function called `self.add_section()`. For example, instead of the previous:

```python
self.sections = list()
self.sections.append({
  'name': 'My Section',
  'anchor': 'my-html-id',
  'content': '<p>Description of what this plot shows.</p>' +
             linegraph.plot(data, pconfig)
})
```

the syntax is now:

```python
self.add_section(
  name = 'My Section',
  anchor = 'my-html-id',
  description = 'Description of what this plot shows.',
  helptext = 'More extensive help text can about how to interpret this.'
  plot = linegraph.plot(data, pconfig)
)
```

Note that content should now be split up into three new keys: `description`, `helptext` and `plot`. This will allow consistent formatting and future developments with improved module help text. Text is wrapped in `<p>` tags by the function, so these are no longer needed. Raw content can still be provided in a `content` string as before if required.

All fields are optional. If `name` is omitted then the end result will be the same as previously done with `self.intro += content`.

#### Updated number formatting
A couple of minor updates to how numbers are handled in tables may affect your configs. Firstly, format strings looking like `{:.1f}` should now be `{:,.1f}` (note the extra comma). This enables customisable number formatting with separated thousand groups.

Secondly, any table columns reporting a read count should use new config options to allow user-configurable multipliers. For example, instead of this:

```python
headers['read_counts'] = {
  'title': 'M Reads',
  'description': 'Read counts (millions)',
  'modify': lambda x: x / 1000000,
  'format': '{:.,2f} M',
  'shared_key': 'read_count'
}
```

you should now use this:

```python
headers['read_counts'] = {
  'title': '{} Reads'.format(config.read_count_prefix),
  'description': 'Total raw sequences ({})'.format(config.read_count_desc),
  'modify': lambda x: x * config.read_count_multiplier,
  'format': '{:,.2f} ' + config.read_count_prefix,
  'shared_key': 'read_count'
}
```

Not as pretty, but allows users to view low depth coverage.

