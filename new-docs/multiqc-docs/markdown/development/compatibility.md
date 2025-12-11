---
title: Breaking changes
description: Notes about major MultiQC updates
---

# Updating after major changes

When releasing new versions of MultiQC we aim to maintain compatibility so that your existing modules and plugins will keep working.
However, in some cases we have to make changes that require code to be modified.
This section describes any major breaking changes in MultiQC releases.

## GitHub Repo and Docker images moved

On December 18th 2023, just after the v1.19 release of MultiQC, the GitHub repo and Docker images were moved.

- GitHub source code:
  - Repo name `ewels/MultiQC` became `MultiQC/MultiQC`
  - Default branch `master` was renamed to `main`
- Docker images renamed:
  - DockerHub: `ewels/multiqc` became `multiqc/multiqc`
  - GitHub Packages: `ghcr.io/ewels/multiqc` became `ghcr.io/multiqc/multiqc`

### GitHub repository name change

The [`ewels/MultiQC` repository](https://github.com/ewels/MultiQC) was moved to [`MultiQC/MultiQC`](https://github.com/MultiQC/MultiQC).
All issues, pull-requests and other GitHub metadata moved with it.
In most cases, GitHub should automatically redirect to the new location and you should not notice any difference.
However, if you have a local clone / fork of the repository it's still a good idea to update the remote name.

Assuming that you have forked the repo and have a local clone, configured with a remote called `upstream` with which you pull in new changes (with `git pull upstream`), you can rename the remote with the following:

```bash
git remote set-url upstream git@github.com:MultiQC/MultiQC.git
```

### Default branch name changed

At the same time as moving the repo, we also changed the default branch name from `master` to `main`. This is inline with changing industry standards.
See [`github/renaming`](https://github.com/github/renaming/) and [this software freedom conservancy blog post](https://sfconservancy.org/news/2020/jun/23/gitbranchname/) for more details.

If you maintain your own fork of MultiQC, it will be unaffected by this change. The branch name is only switched on the main MultiQC repository. All pull-requests have been automatically updated to point to the renamed branch.

If you wish to change the default branch name on your fork, you can do (probably a good idea so as not to get confusing).
First, rename the branch on GitHub.com (_Settings_ -> _Default branch_, see also the [GitHub docs](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-branches-in-your-repository/renaming-a-branch))

Once renamed on GitHub.com, you'll need to rename the branch in your local clone and point to the renamed remote branch:

```bash
git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```

### New Docker images

To coincide with the GitHub repository renaming, the official Docker images have also been renamed.

**The previous `ewels/multiqc` images have _not_ been removed**, so as not to hinder reproducibility.
However, they will no longer be updated and get no new pushes to `:dev`, `:latest` or releases after `:v1.19`.

New docker images have been created [at DockerHub (`multiqc/multiqc`)](https://hub.docker.com/r/multiqc/multiqc/)
and at [GitHub Packages (`ghcr.io/multiqc/multiqc`)](https://github.com/MultiQC/MultiQC/pkgs/container/multiqc).
These do not have old versions, but will be updated from now on with releases v1.20 and onwards.

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

Report sections can be added by the method called `self.add_section()`. For example:

```python
self.add_section(
  name='My Section',
  anchor='my-html-id',
  description='Description of what this plot shows.',
  helptext='More extensive help text can about how to interpret this.',
  plot=linegraph.plot(data, pconfig),
)
```

Text passed as `description` and `helptext` is wrapped in `<p>` tags, and additional raw content can be provided with a `content` or `content_before_plot` string if required.

#### Updated number formatting

A couple of minor updates to how numbers are handled in tables may affect your configs. Firstly, format strings looking like `{:.1f}` should now be `{:,.1f}` (note the extra comma). This enables customisable number formatting with separated thousand groups. For example, `{:,.2f}` will format `1234567.89` as `1,234,567.89`. For decimal numbers, use `{:,d}`.

Secondly, any table columns reporting a read and base counts should use new config options to allow user-configurable multipliers. For example, instead of this:

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
  'title': f'{config.read_count_prefix} Reads',
  'description': f'Total raw sequences ({config.read_count_desc})',
  'modify': lambda x: x * config.read_count_multiplier,
  'format': '{:,.2f} ' + config.read_count_prefix,
  'shared_key': 'read_count'
}
```

Not as pretty, but allows users to view low depth coverage.

Similarly, for base counts:

```python
headers['base_counts'] = {
  'title': f'{config.base_count_prefix} Bases',
  'description': f'Total raw bases ({config.base_count_desc})',
  'modify': lambda x: x * config.base_count_multiplier,
  'format': '{:,.2f} ' + config.base_count_prefix,
  'shared_key': 'base_count'
}
```
