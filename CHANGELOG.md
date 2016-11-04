# MultiQC Version History

#### v0.9dev
* **Cutadapt**
  * Fixed bug in General Stats table number for old versions of cutadapt (pre v1.7)
  * Added support for _really_ old cutadapt logs (eg. v.1.2)
* **RSeQC**
  * Made the gene body coverage plot show the percentage view by default
* **Samtools**
  * New module to show duplicate stats from `rmdup` logs
  * Fixed a couple of niggles in the idxstats plot
* **SnpEff**
  * Fixed swapped axis labels in the Variant Quality plot
* **STAR**
  * Fixed crash when there are 0 unmapped reads.

Core Updates:
* Fixed bar plot bug where missing categories could shift data between samples
* Missing plot IDs added for easier plot exporting
* Stopped giving warnings about skipping directories (now a debug message)

#### [v0.8](https://github.com/ewels/MultiQC/releases/tag/v0.8) - 2016-09-26
Module updates:
* [**GATK**](https://software.broadinstitute.org/gatk/) - new module!
  * Added support for VariantEval reports, only parsing a little of the information
    in there so far, but it's a start.
  * Module originally written by @robinandeer at the [OBF Codefest](https://www.open-bio.org/wiki/Codefest_2016),
    finished off by @ewels
* [**Bcftools**](https://samtools.github.io/bcftools/) - new module!
* [**QUAST**](http://quast.bioinf.spbau.ru/) - new module!
  * QUAST is a tool for assessing de novo assemblies against reference genomes.
* **Bismark** now supports reports from `bam2nuc`, giving Cytosine coverage in General Stats.
* **Bowtie1**
  * Updated to try to find bowtie command before log, handle multiple logs in one file. Same as bowtie2.
* **FastQC**
  * Sample pass/warn/fail lists now display properly even with large numbers of samples
  * Sequence content heatmap display is better with many samples
* **Kallisto**
  * Now supports logs from SE data.
* **Picard**
  * `BaseDistributionByCycle` - new submodule! Written by @mlusignan
  * `RnaSeqMetrics` - new submodule! This one by @ewels ;)
  * `AlignmentSummaryMetrics` - another new submodule!
  * Fixed truncated files crash bug for Python 3 _(#306)_
* **Qualimap RNASeqQC**
  * Fixed parsing bug affecting counts in _Genomic Origin_ plot.
  * Module now works with European style thousand separators (`1.234,56` instead of `1,234.56`)
* **RSeQC**
  * `infer_experiment` - new submodule! Written by @Hammarn
* **Samtools**
  * `stats` submodule now has separate bar graph showing alignment scores
  * `flagstat` - new submodule! Written by @HLWiencko
  * `idxstats` - new submodule! This one by @ewels again

Core updates:
* New `--export`/`-p` option to generate static images plot in `multiqc_plots` (`.png`, `.svg` and `.pdf`)
  * Configurable with `export_plots`, `plots_dir_name` and `export_plot_formats` config options
  * `--flat` option no longer saves plots in `multiqc_data/multiqc_plots`
* New `--comment`/`-b` flag to add a comment to the top of reports.
* New `--dirs-depth`/`-dd` flag to specify how many directories to prepend with `--dirs`/`-d`
  * Specifying a postive number will take that many directories from the end of the path
  * A negative number will take directories from the start of the path.
* Directory paths now appended before cleaning, so `fn_clean_exts` will now affect these names.
* New `custom_logo` attributes to add your own logo to reports.
* New `report_header_info` config option to add arbitrary information to the top of reports.
* New `--pdf` option to create a PDF report
  * Depends on [Pandoc](http://pandoc.org) being installed and is in a beta-stage currently.
  * Note that specifying this will make MultiQC use the `simple` template, giving a HTML report with
    much reduced functionality.
* New `fn_clean_sample_names` config option to turn off sample name cleaning
  * This will print the full filename for samples. Less pretty reports and rows
    on the General Statistics table won't line up, but can prevent overwriting.
* Table header defaults can now be set easily
* General Statistics table now hidden if empty.
* Some new defaults in the sample name cleaning
* Updated the `simple` template.
  * Now has no toolbox or nav, no JavaScript and is better suited for printing / PDFs.
  * New `config.simple_output` config flag so code knows when we're trying to avoid JS.
* Fixed some bugs with config settings (eg. template) being overwritten.
* NFS log file deletion bug fixed by @brainstorm (#265)
* Fixed bug in `--ignore` behaviour with directory names.
* Fixed nasty bug in beeswarm dot plots where sample names were mixed up (#278)
* Beeswarm header text is now more informative (sample count with more info on a tooltip)
* Beeswarm plots now work when reports have > 1000 samples

Many thanks to those at the [OpenBio Codefest 2016](https://www.open-bio.org/wiki/Codefest_2016)
who worked on MultiQC projects.

#### [v0.7](https://github.com/ewels/MultiQC/releases/tag/v0.7) - 2016-07-04
Module updates:
* [**Kallisto**](https://pachterlab.github.io/kallisto/) - new module!
* **Picard**
  * Code refactored to make maintenance and additions easier.
  * Big update to `HsMetrics` parsing - more results shown in report, new plots (by @lpantano)
  * Updated `InsertSizeMetrics` to understand logs generated by `CollectMultipleMetrics` (#215)
  * Newlines in picard output. Fixed by @dakl
* **Samtools**
  * Code refactored
  * Rewrote the `samtools stats` code to display more stats in report with a beeswarm plot.
* **Qualimap**
  * Rewritten to use latest methods and fix bugs.
  * Added _Percentage Aligned_ column to general stats for `BamQC` module.
  * Extra table thresholds added by @avilella (hidden by default)
* **General Statistics**
  * Some tweaks to the display defaults (FastQC, Bismark, Qualimap, SnpEff)
  * Now possible to skip the General Statistics section of the report with `--exclude general_stats`
* **Cutadapt** module updated to recognise logs from old versions of cutadapt (<= v1.6)
* **Trimmomatic**
  * Now handles `,` decimal places in percentage values.
  * Can cope with line breaks in log files (see issue #212)
* **FastQC** refactored
  * Now skips zip files if the sample name has already been found. Speeds up MultiQC execution.
  * Code cleaned up. Parsing and data-structures standardised.
  * New popovers on Pass / Warn / Fail status bars showing sample names. Fast highlighting and hiding.
  * New column in General Stats (hidden by default) showing percentage of FastQC modules that failed.
* **SnpEff**
  * Search pattern now more generic, should match reports from others.
  * _Counts by Effect_ plot removed (had hundreds of categories, was fairly unusable).
  * `KeyError` bug fixed.
* **Samblaster** now gets sample name from `ID` instead of `SM` (@dakl)
* **Bowtie 2**
  * Now parses overall alignment rate as intended.
  * Now depends on even less log contents to work with more inputs.
* **MethylQA** now handles variable spacing in logs
* **featureCounts** now splits columns on tabs instead of whitespace, can handle filenames with spaces

Core updates:
* **Galaxy**: MultiQC now available in Galax! Work by @devengineson / @yvanlebras / @cmonjeau
  * See it in the [Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/view/engineson/multiqc/)
* **Heatmap**: New plot type!
* **Scatter Plot**: New plot type!
* **Download raw data** behind plots in reports! Available in the Export toolbox.
  * Choose from tab-separated, comma-separated and the complete JSON.
* **Table columns can be hidden** on page load (shown through _Configure Columns_)
  * Defaults are configurable using the `table_columns_visible` config option.
* **Beeswarm plot**: Added missing rename / highlight / hiding functionality.
* New `-l` / `--file-list` option: specify a file containing a **list of files** to search.
* **Updated HighCharts** to v4.2.5. Added option to export to JPEG.
* Can now **cancel execution** with a single `ctrl+c` rather than having to button mash
* More granular control of **skipping files** during scan (filename, dirname, path matching)
  * Fixed `--exclude` so that it works with directories as well as files
* **New _Clear_ button** in toolbox to bulk remove highlighting / renaming / hiding filters.
* Improved documentation about behaviour for large sample numbers.
* Handle YAML parsing errors for the config file more gracefully
* Removed empty columns from tables again
* Fixed bug in changing module search patterns, reported by @lweasel
* Added timeout parameter to version check to prevent hang on systems with long defaults
* Fixed table display bug in Firefox
* Fixed bug related to order in which config files are loaded
* Fixed bug that broke the _"Show only"_ toolbox feature with multiple names.
* Numerous other small bugs.


#### [v0.6](https://github.com/ewels/MultiQC/releases/tag/v0.6) - 2016-04-29
Module updates:
* New [Salmon](http://combine-lab.github.io/salmon/) module.
* New [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) module.
* New [Bamtools stats](https://github.com/pezmaster31/bamtools) module.
* New beeswarm plot type. General Stats table replaced with this when many samples in report.
* New RSeQC module: Actually a suite of 8 new modules supporting various outputs from RSeQC
* Rewrote bowtie2 module: Now better at parsing logs and tries to scrape input from wrapper logs.
* Made cutadapt show counts by default instead of obs/exp
* Added percentage view to Picard insert size plot

Core updates:
* Dynamic plots now update their labels properly when changing datasets and to percentages
* Config files now loaded from working directory if present
* Started new docs describing how each module works
* Refactored featureCounts module. Now handles summaries describing multiple samples.
* Stopped using so many hidden files. `.multiqc.log` now called `multiqc.log`
* New `-c`/`--config` command line option to specify a MultiQC configuration file
* Can now load run-specific config files called `multiqc_config.yaml` in working directory
* Large code refactoring - moved plotting code out of `BaseModule` and into new `multiqc.plots` submodules
* Generalised code used to generate the General Stats table so that it can be used by modules
* Removed interactive report tour, replaced with a link to a youtube tutorial
* Made it possible to permanently hide the blue welcome message for all future reports
* New option to smooth data for line plots. Avoids mega-huge plots. Applied to SnpEff, RSeQC, Picard.

Bugfixes:
* Qualimap handles infinity symbol (thanks @chapmanb )
* Made SnpEff less fussy about required fields for making plots
* UTF-8 file paths handled properly in Py2.7+
* Extending two config variables wasn't working. Now fixed.
* Dragging the height bar of plots now works again.
* Plots now properly change y axis limits and labels when changing datasets
* Flat plots now have correct path in `default_dev` template

#### [v0.5](https://github.com/ewels/MultiQC/releases/tag/v0.5) - 2016-03-29
Module updates:
* New [Skewer](https://github.com/relipmoc/skewer) module, written by @dakl
* New [Samblaster](https://github.com/GregoryFaust/samblaster) module, written by @dakl
* New [Samtools stats](http://www.htslib.org/) module, written by @lpantano
* New [HiCUP](http://www.bioinformatics.babraham.ac.uk/projects/hicup/) module
* New [SnpEff](http://snpeff.sourceforge.net/) module
* New [methylQA](http://methylqa.sourceforge.net/) module

Core updates:
* New "Flat" image plots, rendered at run time with MatPlotLib
  * By default, will use image plots if > 50 samples (set in config as `plots_flat_numseries`)
  * Means that _very_ large numbers of samples can be viewed in reports. _eg._ single cell data.
  * Templates can now specify their own plotting functions
  * Use `--flat` and `--interactive` to override this behaviour
* MultiQC added to `bioconda` (with help from @dakl)
* New plugin hook: `config_loaded`
* Plugins can now add new command line options (thanks to @robinandeer)
* Changed default data directory name from `multiqc_report_data` to `multiqc_data`
* Removed support for depreciated MultiQC_OSXApp
* Updated logging so that a verbose `multiqc_data/.multiqc.log` file is always written
* Now logs more stuff in verbose mode - command used, user configs and so on.
* Added a call to multiqc.info to check for new versions. Disable with config `no_version_check`
* Removed general stats manual row sorting.
* Made filename matching use glob unix style filename match patterns
* Everything (including the data directory) is now created in a temporary directory and moved when MultiQC is complete.
* A handful of performance updates for large analysis directories

#### [v0.4](https://github.com/ewels/MultiQC/releases/tag/v0.4) - 2016-02-16
* New `multiqc_sources.txt` which identifies the paths used to collect all report data for each sample
* Export parsed data as tab-delimited text, `JSON` or `YAML` using the new `-k`/`--data-format` command line option
* Updated HighCharts from `v4.2.2` to `v4.2.3`, fixes tooltip hover bug.
* Nicer export button. Now tied to the export toolbox, hopefully more intuitive.
* FastQC: Per base sequence content heatmap can now be clicked to show line graph for single sample
* FastQC: No longer show adapter contamination datasets with <= 0.1% contamination.
* Picard: Added support for `CollectOxoGMetrics` reports.
* Changed command line option `--name` to `--filename`
* `--name` also used for filename if `--filename` not specified.
* Hide samples toolbox now has switch to _show only_ matching samples
* New regex help box with examples added to report
* New button to copy general stats table to the clipboard
* General Stats table 'floating' header now sorts properly when scrolling
* Bugfix: MultiQC default_dev template now copies module assets properly
* Bufgix: General Stats table floating header now resizes properly when page width changes

#### [v0.3.2](https://github.com/ewels/MultiQC/releases/tag/v0.3.2) - 2016-02-08
* All modules now load their log file search parameters from a config
  file, allowing you to overwrite them using your user config file
  * This is useful if your analysis pipeline renames program outputs
* New Picard (sub)modules - Insert Size, GC Bias & HsMetrics
* New Qualimap (sub)module - RNA-Seq QC
* Made Picard MarkDups show percent by default instead of counts
* Added M-Bias plot to Bismark
* New option to stream report HTML to `stdout`
* Files can now be specified as well as directories
* New options to specify whether the parsed data directory should be created
  * command line flags: `--data` / `--no-data`
  * config option name: `make_data_dir`
* Fixed bug with incorrect path to installation dir config YAML file
* New toolbox drawer for bulk-exporting graph images
* Report side navigation can now be hidden to maximise horizontal space
* Mobile styling improved for narrow screen
* More vibrant colours in the general stats table
* General stats table numbers now left aligned
* Settings now saved and loaded to named localstorage locations
  * Simplified interface - no longer global / single report saving
  * Removed static file config. Solves JS error, no-one was doing this
    since we have standalone reports anyway.
* Added support for Python 3.5
* Fixed bug with module specific CSS / JS includes in some templates
* Made the 'ignore files' config use unix style file pattern matching
* Fixed some bugs in the FastQ Screen module
* Fixed some bugs in the FastQC module
* Fixed occasional general stats table bug
* Table sorting on sample names now works after renaming
* Bismark module restructure
  * Each report type now handled independently (alignment / dedup / meth extraction)
  * M-Bias plot now shows R1 and R2
* FastQC GC content plot now has option for counts or percentages
  * Allows comparison between samples with very different read counts
* Bugfix for reports javascript
  * Caused by updated to remotely loaded HighCharts export script
  * Export script now bundled with multiqc, so does not depend on internet connection
  * Other JS errors fixed in this work
* Bugfix for older FastQC reports - handle old style sequence dup data
* Bugfix for varying Tophat alignment report formats
* Bugfix for Qualimap RNA Seq reports with paired end data


#### [v0.3.1](https://github.com/ewels/MultiQC/releases/tag/v0.3.1) - 2015-11-04
* Hotfix patch to fix broken FastQC module (wasn't finding `.zip` files properly)
* General Stats table colours now flat. Should improve browser speed.
* Empty rows now hidden if appear due to column removal in general stats
* FastQC Kmer plot removed until we have something better to show.

#### [v0.3](https://github.com/ewels/MultiQC/releases/tag/v0.3) - 2015-11-04
* Lots of lovely new documentation!
* Child templates - easily customise specific parts of the default report template
* Plugin hooks - allow other tools to execute custom code during MultiQC execution
* New Preseq module
* New design for general statistics table (snazzy new background bars)
* Further development of toolbox
  * New button to clear all filters
  * Warnings when samples are hidden, plus empty plots and table cols are hidden
  * Active toolbar tab buttons are highlighted
* Lots of refactoring by @moonso to please the Pythonic gods
  * Switched to click instead of argparse to handle command line arguments
  * Code generally conforms to best practices better now.
* Now able to supply multiple directories to search for reports
* Logging output improved (now controlled by `-q` and `-v` for quiet and verbose)
* More HTML output dealt with by the base module, less left to the modules
  * Module introduction text
  * General statistics table now much easier to add to (new helper functions)
* Images, CSS and Javascript now included in HTML, meaning that there is a single
  report file to make sharing easier
* More accessible scrolling in the report - styled scrollbars and 'to top' button.
* Modules and templates now use setuptools entry points, facilitating plugins
  by other packages. Allows niche extensions whilst keeping the core codebase clean.
* The general stats table now has a sticky header row when scrolling, thanks to
  some new javascript wizardry...
* General stats columns can have a _shared key_ which allows common colour schemes
  and data ranges. For instance, all columns describing a read count will now share
  their scale across modules.
* General stats columns can be hidden and reordered with a new modal window.
* Plotting code refactored, reports with many samples (>50 by default) don't
  automatically render to avoid freezing the browser.
* Plots with highlighted and renamed samples now honour this when exporting to
  different file types.

#### [v0.2](https://github.com/ewels/MultiQC/releases/tag/v0.2) - 2015-09-18
* Code restructuring for nearly all modules. Common base module
  functions now handle many more functions (plots, config, file import)
  * See the [contributing notes](https://github.com/ewels/MultiQC/blob/master/CONTRIBUTING.md)
    for instructions on how to use these new helpers to make your own module
* New report toolbox - sample highlighting, renaming, hiding
  * Config is autosaved by default, can also export to a file for sharing
  * Interactive tour to help users find their way around
* New Tophat, Bowtie 2 and QualiMap modules
  * Thanks to @guillermo-carrasco for the QualiMap module
* Bowtie module now works
* New command line parameter `-d` prefixes sample names with the directory that
  they were found in. Allows duplicate filenames without being overwritten.
* Introduction walkthrough helps show what can be done in the report
* Now compatible with both Python 2 and Python 3
* Software version number now printed on command line properly, and in reports.
* Bugfix: FastQC doesn't break when only one report found
* Bugfix: FastQC seq content heatmap highlighting
* Many, many small bugfixes

#### [v0.1](https://github.com/ewels/MultiQC/releases/tag/v0.1) - 2015-09-01
* The first public release of MultiQC, after a month of development. Basic
structure in place and modules for FastQC, FastQ Screen, Cutadapt, Bismark, 
STAR, Bowtie, Subread featureCounts and Picard MarkDuplicates. Approaching
stability, though still under fairly heavy development.
