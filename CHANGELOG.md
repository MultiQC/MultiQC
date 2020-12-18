# MultiQC Version History

## MultiQC v1.10dev

#### New MultiQC Features

* `--sample-filters` now also accepts `show_re` and `hide_re` in addition to `show` and `hide`. The `_re` options use regex, while the "normal" options use globbing.
* MultiQC config files now work with `.yml` file extension as well as `.yaml`
  * `.yaml` will take preference if both found.
* Section comments can now also be added for _General Statistics_
  * `section_comments: { general_stats: "My comment" }`
* New CI test looks for git merge markers in files

#### New Modules

* [**HOPS**](https://www.github.com/rhubler/HOPS)
  * Post-alignment ancient DNA analysis tool for MALT
* [**PURPLE**](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator)
  * A purity, ploidy and copy number estimator for whole genome tumor data
* [**Pychopper**](https://github.com/nanoporetech/pychopper)
  * Identify, orient and trim full length Nanopore cDNA reads
* [**qc3C**](https://github.com/cerebis/qc3C)
  * Reference-free QC of Hi-C sequencing data
* [**Sentieon**](https://www.sentieon.com/products/)
  * Submodules added to catch Picard-based QC metrics files
* [**ARTIC**](https://github.com/artic-network/fieldbioinformatics)
  * Viral consensus genome pipeline for nanopore amplicon sequencing data

#### Module updates

* **DRAGEN**
    * Fix issue where missing out fields could crash the module ([#1223](https://github.com/ewels/MultiQC/issues/1223))
* **featureCounts**
    * Add support for output from [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) ([#1022](https://github.com/ewels/MultiQC/issues/1022))
* **hap.py**
    * Updated module to plot both SNP and INDEL stats ([#1241](https://github.com/ewels/MultiQC/issues/1241))
* **Kaiju**
    * Fixed bug affecting inputs with taxa levels other than Phylum ([#1217](https://github.com/ewels/MultiQC/issues/1217))
    * Rework barplot, add top 5 taxons ([#1219](https://github.com/ewels/MultiQC/issues/1219))
* **MALT**
    * Fix y-axis labelling in bargraphs
* **mosdepth**
    * Enable prepending of directory to sample names
* **Picard**
    * Fix `HsMetrics` bait percentage columns ([#1212](https://github.com/ewels/MultiQC/issues/1212))
* **PycoQC**
    * Log10 x-axis for _Read Length_ plot ([#1214](https://github.com/ewels/MultiQC/issues/1214))
* **fgbio**
    * Fix `ErrorRateByReadPosition` to calculate `ymax` not just on the overall `error_rate`, but also specific base errors (ex. `a_to_c_error_rate`, `a_to_g_error_rate`, ...).  ([#1215](https://github.com/ewels/MultiQC/pull/1251))
    * Fix `ErrorRateByReadPosition` plotted line names to no longer concatenate multiple read identifiers and no longer have off-by-one read numbering (ex. `Sample1_R2_R3` -> `Sample1_R2`) ([#[1304](https://github.com/ewels/MultiQC/pull/1304))
* **GATK**
  * Add support for the creation of a "Reported vs Empirical Quality" graph to the Base Recalibration module.
* **Rockhopper**
    * Fix issue with parsing genome names in Rockhopper summary files ([#1333](https://github.com/ewels/MultiQC/issues/1333))
    * Fix issue properly parsing multiple samples within a single Rockhopper summary file

#### New Custom Content features

* General Stats custom content now gives a log message
* If `id` is not set in `JSON` or `YAML` files, it defaults to the sample name instead of just `custom_content`
* Data from `JSON` or `YAML` now has `data` keys (sample names) run through the `clean_s_name()` function to apply sample name cleanup

#### Bug Fixes


## [MultiQC v1.9](https://github.com/ewels/MultiQC/releases/tag/v1.9) - 2020-05-30

#### Dropped official support for Python 2

Python 2 had its [official sunset date](https://www.python.org/doc/sunset-python-2/)
on January 1st 2020, meaning that it will no longer be developed by the Python community.
Part of the [python.org statement](https://www.python.org/doc/sunset-python-2/) reads:

> That means that we will not improve it anymore after that day,
> even if someone finds a security problem in it.
> You should upgrade to Python 3 as soon as you can.

[Very many Python packages no longer support Python 2](https://python3statement.org/)
and it whilst the MultiQC code is currently compatible with both Python 2 and Python 3,
it is increasingly difficult to maintain compatibility with the dependency packages it
uses, such as MatPlotLib, numpy and more.

As of MultiQC version 1.9, **Python 2 is no longer officially supported**.
Automatic CI tests will no longer run with Python 2 and Python 2 specific workarounds
are no longer guaranteed.

Whilst it may be possible to continue using MultiQC with Python 2 for a short time by
pinning dependencies, MultiQC compatibility for Python 2 will now slowly drift and start
to break. If you haven't already, **you need to switch to Python 3 now**.

#### New MultiQC Features

* Now using [GitHub Actions](https://github.com/features/actions) for all CI testing
    * Dropped Travis and AppVeyor, everything is now just on GitHub
    * Still testing on both Linux and Windows, with multiple versions of Python
    * CI tests should now run automatically for anyone who forks the MultiQC repository
* Linting with `--lint` now checks line graphs as well as bar graphs
* New `gathered` template with no tool name sections ([#1119](https://github.com/ewels/MultiQC/issues/1119))
* Added `--sample-filters` option to add _show_/_hide_ buttons at the top of the report ([#1125](https://github.com/ewels/MultiQC/issues/1125))
    * Buttons control the report toolbox Show/Hide tool, filtering your samples
    * Allows reports to be pre-configured based on a supplied list of sample names at report-generation time.
* Line graphs can now have `Log10` buttons (same functionality as bar graphs)
* Importing and running `multiqc` in a script is now a _little_ Better
    * `multiqc.run` now returns the `report` and `config` as well as the exit code. This means that you can explore the MultiQC run time a little in the Python environment.
    * Much more refactoring is needed to make MultiQC as useful in Python scripts as it could be. Watch this space.
* If a custom module `anchor` is set using `module_order`, it's now used a bit more:
    * Prefixed to module _section_ IDs
    * Appended to files saved in `multiqc_data`
    * Should help to prevent duplicates requiring `-1` suffixes when running a module multiple times
* New heatmap plot config options `xcats_samples` and `ycats_samples`
    * If set to `False`, the report toolbox options (_highlight_, _rename_, _show/hide_) do not affect that axis.
    * Means that the _Show only matching samples_ report toolbox option works on FastQC Status Checks, for example ([#1172](https://github.com/ewels/MultiQC/issues/1172))
* Report header time and analysis paths can now be hidden
    * New config options `show_analysis_paths` and `show_analysis_time` ([#1113](https://github.com/ewels/MultiQC/issues/1113))
* New search pattern key `skip: true` to skip specific searches when modules look for a lot of different files (eg. Picard).
* New `--profile-runtime` command line option (`config.profile_runtime`) to give analysis of how long the report takes to be generated
    * Plots of the file search results and durations are added to the end of the MultiQC report as a special module called _Run Time_
    * A summary of the time taken for the major stages of MultiQC execution are printed to the command line log.
* New table config option `only_defined_headers`
    * Defaults to `true`, set to `false` to also show any data columns that are not defined as headers
    * Useful as allows table-wide defaults to be set with column-specific overrides
* New `module` key allowed for `config.extra_fn_clean_exts` and `config.fn_clean_exts`
    * Means you can limit the action of a sample name cleaning pattern to specific MultiQC modules ([#905](https://github.com/ewels/MultiQC/issues/905))

#### New Custom Content features

* Improve support for HTML files - now just end your HTML filename with `_mqc.html`
    * Native handling of HTML snippets as files, no MultiQC config or YAML file required.
    * Also with embedded custom content configuration at the start of the file as a HTML comment.
* Add ability to group custom-content files into report sections
    * Use the new `parent_id`, `parent_name` and `parent_description` config keys to group content together like a regular module ([#1008](https://github.com/ewels/MultiQC/issues/1008))
* Custom Content files can now be configured using `custom_data`, without giving search patterns or data
    * Allows you to set descriptions and nicer titles for images and other 'blunt' data types in reports ([#1026](https://github.com/ewels/MultiQC/issues/1026))
    * Allows configuration of custom content separately from files themselves (`tsv`, `csv`, `txt` formats) ([#1205](https://github.com/ewels/MultiQC/issues/1205))

#### New Modules:

* [**DRAGEN**](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)
    * Illumina Bio-IT Platform that uses FPGA for secondary NGS analysis
* [**iVar**](https://github.com/andersen-lab/ivar)
    * Added support for iVar: a computational package that contains functions broadly useful for viral amplicon-based sequencing.
* [**Kaiju**](http://kaiju.binf.ku.dk/)
    * Fast and sensitive taxonomic classification for metagenomics
* [**Kraken**](https://ccb.jhu.edu/software/kraken2/)
    * K-mer matching tool for taxonomic classification. Module plots bargraph of counts for top-5 hits across each taxa rank. General stats summary.
* [**MALT**](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/malt/)
    * Megan Alignment Tool: Metagenomics alignment tool.
* [**miRTop**](https://github.com/miRTop/mirtop)
    * Command line tool to annotate miRNAs with a standard mirna/isomir naming (mirGFF3)
    * Module started by [@oneillkza](https://github.com/oneillkza/) and completed by [@FlorianThibord](https://github.com/FlorianThibord/)
* [**MultiVCFAnalyzer**](https://github.com/alexherbig/multivcfanalyzer)
    * Combining multiple VCF files into one coherent report and format for downstream analysis.
* **Picard** - new submodules for `QualityByCycleMetrics`, `QualityScoreDistributionMetrics` & `QualityYieldMetrics`
    * See [#1116](https://github.com/ewels/MultiQC/issues/1114)
* [**Rockhopper**](https://cs.wellesley.edu/~btjaden/Rockhopper/)
    * RNA-seq tool for bacteria, includes bar plot showing where features map.
* [**Sickle**](https://github.com/najoshi/sickle)
    * A windowed adaptive trimming tool for FASTQ files using quality
* [**Somalier**](https://github.com/brentp/somalier)
    * Relatedness checking and QC for BAM/CRAM/VCF for cancer, DNA, BS-Seq, exome, etc.
* [**VarScan2**](https://github.com/dkoboldt/varscan)
    * Variant calling and somatic mutation/CNV detection for next-generation sequencing data

#### Module updates:

* **BISCUIT**
    * Major rewrite to work with new BISCUIT QC script (BISCUIT `v0.3.16+`)
        * This change breaks backwards-compatability with previous BISCUIT versions. If you are unable to upgrade BISCUIT, please use MultiQC v1.8.
    * Fixed error when missing data in log files ([#1101](https://github.com/ewels/MultiQC/issues/1101))
* **bcl2fastq**
    * Samples with multiple library preps (i.e barcodes) will now be handled correctly ([#1094](https://github.com/ewels/MultiQC/issues/1094))
* **BUSCO**
    * Updated log search pattern to match new format in v4 with auto-lineage detection option ([#1163](https://github.com/ewels/MultiQC/issues/1163))
* **Cutadapt**
    * New bar plot showing the proportion of reads filtered out for different criteria (eg. _too short_, _too many Ns_) ([#1198](https://github.com/ewels/MultiQC/issues/1198))
* **DamageProfiler**
    * Removes redundant typo in init name. This makes referring to the module's column consistent with other modules when customising general stats table.
* **DeDup**
    * Updates plots to make compatible with 0.12.6
    * Fixes reporting errors - barplot total represents _mapped_ reads, not total reads in BAM file
    * New: Adds 'Post-DeDup Mapped Reads' column to general stats table.
* **FastQC**
    * Fixed tooltip text in _Sequence Duplication Levels_ plot ([#1092](https://github.com/ewels/MultiQC/issues/1092))
    * Handle edge-case where a FastQC report was for an empty file with 0 reads ([#1129](https://github.com/ewels/MultiQC/issues/1129))
* **FastQ Screen**
    * Don't skip plotting `% No Hits` even if it's `0%` ([#1126](https://github.com/ewels/MultiQC/issues/1126))
    * Refactor parsing code. Avoids error with `-0.00 %Unmapped` ([#1126](https://github.com/ewels/MultiQC/issues/1126))
    * New plot for _Bisulfite Reads_, if data is present
    * Categories in main plot are now sorted by the total read count and hidden if 0 across all samples
* **fgbio**
    * New: Plot error rate by read position from `ErrorRateByReadPosition`
    * GroupReadsByUmi plot can now be toggled to show relative percents ([#1147](https://github.com/ewels/MultiQC/pull/1147))
* **FLASh**
    * Logs not reporting innie and outine uncombined pairs now plot combined pairs instead ([#1173](https://github.com/ewels/MultiQC/issues/1173))
* **GATK**
    * Made parsing for VariantEval more tolerant, so that it will work with output from the tool when run in different modes ([#1158](https://github.com/ewels/MultiQC/issues/1158))
* **MTNucRatioCalculator**
    * Fixed misleading value suffix in general stats table
* **Picard MarkDuplicates**
    * **Major change** - previously, if multiple libraries (read-groups) were found then only the first would be used and all others ignored. Now, values from all libraries are merged and `PERCENT_DUPLICATION` and `ESTIMATED_LIBRARY_SIZE` are recalculated. Libraries can be kept as separate samples with a new MultiQC configuration option - `picard_config: markdups_merge_multiple_libraries: False`
    * **Major change** - Updated `MarkDuplicates` bar plot to double the read-pair counts, so that the numbers stack correctly. ([#1142](https://github.com/ewels/MultiQC/issues/1142))
* **Picard HsMetrics**
    * Updated large table to use columns specified in the MultiQC config. See [docs](https://multiqc.info/docs/#hsmetrics). ([#831](https://github.com/ewels/MultiQC/issues/831))
* **Picard WgsMetrics**
    * Updated parsing code to recognise new java class string ([#1114](https://github.com/ewels/MultiQC/issues/1114))
* **QualiMap**
    * Fixed QualiMap mean coverage calculation [#1082](https://github.com/ewels/MultiQC/issues/1082), [#1077](https://github.com/ewels/MultiQC/issues/1082)
* **RSeqC**
    * Support added for output from `geneBodyCoverage2.py` script ([#844](https://github.com/ewels/MultiQC/issues/844))
    * Single sample view in the _"Junction saturation"_ plot now works with the toolbox properly _(rename, hide, highlight)_ ([#1133](https://github.com/ewels/MultiQC/issues/1133))
* **RNASeQC2**
    * Updated to handle the parsing metric files from the [newer rewrite of RNA-SeqQC](https://github.com/broadinstitute/rnaseqc).
* **Samblaster**
    * Improved parsing to handle variable whitespace ([#1176](https://github.com/ewels/MultiQC/issues/1176))
* **Samtools**
    * Removes hardcoding of general stats column names. This allows column names to indicate when a module has been run twice ([https://github.com/ewels/MultiQC/issues/1076](https://github.com/ewels/MultiQC/issues/1076)).
    * Added an observed over expected read count plot for `idxstats` ([#1118](https://github.com/ewels/MultiQC/issues/1118))
    * Added additional (by default hidden) column for `flagstat` that displays number total number of reads in a bam
* **sortmerna**
    * Fix the bug for the latest sortmerna version 4.2.0 ([#1121](https://github.com/ewels/MultiQC/issues/1121))
* **sexdeterrmine**
    * Added a scatter plot of relative X- vs Y-coverage to the generated report.
* **VerifyBAMID**
    * Allow files with column header `FREEMIX(alpha)` ([#1112](https://github.com/ewels/MultiQC/issues/1112))

#### Bug Fixes:

* Added a new test to check that modules work correctly with `--ignore-samples`. A lot of them didn't:
    * `Mosdepth`, `conpair`, `Qualimap BamQC`, `RNA-SeQC`, `GATK BaseRecalibrator`, `SNPsplit`, `SeqyClean`, `Jellyfish`, `hap.py`, `HOMER`, `BBMap`, `DeepTools`, `HiCExplorer`, `pycoQC`, `interop`
    * These modules have now all been fixed and `--ignore-samples` should work as you expect for whatever data you have.
* Removed use of `shutil.copy` to avoid problems with working on multiple filesystems ([#1130](https://github.com/ewels/MultiQC/issues/1130))
* Made folder naming behaviour of `multiqc_plots` consistent with `multiqc_data`
    * Incremental numeric suffixes now added if folder already exists
    * Plots folder properly renamed if using `-n`/`--filename`
* Heatmap plotting function is now compatible with MultiQC toolbox `hide` and `highlight` ([#1136](https://github.com/ewels/MultiQC/issues/1136))
* Plot config `logswitch_active` now works as advertised
* When running MultiQC modules several times, multiple data files are now created instead of overwriting one another ([#1175](https://github.com/ewels/MultiQC/issues/1175))
* Fixed minor bug where tables could report negative numbers of columns in their header text
* Fixed bug where numeric custom content sample names could trigger a `TypeError` ([#1091](https://github.com/ewels/MultiQC/issues/1091))
* Fixed custom content bug HTML data in a config file would trigger a `ValueError` ([#1071](https://github.com/ewels/MultiQC/issues/1071))
* Replaced deprecated 'warn()' with 'warning()' of the logging module
* Custom content now supports `section_extra` config key to add custom HTML after description.
* Barplots with `ymax` set now ignore this when you click the _Percentages_ tab.

## [MultiQC v1.8](https://github.com/ewels/MultiQC/releases/tag/v1.8) - 2019-11-20

#### New Modules:

* [**fgbio**](http://fulcrumgenomics.github.io/fgbio/)
    * Process family size count hist data from `GroupReadsByUmi`
* [**biobambam2**](https://github.com/gt1/biobambam2)
    * Added submodule for `bamsormadup` tool
    * Totally cheating - it uses Picard MarkDuplicates but with a custom search pattern and naming
* [**SeqyClean**](https://github.com/ibest/seqyclean)
    * Adds analysis for seqyclean files
* [**mtnucratio**](https://github.com/apeltzer/MTNucRatioCalculator)
    * Added little helper tool to compute mt to nuclear ratios for NGS data.
* [**mosdepth**](https://github.com/brentp/mosdepth)
    * fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing
* [**SexDetErrmine**](https://github.com/TCLamnidis/Sex.DetERRmine)
    * Relative coverage and error rate of X and Y chromosomes
* [**SNPsplit**](https://github.com/FelixKrueger/SNPsplit)
    * Allele-specific alignment sorting

#### Module updates:
* **bcl2fastq**
    * Added handling of demultiplexing of more than 2 reads
    * Allow bcl2fastq to parse undetermined barcode information in situations when lane indexes do not start at 1
* **BBMap**
    * Support for scafstats output marked as not yet implemented in docs
* **DeDup**
    * Added handling clusterfactor and JSON logfiles
* **damageprofiler**
    * Added writing metrics to data output file.
* **DeepTools**
    * Fixed Python3 bug with int() conversion ([#1057](https://github.com/ewels/MultiQC/issues/1057))
    * Handle varied TES boundary labels in plotProfile ([#1011](https://github.com/ewels/MultiQC/issues/1011))
    * Fixed bug that prevented running on only plotProfile files when no other deepTools files found.
* **fastp**
    * Fix faulty column handling for the _after filtering_ Q30 rate ([#936](https://github.com/ewels/MultiQC/issues/936))
* **FastQC**
    * When including a FastQC section multiple times in one report, the Per Base Sequence Content heatmaps now behave as you would expect.
    * Added heatmap showing FastQC status checks for every section report across all samples
    * Made sequence content individual plots work after samples have been renamed ([#777](https://github.com/ewels/MultiQC/issues/777))
    * Highlighting samples from status - respect chosen highlight colour in the toolbox ([#742](https://github.com/ewels/MultiQC/issues/742))
* **FastQ Screen**
    * When including a FastQ Screen section multiple times in one report, the plots now behave as you would expect.
    * Fixed MultiQC linting errors
* **fgbio**
    * Support the new output format of `ErrorRateByReadPosition` first introduced in version `1.3.0`, as well as the old output format.
* **GATK**
    * Refactored BaseRecalibrator code to be more consistent with MultiQC Python style
    * Handle zero count errors in BaseRecalibrator
* **HiC Explorer**
    * Fixed bug where module tries to parse `QC_table.txt`, a new log file in hicexplorer v2.2.
    * Updated the format of the report to fits the changes which have been applied to the QC report of hicexplorer v3.3
    * Updated code to save parsed results to `multiqc_data`
* **HTSeq**
    * Fixed bug where module would crash if a sample had zero reads ([#1006](https://github.com/ewels/MultiQC/issues/1006))
* **LongRanger**
    * Added support for the LongRanger Align pipeline.
* **miRTrace**
    * Fixed bug where a sample in some plots was missed. ([#932](https://github.com/ewels/MultiQC/issues/932))
* **Peddy**
    * Fixed bug where sample name cleaning could lead to error. ([#1024](https://github.com/ewels/MultiQC/issues/1024))
    * All plots (including _Het Check_ and _Sex Check_) now hidden if no data
* **Picard**
    * Modified `OxoGMetrics` so that it will find files created with GATK `CollectMultipleMetrics` and `ConvertSequencingArtifactToOxoG`.
* **QoRTs**
    * Fixed bug where `--dirs` broke certain input files. ([#821](https://github.com/ewels/MultiQC/issues/821))
* **Qualimap**
    * Added in mean coverage computation for general statistics report
    * Creates now tables of collected data in `multiqc_data`
* **RNA-SeQC**
    * Updated broken URL link
* **RSeQC**
    * Fixed bug where Junction Saturation plot when clicking a single sample was mislabelling the lines.
    * When including a RSeQC section multiple times in one report, clicking Junction Saturation plot now behaves as you would expect.
    * Fixed bug where exported data in `multiqc_rseqc_read_distribution.txt` files had incorrect values for `_kb` fields ([#1017](https://github.com/ewels/MultiQC/issues/1017))
* **Samtools**
    * Utilize in-built `read_count_multiplier` functionality to plot `flagstat` results more nicely
* **SnpEff**
    * Increased the default summary csv file-size limit from 1MB to 5MB.
* **Stacks**
    * Fixed bug where multi-population sum stats are parsed correctly ([#906](https://github.com/ewels/MultiQC/issues/906))
* **TopHat**
    * Fixed bug where TopHat would try to run with files from Bowtie2 or HiSAT2 and crash
* **VCFTools**
    * Fixed a bug where `tstv_by_qual.py` produced invalid json from infinity-values.
* **snpEff**
    * Added plot of effects

#### New MultiQC Features:
* Added some installation docs for windows
* Added some docs about using MultiQC in bioinformatics pipelines
* Rewrote Docker image
    * New base image `czentye/matplotlib-minimal` reduces image size from ~200MB to ~80MB
    * Proper installation method ensures latest version of the code
    * New entrypoint allows easier command-line usage
* Support opening MultiQC on websites with CSP `script-src 'self'` with some sha256 exceptions
    * Plot data is no longer intertwined with javascript code so hashes stay the same
* Made `config.report_section_order` work for module sub-sections as well as just modules.
* New config options `exclude_modules` and `run_modules` to complement `-e` and `-m` cli flags.
* Command line output is now coloured by default :rainbow: (use `--no-ansi` to turn this off)
* Better launch comparability due to code refactoring by [@KerstenBreuer](https://github.com/KerstenBreuer) and [@ewels](https://github.com/ewels)
    * Windows support for base `multiqc` command
    * Support for running as a python module: `python -m multiqc .`
    * Support for running within a script: `import multiqc` and `multiqc.run('/path/to/files')`
* Config option `custom_plot_config` now works for bargraph category configs as well ([#1044](https://github.com/ewels/MultiQC/issues/1044))
* Config `table_columns_visible` can now be given a module namespace and it will hide all columns from that module ([#541](https://github.com/ewels/MultiQC/issues/541))

#### Bug Fixes:
* MultiQC now ignores all `.md5` files
* Use `SafeLoader` for PyYaml load calls, avoiding recent warning messages.
* Hide `multiqc_config_example.yaml` in the `test` directory to stop people from using it without modification.
* Fixed matplotlib background colour issue (@epakarin - [#886](https://github.com/ewels/MultiQC/issues))
* Table rows that are empty due to hidden columns are now properly hidden on page load ([#835](https://github.com/ewels/MultiQC/issues/835))
* Sample name cleaning: All sample names are now truncated to their basename, without a path.
  * This includes for `regex` and `replace` (before was only the default `truncate`).
  * Only affects modules that take sample names from file contents, such as cutadapt.
  * See [#897](https://github.com/ewels/MultiQC/issues/897) for discussion.




## [MultiQC v1.7](https://github.com/ewels/MultiQC/releases/tag/v1.7) - 2018-12-21

#### New Modules:
* [**BISCUIT**](https://github.com/zwdzwd/biscuit)
    * BISuilfite-seq CUI Toolkit
    * Module written by [@zwdzwd](https://github.com/zwdzwd/)
* [**DamageProfiler**](https://github.com/Integrative-Transcriptomics/DamageProfiler)
    * A tool to determine ancient DNA misincorporation rates.
    * Module written by [@apeltzer](https://github.com/apeltzer/)
* [**FLASh**](https://ccb.jhu.edu/software/FLASH/)
    * FLASH (Fast Length Adjustment of SHort reads)
    * Module written by [@pooranis](https://github.com/pooranis/)
* [**MinIONQC**](https://github.com/roblanf/minion_qc)
    * QC of reads from ONT long-read sequencing
    * Module written by [@ManavalanG](https://github.com/ManavalanG)
* [**phantompeakqualtools**](https://www.encodeproject.org/software/phantompeakqualtools)
    * A tool for informative enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data.
    * Module written by [@chuan-wang](https://github.com/chuan-wang/)
* [**Stacks**](http://catchenlab.life.illinois.edu/stacks/)
    * A software for analyzing restriction enzyme-based data (e.g. RAD-seq). Support for Stacks >= 2.1 only.
    * Module written by [@remiolsen](https://github.com/remiolsen/)

#### Module updates:
* **AdapterRemoval**
    * Handle error when zero bases are trimmed. See [#838](https://github.com/ewels/MultiQC/issues/838).
* **Bcl2fastq**
    * New plot showing the top twenty of undetermined barcodes by lane.
    * Informations for R1/R2 are now separated in the General Statistics table.
    * SampleID is concatenate with SampleName because in Chromium experiments several sample have the same SampleName.
* **deepTools**
    * New PCA plots from the `plotPCA` function (written by [@chuan-wang](https://github.com/chuan-wang/))
    * New fragment size distribution plots from `bamPEFragmentSize --outRawFragmentLengths` (written by [@chuan-wang](https://github.com/chuan-wang/))
    * New correlation heatmaps from the `plotCorrelation` function (written by [@chuan-wang](https://github.com/chuan-wang/))
    * New sequence distribution profiles around genes, from the `plotProfile` function (written by [@chuan-wang](https://github.com/chuan-wang/))
    * Reordered sections
* **Fastp**
    * Fixed bug in parsing of empty histogram data. See [#845](https://github.com/ewels/MultiQC/issues/845).
* **FastQC**
    * Refactored _Per Base Sequence Content_ plots to show original underlying data, instead of calculating it from the page contents. Now shows original FastQC base-ranges and fixes 100% GC bug in final few pixels. See [#812](https://github.com/ewels/MultiQC/issues/812).
    * When including a FastQC section multiple times in one report, the summary progress bars now behave as you would expect.
* **FastQ Screen**
    * Don't hide genomes in the simple plot, even if they have zero unique hits. See [#829](https://github.com/ewels/MultiQC/issues/829).
* **InterOp**
    * Fixed bug where read counts and base pair yields were not displaying in tables correctly.
    * Number formatting for these fields can now be customised in the same way as with other modules, as described [in the docs](http://multiqc.info/docs/#number-base-multiplier)
* **Picard**
    * InsertSizeMetrics: You can now configure to what degree the insert size plot should be smoothed.
    * CollectRnaSeqMetrics: Add warning about missing rRNA annotation.
    * CollectRnaSeqMetrics: Add chart for counts/percentage of reads mapped to the correct strand.
    * Now parses VariantCallingMetrics reports. (Similar to GATK module's VariantEval.)
* **phantompeakqualtools**
    * Properly clean sample names
* **Trimmomatic**
    * Updated Trimmomatic module documentation to be more helpful
    * New option to use filenames instead of relying on the command line used. See [#864](https://github.com/ewels/MultiQC/issues/864).

#### New MultiQC Features:
* Embed your custom images with a new Custom Content feature! Just add `_mqc` to the end of the filename for `.png`, `.jpg` or `.jpeg` files.
* Documentation for Custom Content reordered to make it a little more sane
* You can now add or override any config parameter for any MultiQC plot! See [the documentation](http://multiqc.info/docs/#customising-plots) for more info.
* Allow `table_columns_placement` config to work with table IDs as well as column namespaces. See [#841](https://github.com/ewels/MultiQC/issues/841).
* Improved visual spacing between grouped bar plots


#### Bug Fixes:
* Custom content no longer clobbers `col1_header` table configs
* The option `--file-list` that refers to a text file with file paths to analyse will no longer ignore directory paths
* [Sample name directory prefixes](https://multiqc.info/docs/#sample-names-prefixed-with-directories) are now added _after_ cleanup.
* If a module is run multiple times in one report, it's CSS and JS files will only be included once (`default` template)




## [MultiQC v1.6](https://github.com/ewels/MultiQC/releases/tag/v1.6) - 2018-08-04

Some of these updates are thanks to the efforts of people who attended the [NASPM](https://twitter.com/NordicGenomics) 2018 MultiQC hackathon session. Thanks to everyone who attended!

#### New Modules:
* [**fastp**](https://github.com/OpenGene/fastp)
    * An ultra-fast all-in-one FASTQ preprocessor (QC, adapters, trimming, filtering, splitting...)
    * Module started by [@florianduclot](https://github.com/florianduclot/) and completed by [@ewels](https://github.com/ewels/)
* [**hap.py**](https://github.com/Illumina/hap.py)
    * Hap.py is a set of programs based on htslib to benchmark variant calls against gold standard truth datasets
    * Module written by [@tsnowlan](https://github.com/tsnowlan/)
* [**Long Ranger**](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger)
    * Works with data from the 10X Genomics Chromium. Performs sample demultiplexing, barcode processing, alignment, quality control, variant calling, phasing, and structural variant calling.
    * Module written by [@remiolsen](https://github.com/remiolsen/)
* [**miRTrace**](https://github.com/friedlanderlab/mirtrace)
    * A quality control software for small RNA sequencing data.
    * Module written by [@chuan-wang](https://github.com/chuan-wang/)


#### Module updates:
* **BCFtools**
    * New plot showing SNP statistics versus quality of call from bcftools stats ([@MaxUlysse](https://github.com/MaxUlysse) and [@Rotholandus](https://github.com/Rotholandus))
* **BBMap**
    * Support added for BBDuk kmer-based adapter/contaminant filtering summary stats ([@boulund](https://github.com/boulund)
* **FastQC**
    * New read count plot, split into unique and duplicate reads if possible.
    * Help text added for all sections, mostly copied from the excellent FastQC help.
    * Sequence duplication plot rescaled
* **FastQ Screen**
    * Samples in large-sample-number plot are now sorted alphabetically ([@hassanfa](https://github.com/hassanfa)
* **MACS2**
    * Output is now more tolerant of missing data (no plot if no data)
* **Peddy**
    * Background samples now shown in ancestry PCA plot ([@roryk](https://github.com/roryk))
    * New plot showing sex checks versus het ratios, supporting unknowns ([@oyvinev](https://github.com/oyvinev))
* **Picard**
    * New submodule to handle `ValidateSamFile` reports ([@cpavanrun](https://github.com/cpavanrun))
    * WGSMetrics now add the mean and standard-deviation coverage to the general stats table (hidden) ([@cpavanrun](https://github.com/cpavanrun))
* **Preseq**
    * New config option to plot preseq plots with unique old coverage on the y axis instead of read count
    * Code refactoring by [@vladsaveliev](https://github.com/vladsaveliev)
* **QUAST**
    * Null values (`-`) in reports now handled properly. Bargraphs always shown despite varying thresholds. ([@vladsaveliev](https://github.com/vladsaveliev))
* **RNA-SeQC**
    * Don't create the report section for Gene Body Coverage if no data is given
* **Samtools**
    * Fixed edge case bug where MultiQC could crash if a sample had zero count coverage with idxstats.
    * Adds % proper pairs to general stats table
* **Skewer**
    * Read length plot rescaled
* **Tophat**
    * Fixed bug where some samples could be given a blank sample name ([@lparsons](https://github.com/lparsons))
* **VerifyBamID**
    * Change column header help text for contamination to match percentage output ([@chapmanb](https://github.com/chapmanb))

#### New MultiQC Features:
* New config option `remove_sections` to skip specific report sections from modules
* Add `path_filters_exclude` to exclude certain files when running modules multiple times. You could previously only include certain files.
* New `exclude_*` keys for file search patterns
    * Have a subset of patterns to exclude otherwise detected files with, by filename or contents
* Command line options all now use mid-word hyphens (not a mix of hyphens and underscores)
    * Old underscore terms still maintained for backwards compatibility
* Flag `--view-tags` now works without requiring an "analysis directory".
* Removed Python dependency for `enum34` ([@boulund](https://github.com/boulund))
* Columns can be added to `General Stats` table for custom content/module.
* New `--ignore-symlinks` flag which will ignore symlinked directories and files.
* New `--no-megaqc-upload` flag which disables automatically uploading data to MegaQC

#### Bug Fixes
* Fix path_filters for top_modules/module_order configuration only selecting if *all* globs match. It now filters searches that match *any* glob.
* Empty sample names from cleaning are now no longer allowed
* Stop prepend_dirs set in the config from getting clobbered by an unpassed CLI option ([@tsnowlan](https://github.com/tsnowlan))
* Modules running multiple times now have multiple sets of columns in the General Statistics table again, instead of overwriting one another.
* Prevent tables from clobbering sorted row orders.
* Fix linegraph and scatter plots data conversion (sporadically the incorrect `ymax` was used to drop data points) ([@cpavanrun](https://github.com/cpavanrun))
* Adjusted behavior of ceiling and floor axis limits
* Adjusted multiple file search patterns to make them more specific
    * Prevents the wrong module from accidentally slurping up output from a different tool. By [@cpavanrun](https://github.com/cpavanrun) (see [PR #727](https://github.com/ewels/MultiQC/pull/727))
* Fixed broken report bar plots when `-p`/`--export-plots` was specified (see issue [#801](https://github.com/ewels/MultiQC/issues/801))





## [MultiQC v1.5](https://github.com/ewels/MultiQC/releases/tag/v1.5) - 2018-03-15

#### New Modules:
* [**HiCPro**](https://github.com/nservant/HiC-Pro) - New module!
    * HiCPro: Quality controls and processing of Hi-C
    * Module written by [@nservant](https://github.com/nservant),
* [**DeDup**](http://www.github.com/apeltzer/DeDup) - New module!
    * DeDup: Improved Duplicate Removal for merged/collapsed reads in ancient DNA analysis
    * Module written by [@apeltzer](https://github.com/apeltzer),
* [**Clip&Merge**](http://github.com/apeltzer/ClipAndMerge) - New module!
    * Clip&Merge: Adapter clipping and read merging for ancient DNA analysis
    * Module written by [@apeltzer](https://github.com/apeltzer),

#### Module updates:
* **bcl2fastq**
    * Catch `ZeroDivisionError` exceptions when there are 0 reads ([@aledj2](https://github.com/aledj2))
    * Add parsing of `TrimmedBases` and new General Stats column for % bases trimmed ([@matthdsm](https://github.com/matthdsm)).
* **BUSCO**
    * Fixed configuration bug that made all sample names become `'short'`
* **Custom Content**
    * Parsed tables now exported to `multiqc_data` files
* **Cutadapt**
    * Refactor parsing code to collect all length trimming plots
* **FastQC**
    * Fixed starting y-axis label for GC-content lineplot being incorrect.
* **HiCExplorer**
    * Updated to work with v2.0 release.
* **Homer**
    * Made parsing of `tagInfo.txt` file more resilient to variations in file format so that it works with new versions of Homer.
    * Kept order of chromosomes in coverage plot consistent.
* **Peddy**
    * Switch `Sex error` logic to `Correct sex` for better highlighting ([@aledj2](https://github.com/aledj2))
* **Picard**
    * Updated module and search patterns to recognise new output format from Picard version >= 2.16 and GATK output.
* **Qualimap BamQC**
    * Fixed bug where start of _Genome Fraction_ could have a step if target is 100% covered.
* **RNA-SeQC**
    * Added rRNA alignment stats to summary table [@Rolandde](https://github.com/Rolandde)
* **RSeqC**
    * Fixed read distribution plot by adding category for `other_intergenic` (thanks to [@moxgreen](https://github.com/moxgreen))
    * Fixed a dodgy plot title (Read GC content)
* **Supernova**
    * Added support for Supernova 2.0 reports. Fixed a TypeError bug when using txt reports only. Also a bug when parsing empty histogram files.

#### New MultiQC Features:
* Invalid choices for `--module` or `--exclude` now list the available modules alphabetically.
* Linting now checks for presence in `config.module_order` and tags.

#### Bug Fixes
* Excluding modules now works in combination with using module tags.
* Fixed edge-case bug where certain combinations of `output_fn_name` and `data_dir_name` could trigger a crash
* Conditional formatting - values are now longer double-labelled
* Made config option `extra_series` work in scatter plots the same way that it works for line plots
* Locked the `matplotlib` version to `v2.1.0` and below
    * Due to [two](https://github.com/matplotlib/matplotlib/issues/10476) [bugs](https://github.com/matplotlib/matplotlib/issues/10784) that appeared in `v2.2.0` - will remove this constraint when there's a new release that works again.





## [MultiQC v1.4](https://github.com/ewels/MultiQC/releases/tag/v1.4) - 2018-01-11

A slightly earlier-than-expected release due to a new problem with dependency packages that is breaking MultiQC installations since 2018-01-11.

#### New Modules:
* [**Sargasso**](http://statbio.github.io/Sargasso/)
    * Parses output from Sargasso - a tool to separate mixed-species RNA-seq reads according to their species of origin
    * Module written by [@hxin](https://github.com/hxin/)
* [**VerifyBAMID**](https://genome.sph.umich.edu/wiki/VerifyBamID)
    * Parses output from VerifyBAMID - a tool to detect contamination in BAM files.
    * Adds the `CHIPMIX` and `FREEMIX` columns to the general statistics table.
    * Module written by [@aledj2](https://github.com/aledj2/)

#### Module updates:
* **MACS2**
    * Updated to work with output from older versions of MACS2 by [@avilella](https://github.com/avilella/)
* **Peddy**
    * Add het check plot to suggest potential contamination by [@aledj2](https://github.com/aledj2)
* **Picard**
    * Picard HsMetrics `HS_PENALTY` plot now has correct axis labels
    * InsertSizeMetrics switches commas for points if it can't convert floats. Should help some european users.
* **QoRTs**
    * Added support for new style of output generated in the v1.3.0 release
* **Qualimap**
    * New `Error rate` column in General Statistics table, added by [@Cashalow](https://github.com/Cashalow/)
        * Hidden by default - customise your MultiQC config to always show this column (see [docs](http://multiqc.info/docs/#hiding-columns))
* **QUAST**
    * New option to customise the default display of contig count and length (eg. `bp` instead of `Mbp`).
    * See [documentation](http://multiqc.info/docs/#quast). Written by [@ewels](https://github.com/ewels/) and [@Cashalow](https://github.com/Cashalow/)
* **RSeQC**
    * Removed normalisation in Junction Saturation plot. Now raw counts instead of % of total junctions.

#### New MultiQC Features:
* Conditional formatting / highlighting of cell contents in tables
    * If you want to make values that match a criteria stand out more, you can now write custom rules and formatting instructions for tables.
    * For instructions, see [the documentation](http://multiqc.info/docs/#conditional-formatting)
* New `--lint` option which is strict about best-practices for writing new modules
    * Useful when writing new modules and code as it throws warnings
    * Currently only implemented for bar plots and a few other places. More linting coming soon...
* If MultiQC breaks and shows am error message, it now reports the filename of the last log it found
    * Hopefully this will help with debugging / finding dodgy input data

#### Bug Fixes
* Addressed new dependency error with conflicting package requirements
    * There was a conflict between the `networkx`, `colormath` and `spectra` releases.
    * I previously forced certain software versions to get around this, but `spectra` has now updated with the unfortunate effect of introducing a new dependency clash that halts installation.
* Fixed newly introduced bug where Custom Content MultiQC config file search patterns had been broken
* Updated pandoc command used in `--pdf` to work with new releases of Pandoc
* Made config `table_columns_visible` module name key matching case insensitive to make less frustrating





## [MultiQC v1.3](https://github.com/ewels/MultiQC/releases/tag/v1.3) - 2017-11-03

#### Breaking changes - custom search patterns
Only for users with custom search patterns for the `bowtie` or `star`: you will
need to update your config files - the `bowtie` search key is now `bowtie1`,
`star_genecounts` is now `star/genecounts`.

For users with custom modules - search patterns _must_ now conform to the search
pattern naming convention: `modulename` or `modulename/anything` (the search pattern
string beginning with the name of your module, anything you like after the first `/`).

#### New Modules:
* [**10X Supernova**](https://support.10xgenomics.com/de-novo-assembly/software/overview/welcome)
    * Parses statistics from the _de-novo_ Supernova software.
    * Module written by [@remiolsen](https://github.com/remiolsen/)
* [**BBMap**](https://sourceforge.net/projects/bbmap/)
    * Plot metrics from a number of BBMap tools, a suite of DNA/RNA mapping tools and utilities
    * Module written by [@boulund](https://github.com/boulund/) and [@epruesse](https://github.com/epruesse/)
* [**deepTools**](https://github.com/fidelram/deepTools) - new module!
    * Parse text output from `bamPEFragmentSize`, `estimateReadFiltering`, `plotCoverage`, `plotEnrichment`, and `plotFingerprint`
    * Module written by [@dpryan79](https://github.com/dpryan79/)
* [**Homer Tag Directory**](http://homer.ucsd.edu/homer/ngs/tagDir.html) - new submodule!
    * Module written by [@rdali](https://github.com/rdali/)
* [**illumina InterOp**](http://illumina.github.io/interop/index.html)
    * Module to parse metrics from illumina sequencing runs and demultiplexing, generated by the InterOp package
    * Module written by [@matthdsm](https://github.com/matthdsm/)
* [**RSEM**](https://deweylab.github.io/RSEM/) - new module!
    * Parse `.cnt` file comming from rsem-calculate-expression and plot read repartitions (Unalignable, Unique, Multi ...)
    * Module written by [@noirot](https://github.com/noirot/)
* [**HiCExplorer**](https://github.com/maxplanck-ie/HiCExplorer)
    * New module to parse the log files of `hicBuildMatrix`.
    * Module written by [@joachimwolff](https://github.com/joachimwolff/)

#### Module updates:
* **AfterQC**
    * Handle new output format where JSON summary key changed names.
* **bcl2fastq**
    * Clusters per sample plot now has tab where counts are categoried by lane.
* **GATK**
    * New submodule to handle Base Recalibrator stats, written by [@winni2k](https://github.com/winni2k/)
* **HiSAT2**
    * Fixed bug where plot title was incorrect if both SE and PE bargraphs were in one report
* **Picard HsMetrics**
    * Parsing code can now handle commas for decimal places
* **Preseq**
    * Updated odd file-search pattern that limited input files to 500kb
* **QoRTs**
    * Added new plots, new helptext and updated the module to produce a lot more output.
* **Qualimap BamQC**
    * Fixed edge-case bug where the refactored coverage plot code could raise an error from the `range` call.
* Documentation and link fixes for Slamdunk, GATK, bcl2fastq, Adapter Removal, FastQC and main docs
    * Many of these spotted and fixed by [@juliangehring](https://github.com/juliangehring/)
* Went through all modules and standardised plot titles
    * All plots should now have a title with the format _Module name: Plot name_

#### New MultiQC Features:
* New MultiQC docker image
    * Ready to use docker image now available at https://hub.docker.com/r/ewels/multiqc/ (200 MB)
    * Uses automated builds - pull `:latest` to get the development version, future releases will have stable tags.
    * Written by [@MaxUlysse](https://github.com/MaxUlysse/)
* New `module_order` config options allow modules to be run multiple times
    * Filters mean that a module can be run twice with different sets of files (eg. before and after trimming)
    * Custom module config parameters can be passed to module for each run
* File search refactored to only search for running modules
    * Makes search much faster when running with lots of files and limited modules
    * For example, if using `-m star` to only use the STAR module, all other file searches now skipped
* File search now warns if an unrecognised search type is given
* MultiQC now saves nearly all parsed data to a structured output file by default
    * See `multiqc_data/multiqc_data.json`
    * This can be turned off by setting `config.data_dump_file: false`
* Verbose logging when no log files found standardised. Less duplication in code and logs easier to read!
* New documentation section describing how to use MultiQC with Galaxy
* Using `shared_key: 'read_counts'` in table header configs now applies relevant defaults

#### Bug Fixes
* Installation problem caused by changes in upstream dependencies solved by stricter installation requirements
* Minor `default_dev` directory creation bug squashed
* Don't prepend the directory separator (`|`) to sample names with `-d` when there are no subdirs
* `yPlotLines` now works even if you don't set `width`





## [MultiQC v1.2](https://github.com/ewels/MultiQC/releases/tag/v1.2) - 2017-08-16

#### CodeFest 2017 Contributions
We had a fantastic group effort on MultiQC at the [2017 BOSC CodeFest](https://www.open-bio.org/wiki/Codefest_2017).
Many thanks to those involved!

#### New Modules:
* [**AfterQC**](https://github.com/OpenGene/AfterQC) - New module!
    * Added parsing of the _AfterQC_ json file data, with a plot of filtered reads.
    * Work by [@raonyguimaraes](https://github.com/raonyguimaraes)
* [**bcl2fastq**](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
    * bcl2fastq can be used to both demultiplex data and convert BCL files to FASTQ file formats for downstream analysis
    * New module parses JSON output from recent versions and summarises some key statistics from the demultiplexing process.
    * Work by [@iimog](https://github.com/iimog) (with a little help from [@tbooth](https://github.com/tbooth) and [@ewels](https://github.com/ewels))
* [**leeHom**](https://github.com/grenaud/leeHom)
    * leeHom is a program for the Bayesian reconstruction of ancient DNA
* [**VCFTools**](https://vcftools.github.io)
    * Added initial support for VCFTools `relatedness2`
    * Added support for VCFTools `TsTv-by-count` `TsTv-by-qual` `TsTv-summary`
    * Module written by [@mwhamgenomics](https://github.com/mwhamgenomics)

#### Module updates:
* **FastQ Screen**
    * Gracefully handle missing data from very old FastQ Screen versions.
* **RNA-SeQC**
    * Add new transcript-associated reads plot.
* **Picard**
    * New submodule to handle output from `TargetedPcrMetrics`
* **Prokka**
    * Added parsing of the `# CRISPR arrays` data from Prokka when available ([@asetGem](https://github.com/asetGem))
* **Qualimap**
    * Some code refactoring to radically improve performance and run times, especially with high coverage datasets.
    * Fixed bug where _Cumulative coverage genome fraction_ plot could be truncated.

#### New MultiQC Features:
* New module help text
    * Lots of additional help text was written to make MultiQC report plots easier to interpret.
    * Updated modules:
        * Bowtie
        * Bowtie 2
        * Prokka
        * Qualimap
        * SnpEff
    * Elite team of help-writers:
        * [@tabwalsh](https://github.com/tabwalsh)
        * [@ddesvillechabrol](https://github.com/tabwalsh)
        * [@asetGem](https://github.com/asetGem)
* New config option `section_comments` allows you to add custom comments above specific sections in the report
* New `--tags` and `--view_tags` command line options
    * Modules can now be given tags (keywords) and filtered by those. So running `--tags RNA` will only run MultiQC modules related to RNA analysis.
    * Work by [@Hammarn](https://github.com/Hammarn)
* Back-end configuration options to specify the order of table columns
    * Modules and user configs can set priorities for columns to customise where they are displayed
    * Work by [@tbooth](https://github.com/tbooth)
* Added framework for proper unit testing
    * Previous start on unit tests tidied up, new blank template and tests for the `clean_sample_name` functionality.
    * Added to Travis and Appveyor for continuous integration testing.
    * Work by [@tbooth](https://github.com/tbooth)
* Bug fixes and refactoring of report configuration saving / loading
    * Discovered and fixed a bug where a report config could only be loaded once
    * Work by [@DennisSchwartz](https://github.com/DennisSchwartz)
* Table column row headers (sample names) can now be numeric-only.
    * Work by [@iimog](https://github.com/iimog)
* Improved sample name cleaning functionality
    * Added option `regex_keep` to clean filenames by _keeping_ the matching part of a pattern
    * Work by [@robinandeer](https://github.com/robinandeer)
* Handle error when invalid regexes are given in reports
    * Now have a nice toast error warning you and the invalid regexes are highlighted
    * Previously this just crashed the whole report without any warning
    * Work by [@robinandeer](https://github.com/robinandeer)
* Command line option `--dirs-depth` now sets `-d` to `True` (so now works even if `-d` isn't also specified).
* New config option `config.data_dump_file` to export as much data as possible to `multiqc_data/multiqc_data.json`
* New code to send exported JSON data to a a web server
    * This is in preparation for the upcoming MegaQC project. Stay tuned!

#### Bug Fixes:
* Specifying multiple config files with `-c`/`--config` now works as expected
    * Previously this would only read the last specified
* Fixed table rendering bug that affected Chrome v60 and IE7-11
    * Table cell background bars weren't showing up. Updated CSS to get around this rendering error.
* HTML ID cleanup now properly cleans strings so that they work with jQuery as expected.
* Made bar graph sample highlighting work properly again
* Config `custom_logo` paths can now be relative to the config file (or absolute as before)
* Report doesn't keep annoyingly telling you that toolbox changes haven't been applied
    * Now uses more subtle _toasts_ and only when you close the toolbox (not every click).
* Switching report toolbox options to regex mode now enables the _Apply_ button as it should.
* Sorting table columns with certain suffixes (eg. `13X`) no works properly (numerically)
* Fixed minor bug in line plot data smoothing (now works with unsorted keys)

---

## [MultiQC v1.1](https://github.com/ewels/MultiQC/releases/tag/v1.1) - 2017-07-18

#### New Modules:

* [**BioBloom Tools**](https://github.com/bcgsc/biobloom)
    * Create Bloom filters for a given reference and then to categorize sequences
* [**Conpair**](https://github.com/nygenome/Conpair)
    * Concordance and contamination estimator for tumor–normal pairs
* [**Disambiguate**](https://github.com/AstraZeneca-NGS/disambiguate)
    * Bargraph displaying the percentage of reads aligning to two different reference genomes.
* [**Flexbar**](https://github.com/seqan/flexbar)
    * Flexbar is a tool for flexible barcode and adapter removal.
* [**HISAT2**](https://ccb.jhu.edu/software/hisat2/)
    * New module for the HISAT2 aligner.
    * Made possible by updates to HISAT2 logging by @infphilo (requires `--new-summary` HISAT2 flag).
* [**HOMER**](http://homer.ucsd.edu/homer/)
    * Support for summary statistics from the `findPeaks` tool.
* [**Jellyfish**](http://www.cbcb.umd.edu/software/jellyfish/)
    * Histograms to estimate library complexity and coverage from k-mer content.
    * Module written by @vezzi
* [**MACS2**](https://github.com/taoliu/MACS)
    * Summary of redundant rate from MACS2 peak calling.
* [**QoRTs**](http://hartleys.github.io/QoRTs/)
    * QoRTs is toolkit for analysis, QC and data management of RNA-Seq datasets.
* [**THetA2**](http://compbio.cs.brown.edu/projects/theta/)
    * THeTA2 _(Tumor Heterogeneity Analysis)_ estimates tumour purity and clonal / subclonal copy number.

#### Module updates:

* **BCFtools**
    * Option to collapse complementary changes in substitutions plot, useful for non-strand specific experiments (thanks to @vladsaveliev)
* **Bismark**
    * M-Bias plots no longer show read 2 for single-end data.
* **Custom Content**
    * New option to print raw HTML content to the report.
* **FastQ Screen**
    * Fixed edge-case bug where many-sample plot broke if total number of reads was less than the subsample number.
    * Fixed incorrect logic of config option `fastqscreen_simpleplot` (thanks to @daler)
    * Organisms now alphabetically sorted in fancy plot so that order is nonrandom (thanks to @daler)
    * Fixed bug where `%No Hits` was missed in logs from recent versions of FastQ Screen.
* **HTSeq Counts**
    * Fixed but so that module still works when `--additional-attr` is specified in v0.8 HTSeq above (thanks to @nalcala)
* **Picard**
    * CollectInsertSize: Fixed bug that could make the General Statistics _Median Insert Size_ value incorrect.
    * Fixed error in sample name regex that left trailing `]` characters and was generally broken (thanks to @jyh1 for spotting this)
* **Preseq**
    * Improved plots display (thanks to @vladsaveliev)
* **Qualimap**
    * Only calculate bases over target coverage for values in General Statistics. Should give a speed increase for very high coverage datasets.
* **QUAST**
    * Module is now compatible with runs from [MetaQUAST](http://quast.sourceforge.net/metaquast) (thanks to @vladsaveliev)
* **RSeQC**
    * Changed default order of sections
    * Added config option to reorder and hide module report sections

#### New MultiQC features:

* If a report already exists, execution is no longer halted.
    * `_1` is appended to the filename, iterating if this also exists.
    * `-f`/`--force` still overwrites existing reports as before
    * Feature written by [@Hammarn](https://github.com/Hammarn)
* New ability to run modules multiple times in a single report
    * Each run can be given different configuration options, including filters for input files
    * For example, have FastQC after trimming as well as FastQC before trimming.
    * See the relevant [documentation](http://multiqc.info/docs/#order-of-modules) for more instructions.
* New option to customise the order of report _sections_
    * This is in addition / alternative to changing the order of module execution
    * Allows one module to have sections in multiple places (eg. Custom Content)
* Tables have new column options `floor`, `ceiling` and `minRange`.
* Reports show warning if JavaScript is disabled
* Config option `custom_logo` now works with file paths relative to config file directory and cwd.

#### Bug Fixes:

* Table headers now sort columns again after scrolling the table
* Fixed buggy table header tooltips
* Base `clean_s_name` function now strips excess whitespace.
* Line graphs don't smooth lines if not needed (number of points < maximum number allowed)
* PDF output now respects custom output directory.

---

## [MultiQC v1.0](https://github.com/ewels/MultiQC/releases/tag/v1.0) - 2017-05-17
Version 1.0! This release has been a long time coming and brings with it some fairly
major improvements in speed, report filesize and report performance. There's also
a bunch of new modules, more options, features and a whole lot of bug fixes.

The version number is being bumped up to 1.0 for a couple of reasons:

1. MultiQC is now _(hopefully)_ relatively stable. A number of facilities and users
   are now using it in a production setting and it's published. It feels like it
   probably deserves v1 status now somehow.
2. This update brings some fairly major changes which will break backwards
   compatibility for plugins. As such, semantic versioning suggests a change in
   major version number.

### Breaking Changes
For most people, you shouldn't have any problems upgrading. There are two
scenarios where you may need to make changes with this update:

#### 1. You have custom file search patterns
Search patterns have been flattened and may no longer have arbitrary depth.
For example, you may need to change the following:
```yaml
fastqc:
    data:
        fn: 'fastqc_data.txt'
    zip:
        fn: '*_fastqc.zip'
```
to this:
```yaml
fastqc/data:
    fn: 'fastqc_data.txt'
fastqc/zip:
    fn: '*_fastqc.zip'
```
See the [documentation](http://multiqc.info/docs/#step-1-find-log-files) for instructions on how to write the new file search syntax.

See [`search_patterns.yaml`](multiqc/utils/search_patterns.yaml) for the new module search keys
and more examples.

####  2. You have custom plugins / modules / external code
To see what changes need to applied to your custom plugin code, please see the [MultiQC docs](http://multiqc.info/docs/#v1.0-updates).

#### New Modules:

* [**Adapter Removal**](https://github.com/mikkelschubert/adapterremoval)
    * AdapterRemoval v2 - rapid adapter trimming, identification, and read merging
* [**BUSCO**](http://busco.ezlab.org/)
    * New module for the `BUSCO v2` tool, used for assessing genome assembly and annotation completeness.
* [**Cluster Flow**](http://clusterflow.io)
    * Cluster Flow is a workflow tool for bioinformatics pipelines. The new module parses executed tool commands.
* [**RNA-SeQC**](http://archive.broadinstitute.org/cancer/cga/rna-seqc)
    * New module to parse output from RNA-SeQC, a java program which computes a series
    of quality control metrics for RNA-seq data.
* [**goleft indexcov**](https://github.com/brentp/goleft/tree/master/indexcov)
    * [goleft indexcov](https://github.com/brentp/goleft/tree/master/indexcov) uses the PED and ROC
    data files to create diagnostic plots of coverage per sample, helping to identify sample gender and coverage issues.
    * Thanks to @chapmanb and @brentp
* [**SortMeRNA**](http://bioinfo.lifl.fr/RNA/sortmerna/)
    * New module for `SortMeRNA`, commonly used for removing rRNA contamination from datasets.
    * Written by @bschiffthaler

#### Module updates:

* **Bcftools**
    * Fixed bug with display of indels when only one sample
* **Cutadapt**
    * Now takes the filename if the sample name is `-` (stdin). Thanks to @tdido
* **FastQC**
    * Data for the Sequence content plot can now be downloaded from reports as a JSON file.
* **FastQ Screen**
    * Rewritten plotting method for high sample numbers plot (~ > 20 samples)
    * Now shows counts for single-species hits and bins all multi-species hits
    * Allows plot to show proper percentage view for each sample, much easier to interpret.
* **HTSeq**
    * Fix bug where header lines caused module to crash
* **Picard**
    * New `RrbsSummaryMetrics` Submodule!
    * New `WgsMetrics` Submodule!
    * `CollectGcBiasMetrics` module now prints summary statistics to `multiqc_data` if found. Thanks to @ahvigil
* **Preseq**
    * Now trims the x axis to the point that meets 90% of `min(unique molecules)`.
  	Hopefully prevents ridiculous x axes without sacrificing too much useful information.
    * Allows to show estimated depth of coverage instead of less informative molecule counts
  	(see [details](http://multiqc.info/docs/#preseq)).
    * Plots dots with externally calculated real read counts (see [details](http://multiqc.info/docs/#preseq)).
* **Qualimap**
    * RNASeq Transcript Profile now has correct axis units. Thanks to @roryk
    * BamQC module now doesn't crash if reports don't have genome gc distributions
* **RSeQC**
    * Fixed Python3 error in Junction Saturation code
    * Fixed JS error for Junction Saturation that made the single-sample combined plot only show _All Junctions_

#### Core MultiQC updates:
* Change in module structure and import statements (see [details](http://multiqc.info/docs/#v1.0-updates)).
* Module file search has been rewritten (see above changes to configs)
    * Significant improvement in search speed (test dataset runs in approximately half the time)
    * More options for modules to find their logs, eg. filename and contents matching regexes (see the [docs](http://multiqc.info/docs/#step-1-find-log-files))
* Report plot data is now compressed, significantly reducing report filesizes.
* New `--ignore-samples` option to skip samples based on parsed sample name
    * Alternative to filtering by input filename, which doesn't always work
    * Also can use config vars `sample_names_ignore` (glob patterns) and `sample_names_ignore_re` (regex patterns).
* New `--sample-names` command line option to give file with alternative sample names
    * Allows one-click batch renaming in reports
* New `--cl_config` option to supply MultiQC config YAML directly on the command line.
* New config option to change numeric multiplier in General Stats
    * For example, if reports have few reads, can show `Thousands of Reads` instead of `Millions of Reads`
    * Set config options `read_count_multiplier`, `read_count_prefix` and `read_count_desc`
* Config options `decimalPoint_format` and `thousandsSep_format` now apply to tables as well as plots
    * By default, thosands will now be separated with a space and `.` used for decimal places.
* Tables now have a maximum-height by default and scroll within this.
    * Speeds up report rendering in the web browser and makes report less stupidly long with lots of samples
    * Button beneath table toggles full length if you want a zoomed-out view
    * Refactored and removed previous code to make the table header "float"
    * Set `config.collapse_tables` to `False` to disable table maximum-heights
* Bar graphs and heatmaps can now be zoomed in on
    * Interactive plots sometimes hide labels due to lack of space. These can now be zoomed in on to see specific samples in more detail.
* Report plots now load sequentially instead of all at once
    * Prevents the browser from locking up when large reports load
* Report plot and section HTML IDs are now sanitised and checked for duplicates
* New template available (called _sections_) which has faster loading
    * Only shows results from one module at a time
    * Makes big reports load in the browser much more quickly, but requires more clicking
    * Try it out by specifying `-t sections`
* Module sections tidied and refactored
    * New helper function `self.add_section()`
    * Sections hidden in nav if no title (no more need for the hacky `self.intro += `)
    * Content broken into `description`, `help` and `plot`, with automatic formatting
    * Empty module sections are now skipped in reports. No need to check if a plot function returns `None`!
    * Changes should be backwards-compatible
* Report plot data export code refactored
    * Now doesn't export hidden samples (uses HighCharts [export-csv](https://github.com/highcharts/export-csv) plugin)
* Handle error when `git` isn't installed on the system.
* Refactored colouring of table cells
    * Was previously done in the browser using [chroma.js](http://gka.github.io/chroma.js/)
    * Now done at report generation time using the [spectra](https://pypi.python.org/pypi/spectra) package
    * Should helpfully speed up report rendering time in the web browser, especially for large reports
* Docs updates (thanks to @varemo)
* Previously hidden log file `.multiqc.log` renamed to `multiqc.log` in `multiqc_data`
* Added option to load MultiQC config file from a path specified in the environment variable `MULTIQC_CONFIG_PATH`
* New table configuration options
    * `sortRows: False` prevents table rows from being sorted alphabetically
    * `col1_header` allows the default first column header to be changed from "Sample Name"
* Tables no longer show _Configure Columns_ and _Plot_ buttons if they only have a single column
* Custom content updates
    * New `custom_content`/`order` config option to specify order of Custom Content sections
    * Tables now use the header for the first column instead of always having `Sample Name`
    * JSON + YAML tables now remember order of table columns
    * Many minor bugfixes
* Line graphs and scatter graphs axis limits
    * If limits are specified, data exceeding this is no longer saved in report
    * Visually identical, but can make report file sizes considerable smaller in some cases
* Creating multiple plots without a config dict now works (previously just gave grey boxes in report)
* All changes are now tested on a Windows system, using [AppVeyor](https://ci.appveyor.com/project/ewels/multiqc/)
* Fixed rare error where some reports could get empty General Statistics tables when no data present.
* Fixed minor bug where config option `force: true` didn't work. Now you don't have to always specify `-f`!


---

## [MultiQC v0.9](https://github.com/ewels/MultiQC/releases/tag/v0.9) - 2016-12-21
A major new feature is released in v0.9 - support for _custom content_. This means
that MultiQC can now easily include output from custom scripts within reports without
the need for a new module or plugin. For more information, please see the
[MultiQC documentation](http://multiqc.info/docs/#custom-content).

#### New Modules:

* [**HTSeq**](http://www-huber.embl.de/HTSeq/doc/count.html)
    * New module for the `htseq-count` tool, often used in RNA-seq analysis.
* [**Prokka**](http://www.vicbioinformatics.com/software.prokka.shtml)
    * Prokka is a software tool for the rapid annotation of prokaryotic genomes.
* [**Slamdunk**](http://t-neumann.github.io/slamdunk/)
    * Slamdunk is a software tool to analyze SLAMSeq data.
* [**Peddy**](https://github.com/brentp/peddy)
    * Peddy calculates genotype :: pedigree correspondence checks, ancestry checks and sex checks using VCF files.

#### Module updates:

* **Cutadapt**
    * Fixed bug in General Stats table number for old versions of cutadapt (pre v1.7)
    * Added support for _really_ old cutadapt logs (eg. v.1.2)
* **FastQC**
    * New plot showing total overrepresented sequence percentages.
    * New option to parse a file containing a theoretical GC curve to display in the background.
        * Human & Mouse Genome / Transcriptome curves bundled, or make your own using
          [fastqcTheoreticalGC](https://github.com/mikelove/fastqcTheoreticalGC). See the
          [MultiQC docs](http://multiqc.info/docs/#fastqc) for more information.
* **featureCounts**
    * Added parsing checks and catch failures for when non-featureCounts files are picked up by accident
* **GATK**
    * Fixed logger error in VariantEval module.
* **Picard**
    * Fixed missing sample overwriting bug in `RnaSeqMetrics`
    * New feature to customise coverage shown from `HsMetrics` in General Statistics table
    see the [docs](http://multiqc.info/docs/#picard) for info).
    * Fixed compatibility problem with output from `CollectMultipleMetrics` for `CollectAlignmentSummaryMetrics`
* **Preseq**
    * Module now recognises output from `c_curve` mode.
* **RSeQC**
    * Made the gene body coverage plot show the percentage view by default
    * Made gene body coverage properly handle sample names
* **Samtools**
    * New module to show duplicate stats from `rmdup` logs
    * Fixed a couple of niggles in the idxstats plot
* **SnpEff**
    * Fixed swapped axis labels in the Variant Quality plot
* **STAR**
    * Fixed crash when there are 0 unmapped reads.
    * Sample name now taken from the directory name if no file prefix found.
* **Qualimap BamQC**
    * Add a line for pre-calculated reference genome GC content
    * Plot cumulative coverage for values above 50x, align with the coverage histogram.
    * New ability to customise coverage thresholds shown in General Statistics table
    (see the [docs](http://multiqc.info/docs/#qualimap) for info).

#### Core MultiQC updates:
* Support for _custom content_ (see top of release notes).
* New ninja report tool: make scatter plots of any two table columns!
* Plot data now saved in `multiqc_data` when 'flat' image plots are created
    * Allows you easily re-plot the data (eg. in Excel) for further downstream investigation
* Added _'Apply'_ button to Highlight / Rename / Hide.
    * These tools can become slow with large reports. This means that you can enter several
    things without having to wait for the report to replot each change.
* Report heatmaps can now be sorted by highlight
* New config options `decimalPoint_format` and `thousandsSep_format`
    * Allows you to change the default `1 234.56` number formatting for plots.
* New config option `top_modules` allows you to specify modules that should come at the top of the report
* Fixed bar plot bug where missing categories could shift data between samples
* Report title now printed in the side navigation
* Missing plot IDs added for easier plot exporting
* Stopped giving warnings about skipping directories (now a debug message)
* Added warnings in report about missing functionality for flat plots (exporting and toolbox)
* Export button has contextual text for images / data
* Fixed a bug where user config files were loaded twice
* Fixed bug where module order was random if `--module` or `--exclude` was used.
* Refactored code so that the order of modules can be changed in the user config
* Beefed up code + docs in scatter plots back end and multiple bar plots.
* Fixed a few back end nasties for Tables
    * Shared-key columns are no longer forced to share colour schemes
    * Fixed bug in lambda modified values when format string breaks
    * Supplying just data with no header information now works as advertised
* Improvements to back end code for bar plots
    * New `tt_decimals` and `tt_suffix` options for bar plots
    * Bar plots now support `yCeiling`, `yFloor` and `yMinRange`, as with line plots.
    * New option `hide_zero_cats:False` to force legends to be shown even when all data is 0
* General Stats _Showing x of y_ columns count is fixed on page load.
* Big code whitespace cleanup

---

## [MultiQC v0.8](https://github.com/ewels/MultiQC/releases/tag/v0.8) - 2016-09-26

#### New Modules:

* [**GATK**](https://software.broadinstitute.org/gatk/)
    * Added support for VariantEval reports, only parsing a little of the information
    in there so far, but it's a start.
    * Module originally written by @robinandeer at the [OBF Codefest](https://www.open-bio.org/wiki/Codefest_2016),
    finished off by @ewels
* [**Bcftools**](https://samtools.github.io/bcftools/)
* [**QUAST**](http://quast.bioinf.spbau.ru/)
    * QUAST is a tool for assessing de novo assemblies against reference genomes.

#### Module updates:

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

#### Core MultiQC updates:
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
* Fixed some buggy behaviour in saving / loading report highlighting + renaming configs (#354)

Many thanks to those at the [OpenBio Codefest 2016](https://www.open-bio.org/wiki/Codefest_2016)
who worked on MultiQC projects.

---

## [MultiQC v0.7](https://github.com/ewels/MultiQC/releases/tag/v0.7) - 2016-07-04
#### Module updates:
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

#### Core MultiQC updates:
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


---

## [MultiQC v0.6](https://github.com/ewels/MultiQC/releases/tag/v0.6) - 2016-04-29
#### Module updates:
* New [Salmon](http://combine-lab.github.io/salmon/) module.
* New [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) module.
* New [Bamtools stats](https://github.com/pezmaster31/bamtools) module.
* New beeswarm plot type. General Stats table replaced with this when many samples in report.
* New RSeQC module: Actually a suite of 8 new modules supporting various outputs from RSeQC
* Rewrote bowtie2 module: Now better at parsing logs and tries to scrape input from wrapper logs.
* Made cutadapt show counts by default instead of obs/exp
* Added percentage view to Picard insert size plot

#### Core MultiQC updates:
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

---

## [MultiQC v0.5](https://github.com/ewels/MultiQC/releases/tag/v0.5) - 2016-03-29
#### Module updates:
* New [Skewer](https://github.com/relipmoc/skewer) module, written by @dakl
* New [Samblaster](https://github.com/GregoryFaust/samblaster) module, written by @dakl
* New [Samtools stats](http://www.htslib.org/) module, written by @lpantano
* New [HiCUP](http://www.bioinformatics.babraham.ac.uk/projects/hicup/) module
* New [SnpEff](http://snpeff.sourceforge.net/) module
* New [methylQA](http://methylqa.sourceforge.net/) module

#### Core MultiQC updates:
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

---

## [MultiQC v0.4](https://github.com/ewels/MultiQC/releases/tag/v0.4) - 2016-02-16
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

---

## [MultiQC v0.3.2](https://github.com/ewels/MultiQC/releases/tag/v0.3.2) - 2016-02-08
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


---

## [MultiQC v0.3.1](https://github.com/ewels/MultiQC/releases/tag/v0.3.1) - 2015-11-04
* Hotfix patch to fix broken FastQC module (wasn't finding `.zip` files properly)
* General Stats table colours now flat. Should improve browser speed.
* Empty rows now hidden if appear due to column removal in general stats
* FastQC Kmer plot removed until we have something better to show.

---

## [MultiQC v0.3](https://github.com/ewels/MultiQC/releases/tag/v0.3) - 2015-11-04
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

---

## [MultiQC v0.2](https://github.com/ewels/MultiQC/releases/tag/v0.2) - 2015-09-18
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

---

## [MultiQC v0.1](https://github.com/ewels/MultiQC/releases/tag/v0.1) - 2015-09-01
* The first public release of MultiQC, after a month of development. Basic
structure in place and modules for FastQC, FastQ Screen, Cutadapt, Bismark,
STAR, Bowtie, Subread featureCounts and Picard MarkDuplicates. Approaching
stability, though still under fairly heavy development.
