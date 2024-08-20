# MultiQC Version History

## [MultiQC v1.24.1](https://github.com/MultiQC/MultiQC/releases/tag/v1.24.1) - 2024-08-21

A bug fix release mainly to restore compatibility with Python 3.8. Aside from that, few other minor bug fixes:

- FastQC: fix long-standing issue misplacing status labels when `anchor` is specified in the custom config ([#2790](https://github.com/MultiQC/MultiQC/pull/2790))
- Freyja: handle empty inputs, and ensure deterministic sample order ([#2788](https://github.com/MultiQC/MultiQC/pull/2788))
- Allow numeric xcats and ycats for the heatmap plot ([#2787](https://github.com/MultiQC/MultiQC/pull/2787))
- Make sure that config's `extra_fn_clean_exts` and `fn_clean_exts` don't conflict when both specified ([#2783](https://github.com/MultiQC/MultiQC/pull/2783))

## [MultiQC v1.24](https://github.com/MultiQC/MultiQC/releases/tag/v1.24) - 2024-08-19

Mostly a maintenance release, containing several bugfixes, performance improvements,
plus 6 new modules, along with improvements of the existing modules.

The most significant performance boost got the Kraken and Mosdepth modules, that now don't
take way more memory and CPU than any other typical module:

| Tool     | Data Set       | Memory - Before | CPU - Before | Memory - After | CPU - After |
| -------- | -------------- | --------------- | ------------ | -------------- | ----------- |
| mosdepth | 1 set of files | 196 Mb          | 3.04s        | 129 Mb         | 2.27s       |
|          | 10             | 464 Mb          | 8.48s        | 131 Mb         | 6.11s       |
|          | 100            | 3,719 Mb        | 63.19s       | 172 Mb         | 43.12s      |
| kraken   | 1 set of files | 155 Mb          | 2.07s        | 132 Mb         | 2.20s       |
|          | 10             | 606 Mb          | 8.39s        | 180 Mb         | 3.47s       |
|          | 100            | 4,970 Mb        | 71.89s       | 809 Mb         | 14.53s      |

Large plots that may hang browser are now not loaded by default, and the user can click
a button to load, so the heavy plots don't slow down the initial report rendering. This
is controlled by the `config.plots_defer_loading_numseries: 100` option.

### Updates

- Search patterns: allow multiple values for `contents` ([#2696](https://github.com/MultiQC/MultiQC/pull/2696))
- Custom content:
  - Allow dict input to heatmap ([#2761](https://github.com/MultiQC/MultiQC/pull/2761))
  - Allow multiple custom content to general stats table ([#2727](https://github.com/MultiQC/MultiQC/pull/2727))
- Plots:
  - Defer render of plots if number of samples > `config.plots_defer_loading_numseries` ([#2759](https://github.com/MultiQC/MultiQC/pull/2759), [#2777](https://github.com/MultiQC/MultiQC/pull/2777), [#2773](https://github.com/MultiQC/MultiQC/pull/2773), [#2774](https://github.com/MultiQC/MultiQC/pull/2774))
  - Line plot: show markers when num of data points < `config.lineplot_number_of_points_to_hide_markers` (=50) ([#2760](https://github.com/MultiQC/MultiQC/pull/2760)). As a nice consequence, trivial lines of a single data point become visible.
  - Line plot: smooth by default to 500 points on the X axis to avoid inflating the report file size ([#2776](https://github.com/MultiQC/MultiQC/pull/2776))
  - Allow to configure the scale of the exported plot fonts through the config option `config.plots_export_font_scale` ([#2758](https://github.com/MultiQC/MultiQC/pull/2758))
  - Improve the performance of loading large tables in browser ([#2737](https://github.com/MultiQC/MultiQC/pull/2737))
  - Fix the toolbox highlight of the line plots ([#2724](https://github.com/MultiQC/MultiQC/pull/2724))
- The function that returns built plots in an interactive session now uses the module anchor (or lowercase module name) to key the results ([#2741](https://github.com/MultiQC/MultiQC/pull/2741))
- More helpful config validation error: print the parent model name, if applicable ([#2709](https://github.com/MultiQC/MultiQC/pull/2709))

### New modules

- [VG](https://github.com/vgteam/vg) ([#2690](https://github.com/MultiQC/MultiQC/pull/2690)), a toolkit to manipulate graphical genomes. The module parses [vg-stats](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe#evaluating-with-vg-stats) reports that summarize alignment stats from a GAM file.
- [ngs-bits](<(https://github.com/imgag/ngs-bits)>) ([#2231](https://github.com/MultiQC/MultiQC/pull/2231)). A tool that calculating statistics from FASTQ, BAM, and VCF files. The module parses XML output generated for two tools in the ngs-bits collection:
  - [ReadQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/ReadQC.md) for FastQ file stats
  - [MappingQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC.md) for BAM file stats
- [Pairtools](https://github.com/mirnylab/pairtools) ([#1148](https://github.com/MultiQC/MultiQC/pull/1148)). A toolkit for Chromatin Conformation Capture experiments. Handles short-reads paired reference alignments, extracts 3C-specific information, and perform common tasks such as sorting, filtering, and deduplication. The module parses summary statistics generated by pairtools's `dedup` and `stats` tools.
- [nanoq](https://github.com/nerdna/nanoq/) ([#2723](https://github.com/MultiQC/MultiQC/pull/2723)). A tool that reports read quality and length from nanopore sequencing data.
- [Ganon](https://pirovc.github.io/ganon/) ([#1935](https://github.com/MultiQC/MultiQC/pull/1935)). A tool for metagenomics classification: quickly assigns sequence fragments to their closest reference among thousands of references via Interleaved Bloom Filters of k-mer/minimizers

### Fixes

- Fix `--pdf` option to generate `multiqc_report.pdf` ([#2733](https://github.com/MultiQC/MultiQC/pull/2733))
- Fix saving table plots to file ([#2735](https://github.com/MultiQC/MultiQC/pull/2735))
- Fix adding software versions when `config.run_modules` is set ([#2755](https://github.com/MultiQC/MultiQC/pull/2755))
- Fix toolbox highlight in line plots ([#2724](https://github.com/MultiQC/MultiQC/pull/2724))
- Refactor `write_results` to avoid dynamically overriding `config`, fixes `module.write_data_file` ([#2722](https://github.com/MultiQC/MultiQC/pull/2722))
- Search stats: do not double-count ignored files ([#2708](https://github.com/MultiQC/MultiQC/pull/2708))
- Escape values passed to HTML properties (e.g. `val` in tables) ([#2706](https://github.com/MultiQC/MultiQC/pull/2706))
- Fix re-loading explicit user configs in interactive sessions ([#2704](https://github.com/MultiQC/MultiQC/pull/2704))
- Fix file search performance regression ([#2762](https://github.com/MultiQC/MultiQC/pull/2762))
- Fix handling module href string ([#2739](https://github.com/MultiQC/MultiQC/pull/2739))
- Custom content:
  - Cast names to strings instead of asserting (allows numerical sample names) ([#2769](https://github.com/MultiQC/MultiQC/pull/2769))
  - Fix parsing custom content tables with numerical samples ([#2705](https://github.com/MultiQC/MultiQC/pull/2705))
  - Section order and custom content order: allow skip the `*-module` suffix from anchors ([#2770](https://github.com/MultiQC/MultiQC/pull/2770))

### Module updates

- Kraken: optimize memory and runtime ([#2756](https://github.com/MultiQC/MultiQC/pull/2756))
- Mosdepth: optimize memory and runtime ([#2748](https://github.com/MultiQC/MultiQC/pull/2748), [#2749](https://github.com/MultiQC/MultiQC/pull/2749))
- Abstract `config.get_cov_thresholds` function for `mosdepth` and `qualimap` ([#2707](https://github.com/MultiQC/MultiQC/pull/2707))
- Anglerfish: adjust for version 0.6.1 ([#2757](https://github.com/MultiQC/MultiQC/pull/2757))
- Umitools: prefer output for sample name, handle the `<stdin>`/`<stdout>` placeholders ([#2698](https://github.com/MultiQC/MultiQC/pull/2698))
- Kraken: fix top % calculation, more efficient total read count calculation ([#2744](https://github.com/MultiQC/MultiQC/pull/2744))
- Bracken: when printing number of found samples, indicate that running Bracken not Kraken ([#2743](https://github.com/MultiQC/MultiQC/pull/2743))
- Nonpareil: Update docs about new version that generates the JSON file ([#2734](https://github.com/MultiQC/MultiQC/pull/2734))
- STAR: add all alignment summary metrics into a new separate table ([#1828](https://github.com/MultiQC/MultiQC/pull/1828))
- Pairtools: fix typos and grammar, remove redundancies ([#2711](https://github.com/MultiQC/MultiQC/pull/2711))
- Peddy sex plot: color predicted sex ([#2778](https://github.com/MultiQC/MultiQC/pull/2778))
- FastQC: for plots with bp on the X axis, use the interval start instead of average ([#2790](https://github.com/MultiQC/MultiQC/pull/2790))

### Module fixes

- Cellranger: fix for missing `analysis_tab` data ([#2771](https://github.com/MultiQC/MultiQC/pull/2771))
- Fix setting coverage thresholds for `mosdepth` and `qualimap` ([#2754](https://github.com/MultiQC/MultiQC/pull/2754))
- Nonpareil: fix running with >12 samples ([#2752](https://github.com/MultiQC/MultiQC/pull/2752))
- Bracken: fix bug when direct reads not classified ([#2738](https://github.com/MultiQC/MultiQC/pull/2738))
- ngsderive: fix ValueError in the `encoding` submodule ([#2740](https://github.com/MultiQC/MultiQC/pull/2740))
- RSeQC: fix duplicated namespace ([#2732](https://github.com/MultiQC/MultiQC/pull/2732))
- Glimpse: fix parsing data, add proper type hints ([#2721](https://github.com/MultiQC/MultiQC/pull/2721))
- Fix ignoring samples in `spaceranger`, `ngsbits`, `isoseq`, `dragen coverage` ([#2717](https://github.com/MultiQC/MultiQC/pull/2717))
- Glimpse: clean and fix filtering samples ([#2716](https://github.com/MultiQC/MultiQC/pull/2716))
- Spaceranger: fix ignoring samples ([#2714](https://github.com/MultiQC/MultiQC/pull/2714))

### Refactoring

- Refactor `write_results` to avoid dynamically overriding `config`, fixes `module.write_data_file` ([#2722](https://github.com/MultiQC/MultiQC/pull/2722))
- Modules:
  - isoseq and odgi: fix module warnings and error handling ([#2718](https://github.com/MultiQC/MultiQC/pull/2718))
  - Qualimap: refactor and add type hints ([#2707](https://github.com/MultiQC/MultiQC/pull/2707))
  - FastQC: refactor and add type hints ([#2763](https://github.com/MultiQC/MultiQC/pull/2763))
  - Kraken: refactor and add type hints ([#2744](https://github.com/MultiQC/MultiQC/pull/2744))
  - Cell Ranger: refactor, add type hints, get rid of module mixins and fields ([#2775](https://github.com/MultiQC/MultiQC/pull/2775))
  - Spaceranger: refactor, add type hints, get rid of module mixins and fields ([#2714](https://github.com/MultiQC/MultiQC/pull/2714))
  - RSeQC: refactor, add type hints, and remove `multiqc_rseqc.js` ([#2710](https://github.com/MultiQC/MultiQC/pull/2710))

### Infrastructure

- Split up the tests for sample versions discovery ([#2751](https://github.com/MultiQC/MultiQC/pull/2751))
- Move integration test variations into unit tests, actually test them ([#2713](https://github.com/MultiQC/MultiQC/pull/2713))
- Move module docs to the docstrings & generate `docs/modules/*.md` from the docstrings using a separate script ([#2703](https://github.com/MultiQC/MultiQC/pull/2703))
- Embed the search patterns into the module docs ([#2765](https://github.com/MultiQC/MultiQC/pull/2765))
- Dockerfile: use `COPY` instead of `ADD` to copy only relevant files, update base image to Python 3.12 ([#2700](https://github.com/MultiQC/MultiQC/pull/2700))

## [MultiQC v1.23](https://github.com/MultiQC/MultiQC/releases/tag/v1.23) - 2024-07-09

Bug fixes, integration of `pytest` and `mypy`, and one new module.

From the user perspective, this is mostly a maintenance release, containing several
important bugfixes, plus minor improvements and a new module - Glimpse.

For developers, there are two significant additions to the CI workflow:

- [pytest](https://docs.pytest.org/), along with unit tests covering the core library,
- and [mypy](https://mypy-lang.org/), along with ensuring that the core codebase is fully
  type-annotated.

The core unit tests are located in the `multiqc/tests` folder, and the module tests are
located in the corresponding `multiqc/modules/*/tests` subfolders.
The [CI workflows](https://github.com/MultiQC/MultiQC/tree/main/.github/workflows)
are refactored to separate the integration tests and the unit tests, to improve the
granularity and parallelization. The tests are discovered and executed with pytest,
and the coverage is reported by [codecov](https://app.codecov.io/gh/MultiQC/MultiQC).

The `multiqc/tests` subfolder has several test files that cover most of the core library.
It also has a [test_modules_run.py](https://github.com/MultiQC/MultiQC/blob/main/tests/test_modules_run.py)
tests that checks that every module didn't crash when being run on the corresponding data
in [test-data](https://github.com/MultiQC/test-data), and added _something_ into the report.
That is somewhat of a blanket test for modules, that doesn't check if the modules logic
worked correctly. For that reason, the users are encouraged to write more comprehensive
tests that take the specific module logic into account, and place them in
`multiqc/modules/*/tests`. For some initial examples, consider checking:

- The [samtools flagstat](https://github.com/MultiQC/MultiQC/blob/main/multiqc/modules/samtools/tests/test_flagstat.py)
  test that verifies some logic in the `flagstat` submodule of the `samtools` module;
- The [picard tools](https://github.com/MultiQC/MultiQC/blob/main/multiqc/modules/picard/tests/test_picard.py)
  test that checks that every submodule for each Picard tool worked correctly.

### Other changes & fixes

#### Fixes

- Custom content"
  - Multiple fixes of the custom content parsing logic ([#2674](https://github.com/MultiQC/MultiQC/pull/2674))
  - Fix parsing custom content submodules with a custom config ([#2654](https://github.com/MultiQC/MultiQC/pull/2654))
  - Line plot: fix spreading `extra_series` between multiple datasets ([#2684](https://github.com/MultiQC/MultiQC/pull/2684))
  - Line plot series: allow pass lists instead of x/y tuples to work properly with YAML ([#2683](https://github.com/MultiQC/MultiQC/pull/2683))
- Re-enabling the `software_version` module section ([#2670](https://github.com/MultiQC/MultiQC/pull/2670))
- When `--no-ansi` is set, disable colors in `rich_click` too ([#2678](https://github.com/MultiQC/MultiQC/pull/2678))
- Support CWD path filters ( `./path/...`) in config ([#2676](https://github.com/MultiQC/MultiQC/pull/2676))
- Fix writing report to stdout with `--filename stdout`, log to stderr ([#2672](https://github.com/MultiQC/MultiQC/pull/2672))
- Interactive use:
  - Reset all config values in `config.reset()`, even those that are not in `config_default.yaml` ([#2660](https://github.com/MultiQC/MultiQC/pull/2660))
  - Reset config between `multiqc.run()` calls ([#2655](https://github.com/MultiQC/MultiQC/pull/2655))
  - Better handling of calling `write_report` twice ([#2688](https://github.com/MultiQC/MultiQC/pull/2688))

#### Updates

- Run mypy on core library ([#2665](https://github.com/MultiQC/MultiQC/pull/2665))
- Add tests for plot export ([#2682](https://github.com/MultiQC/MultiQC/pull/2682))
- Add tests for command line use, including for passing `TMPDIR` ([#2677](https://github.com/MultiQC/MultiQC/pull/2677))
- Custom content: allow hash-fenced table columns ([#2649](https://github.com/MultiQC/MultiQC/pull/2649))
- Software versions: parse for sorting, but preserve the original strings ([#2671](https://github.com/MultiQC/MultiQC/pull/2671))
- Allow both table-level and column-level custom plot config for table ([#2662](https://github.com/MultiQC/MultiQC/pull/2662))

#### New modules

- Glimpse ([#2492](https://github.com/MultiQC/MultiQC/pull/2492))

#### Module fixes

- Fix parsing kraken vs. bracken: respect `num_lines` in search patterns ([#2657](https://github.com/MultiQC/MultiQC/pull/2657))
- Fix the `bbmap/qchist` search pattern ([#2661](https://github.com/MultiQC/MultiQC/pull/2661))

#### Module updates

- Picard HsMetrics: support any custom X coverage metrics ([#2663](https://github.com/MultiQC/MultiQC/pull/2663))
- Samtools coverage: avoid hard crash for invalid file contents ([#2664](https://github.com/MultiQC/MultiQC/pull/2664))

#### Refactoring

- Abstract code related to temporary directory creation into a separate module ([#2675](https://github.com/MultiQC/MultiQC/pull/2675))

#### Infrastructure

- Use pull-request labels and milestones for changelog generation ([#2691](https://github.com/MultiQC/MultiQC/pull/2691))

## [MultiQC v1.22.3](https://github.com/MultiQC/MultiQC/releases/tag/v1.22.3) - 2024-06-22

Contains fixes of multiple bugs collected after the last release, along with few minor improvements.

### MultiQC fixes

- Fix the `re_contents` search patterns when pattern is found in the middle of the file. Fixes finding logs from several Picard submodules, like `CollectRnaSeqMetrics` and `CollectWgsMetrics` in some cases ([#2610](https://github.com/MultiQC/MultiQC/pull/2610))
- Fixes the `run_modules` option use when the module anchor doesn't match the module entry point ID (e.g. `DRAGEN` and `dragen`) ([#2633](https://github.com/MultiQC/MultiQC/pull/2633))
- Fix use of custom search patterns for custom content ([#2647](https://github.com/MultiQC/MultiQC/pull/2647))
- Fix plot export with `export_plots: true` or `--export` ([#2637](https://github.com/MultiQC/MultiQC/pull/2637))
- Correctly handle old-style `label` sections in `x_lines` or `y_lines` in line plot configs ([#2648](https://github.com/MultiQC/MultiQC/pull/2648))
- Fix disabling `sort_rows` in custom content by subclassing `TableConfig` from `ValidatedConfig` and use deprecated ([#2604](https://github.com/MultiQC/MultiQC/pull/2604))
- When user provides a search pattern dictionary in config, recursively update instead of replacing ([#2620](https://github.com/MultiQC/MultiQC/pull/2620))
- Fix config update when dict replaced with list, e.g. a `search_patterns` item is a list that's replaced with a dict (https://github.com/MultiQC/MultiQC/commit/c388178bb6d9f143c6f8e8b0146647a067021ea4)

### MultiQC updates

- Add unit tests for core and some modules (see `picard` or `samtools`), as well as `codecov` report ([#2624](https://github.com/MultiQC/MultiQC/pull/2624))
  - Now MultiQC checks if every module does something productive with the provided test data in `test-data`.
  - For modules with many submodules (picard, dragen), additionally check if every submodule parses the expected number of samples from `test-data` files.
  - Users can put module tests in `test` subfolders, e.g. https://github.com/MultiQC/MultiQC/tree/main/multiqc/modules/picard/tests
  - Use `pytest` for all core unit tests ([#2623](https://github.com/MultiQC/MultiQC/pull/2623))
  - Move unit tests from the `test-data` repo into `tests` folder ([#2622](https://github.com/MultiQC/MultiQC/pull/2622))
- Plot config validation:
  - Validate line plot `x_lines`, `x_bands`, etc. with a Pydantic model, including `label` subsections ([#2648](https://github.com/MultiQC/MultiQC/pull/2648))
  - Validate line plot series and `extra_series` with a Pydantic model ([#2573](https://github.com/MultiQC/MultiQC/pull/2573))
  - Validate table config ([#2604](https://github.com/MultiQC/MultiQC/pull/2604))
  - Make the "unrecognised field" error a warning
  - Rename deprecated plot config fields in internal modules ([#2636](https://github.com/MultiQC/MultiQC/pull/2636))
- Show progress bar for exporting flat plot images ([#2639](https://github.com/MultiQC/MultiQC/pull/2639))
- Better error message for incorrect `run_modules` ([#2635](https://github.com/MultiQC/MultiQC/pull/2635))
- Increase flat plots sample number threshold to 1000 ([#2615](https://github.com/MultiQC/MultiQC/pull/2615))
- Small speed-up of the line block iterator ([#2588](https://github.com/MultiQC/MultiQC/pull/2588))
- Update README logos for better compatibility ([#2603](https://github.com/MultiQC/MultiQC/pull/2603))
- Docs: don't use raw markdown links ([#2642](https://github.com/MultiQC/MultiQC/pull/2642))
- Allow to override `showlegend` for line config plots. Default to not-show for large datasets to avoid bloated legends ([#2615](https://github.com/MultiQC/MultiQC/pull/2615))
- Show error message if failed to parse custom content header (https://github.com/MultiQC/MultiQC/commit/d736846a0c23410243f80b2bdca984363211ffc3)
- Load every found config file once https://github.com/MultiQC/MultiQC/commit/422b39bc787720cefea81a85dabbf6411b3421ac

### Module fixes and updates

- **Picard**
  - Fix finding `CollectRnaSeqMetrics` and `CollectWgsMetrics` logs by fixing the `re_contents` search patterns ([#2610](https://github.com/MultiQC/MultiQC/pull/2610))
- **biobambam2**
  - Fix parsing `markdups` logs
- **DRAGEN**
  - Coverage histograms: fix duplicated label suffix ([#2619](https://github.com/MultiQC/MultiQC/pull/2619))
  - Fix the `gc_metrics` submodule ([#2629](https://github.com/MultiQC/MultiQC/pull/2629))
  - `vc_metrics`: pre-filter numbers can be zero ([#2618](https://github.com/MultiQC/MultiQC/pull/2618))
- **FastQC**
  - Default to `showlegend: false`, as we don't distinguish the sample colors, unless `fastqc_config: status_checks: false'` is set ([#2615](https://github.com/MultiQC/MultiQC/pull/2615))
- **BBTools**
  - Fix incorrect calculation of % Q30 Bases ([#2628](https://github.com/MultiQC/MultiQC/pull/2628))
- **Samtools**
  - `markdup`: resolve inconsistent non-optical pair duplicate variable name in samtools markdup module ([#2626](https://github.com/MultiQC/MultiQC/pull/2626))
- **NanoStat**
  - Support different `Q` cutoffs ([#2645](https://github.com/MultiQC/MultiQC/pull/2645))
- **Salmon**
  - Fix ignored parsed `library_types` when its type is list ([#2617](https://github.com/MultiQC/MultiQC/pull/2617))
- **UMI-tools**
  - Improve `extract` plots ([#2614](https://github.com/MultiQC/MultiQC/pull/2614))
- **BCL Convert**
  - Fix 'pecent' typo ([#2612](https://github.com/MultiQC/MultiQC/pull/2612))

## [MultiQC v1.22.2](https://github.com/MultiQC/MultiQC/releases/tag/v1.22.2) - 2024-05-31

Bug fix release. Two major issues are fixed here:

- Fixed running same module twice with `path_filters` (e.g. trimmed vs. raw FastQC),
- The raw data `report_saved_raw_data` is re-added in multiqc_data.json by default.

### Fixes

- Fix running same module multiple times in the report (e.g. trimmed vs. raw FastQC) ([#2592](https://github.com/MultiQC/MultiQC/pull/2592))
- Preserve `report_saved_raw_data` in multiqc_data.json by keeping `preserve_module_raw_data: false` by default ([#2591](https://github.com/MultiQC/MultiQC/pull/2591))
- Table headers: do not set namespace to `None` when there is single namespace ([#2590](https://github.com/MultiQC/MultiQC/pull/2590))
- Re-enable falling back to flat plots for large datasets ([#2580](https://github.com/MultiQC/MultiQC/pull/2580))
- Reset in `multiqc.run(*)` to allow running it twice interactively ([#2598](https://github.com/MultiQC/MultiQC/pull/2598))
- Fix scatter plot in `--flat` mode when there are categorical axes ([#2600](https://github.com/MultiQC/MultiQC/pull/2600))
- Fix hiding table column with all empty values in custom content ([#2599](https://github.com/MultiQC/MultiQC/pull/2599))
- Table "Copy" button: include headers ([#2594](https://github.com/MultiQC/MultiQC/pull/2594))

### Module fixes and updates

- **QUAST**
  - Underscore attributes captured by lambdas to avoid wiping them after module finished ([#2581](https://github.com/MultiQC/MultiQC/pull/2581))
- **Cell Ranger**
  - Handle missing `vdj_annotation` and `vdj_enrichment` sections ([#2579](https://github.com/MultiQC/MultiQC/pull/2579))
- **fgbio**
  - Fix links in fgbio.md ([#2586](https://github.com/MultiQC/MultiQC/pull/2586))
- **Glimpse**:
  - Add support for Glimpse concordance metrics [#2491](https://github.com/MultiQC/MultiQC/issues/2491)
- **Custom content**
  - Support DOI for custom content ([#2582](https://github.com/MultiQC/MultiQC/pull/2582))

## [MultiQC v1.22.1](https://github.com/MultiQC/MultiQC/releases/tag/v1.22) - 2024-05-17

### Fix run as a Nextflow job

This bug fix release addresses the file search problem when MultiQC executed as a typical Nextflow job. See [#2575](https://github.com/MultiQC/MultiQC/pull/2575) for detail.

## [MultiQC v1.22](https://github.com/MultiQC/MultiQC/releases/tag/v1.22) - 2024-05-14

### Highlights - notebooks and performance

Version 1.22 brings some major behind-the-scenes refactoring to MultiQC. This unlocks a number of new features, such as the ability to use MultiQC as a Python library in scripts / notebooks, and run-time validation of plot config attributes.

This release also introduces some huge performance improvements thanks to [@rhpvorderman](https://github.com/rhpvorderman). Compared to v1.21, a typical v1.22 run is **53% faster** and has a **6x smaller peak-memory footprint** - well worth updating!

Finally, support for the depreciated HighCharts plotting library is fully removed in v1.22, bringing to a close a long standing project to migrate to Plotly.

For more information, please see the MultiQC release blog article on the Seqera website: <https://seqera.io/blog/>

### MultiQC updates

- Remove the `highcharts` template and Highcharts and Matplotlib dependencies ([#2409](https://github.com/MultiQC/MultiQC/pull/2409))
- Remove CSP.txt and the linting check, move the script that prints missing hashes under `scripts`. Admins of servers with Content Security Policy can use it to print missing hashes when they install a new MultiQC version with: `python scripts/print_missing_csp.py --report full_report.html` ([#2421](https://github.com/MultiQC/MultiQC/pull/2421))
- Do not maintain change log between releases ([#2427](https://github.com/MultiQC/MultiQC/pull/2427))
- Use native clipboard API ([#2419](https://github.com/MultiQC/MultiQC/pull/2419))
- Profile runtime: visualize per-module memory and run time ([#2548](https://github.com/MultiQC/MultiQC/pull/2548), [#2547](https://github.com/MultiQC/MultiQC/pull/2547))
- Refactoring for performance:
  - Search file blocks rather than individual lines for faster results ([#2513](https://github.com/MultiQC/MultiQC/pull/2513))
  - Refactor file content search for a 40% speed increase ([#2505](https://github.com/MultiQC/MultiQC/pull/2505))
  - Sort `filepatterns` for faster searching ([#2506](https://github.com/MultiQC/MultiQC/pull/2506))
  - Use `array.array` for in-memory plot data, stream to render Jinja and dump JSON to reduce memory requirement ([#2515](https://github.com/MultiQC/MultiQC/pull/2515))
  - Speed up all modules by caching `spectra.scale` and using sets instead of lists ([#2509](https://github.com/MultiQC/MultiQC/pull/2509))
  - Stream json data to a file to save 30% of the memory ([#2510](https://github.com/MultiQC/MultiQC/pull/2510))
  - Do `replace_nan` in place rather than creating a new object ([#2529](https://github.com/MultiQC/MultiQC/pull/2529))
  - Use gzip rather than lzstring for compression and decompression of the plot data ([#2504](https://github.com/MultiQC/MultiQC/pull/2504))
  - Use gzip level 6 for faster json compression ([#2553](https://github.com/MultiQC/MultiQC/pull/2553))
  - Clean up module raw data after running each module, significantly reduces the memory footprint ([#2551](https://github.com/MultiQC/MultiQC/pull/2551))
- Refactoring for interactivity and validation:
  - Top-level functions for MultiQC use as a library ([#2442](https://github.com/MultiQC/MultiQC/pull/2442))
  - Pydantic models for plots and datasets ([#2442](https://github.com/MultiQC/MultiQC/pull/2442))
  - Validating plot configs with Pydantic ([#2534](https://github.com/MultiQC/MultiQC/pull/2534))
  - Use dataclasses for table and violin columns ([#2546](https://github.com/MultiQC/MultiQC/pull/2546))
  - Break up the main run function into submodules ([#2446](https://github.com/MultiQC/MultiQC/pull/2446))
  - Deprecate `multiqc.utils.config` and `multiqc.utils.report` in favour of `multiqc.config` and `multiqc.report` ([#2542](https://github.com/MultiQC/MultiQC/pull/2542))
  - Static typing of the report and config modules ([#2445](https://github.com/MultiQC/MultiQC/pull/2445))
  - Add type hints into core codebase ([#2434](https://github.com/MultiQC/MultiQC/pull/2434))
  - Consistent config options: rename `decimalPlaces` to `tt_decimals` ([#2451](https://github.com/MultiQC/MultiQC/pull/2451))
  - Remove encoding and shebang headers from module files ([#2425](https://github.com/MultiQC/MultiQC/pull/2425))
  - Refactor line plot categories: keep boolean throughout the code, and data points as pairs for simplicity ([#2418](https://github.com/MultiQC/MultiQC/pull/2418))
- Fixes:
  - Fix error when using default sort ([#2544](https://github.com/MultiQC/MultiQC/pull/2544))
  - Do not attempt to render flat plot when no data ([#2490](https://github.com/MultiQC/MultiQC/pull/2490))
  - Fix export plots with `--export` and always export data ([#2489](https://github.com/MultiQC/MultiQC/pull/2489))
  - Fix: make sure `modify` lambda not present in JSON dump ([#2455](https://github.com/MultiQC/MultiQC/pull/2455))
  - Enable `--export` even when writing interactive plots ([#2444](https://github.com/MultiQC/MultiQC/pull/2444))
  - Replace `NaN` with `null` in exported JSON ([#2432](https://github.com/MultiQC/MultiQC/pull/2432))
  - Fix `y_minrange` option ([#2415](https://github.com/MultiQC/MultiQC/pull/2415))
- Reduce report size: exclude plot data for sections in `remove_sections` ([#2460](https://github.com/MultiQC/MultiQC/pull/2460))
- Add `ge` and `le` to `cond_formatting_rules` ([#2494](https://github.com/MultiQC/MultiQC/pull/2494))
- CI: use `uv pip` ([#2352](https://github.com/MultiQC/MultiQC/pull/2352))
- Lint check for use of `f["content_lines"]` ([#2485](https://github.com/MultiQC/MultiQC/pull/2485))
- Allow to set style of line graph (`lines` or `lines+markers`) per plot ([#2413](https://github.com/MultiQC/MultiQC/pull/2413))
- Add `CMD` to `Dockerfile` so a default run without any parameters displays the `--help` ([#2279](https://github.com/MultiQC/MultiQC/pull/2279))
- Custom content tables are sorted by key by default (unless `sort_rows: false` is set in config), to harmonize with tables in modules.

### New modules

- [**Hostile**](https://github.com/bede/hostile) ([#2501](https://github.com/MultiQC/MultiQC/pull/2501))
  - New module: Hostile is a short and long host reads removal tool
- [**Sequali**](https://github.com/rhpvorderman/sequali) ([#2441](https://github.com/MultiQC/MultiQC/pull/2441))
  - New module: Sequali Universal sequencing QC

### Module updates

- **Adapter Removal**
  - Standardize module names: use the came case ([#2433](https://github.com/MultiQC/MultiQC/pull/2433))
- **Bamdst**
  - Fix chromosome reports when contig data labels are missing ([#2479](https://github.com/MultiQC/MultiQC/pull/2479))
  - Fix for the case when `chromosomes.report` is not provided ([#2477](https://github.com/MultiQC/MultiQC/pull/2477))
  - Stress file name requirements for chromosomes report ([#2478](https://github.com/MultiQC/MultiQC/pull/2478))
- **BBTools**
  - Set missing values to `None` for `bbmap qahist` ([#2411](https://github.com/MultiQC/MultiQC/pull/2411))
- **Bcftools**
  - Stats: add multialleic sites column ([#2414](https://github.com/MultiQC/MultiQC/pull/2414))
- **BCL Convert**
  - Show message when no undetermined reads instead of error ([#2526](https://github.com/MultiQC/MultiQC/pull/2526))
  - Fix for absent index reads ([#2511](https://github.com/MultiQC/MultiQC/pull/2511))
  - Add all file types to sources ([#2456](https://github.com/MultiQC/MultiQC/pull/2456))
- **Busco**
  - Fix barplot colors ([#2453](https://github.com/MultiQC/MultiQC/pull/2453))
- **Cell Ranger**
  - Fix parsing antibody tab without `antibody_treemap_plot` ([#2525](https://github.com/MultiQC/MultiQC/pull/2525))
- **Cutadapt**
  - Speed up module by caching parsing versions ([#2528](https://github.com/MultiQC/MultiQC/pull/2528))
- **DRAGEN**
  - Add ploidy estimation table ([#2496](https://github.com/MultiQC/MultiQC/pull/2496))
- **fastp**
  - When could not parse sample name from command (i.e. `stdin`), use filename and proceed ([#2536](https://github.com/MultiQC/MultiQC/pull/2536))
- **FastQC**
  - Skip per tile sequence quality section in FastQC reports for better performance ([#2552](https://github.com/MultiQC/MultiQC/pull/2552))
  - Fix a `ZeroDivisionError` error ([#2462](https://github.com/MultiQC/MultiQC/pull/2462))
  - Fix memory leak to make 7 times faster and use 10 times less memory ([#2552](https://github.com/MultiQC/MultiQC/pull/2552))
  - Do not keep intermediate data in memory to reduce memory footprint further ([#2516](https://github.com/MultiQC/MultiQC/pull/2516) )
  - Add option to ignore FastQC quality thresholds ([#2486](https://github.com/MultiQC/MultiQC/pull/2486))
- **goleft indexcov**
  - Work correctly even if no valid contigs in input ([#2540](https://github.com/MultiQC/MultiQC/pull/2540))
- **mosdepth**
  - Fix absolute coverage plot ([#2488](https://github.com/MultiQC/MultiQC/pull/2488))
- **nonpareil**
  - Change write_data_file label to be consistent with other modules ([#2472](https://github.com/MultiQC/MultiQC/pull/2472))
- **Picard**
  - WgsMetrics: coverage plot: show % based ≥x, not >x ([#2473](https://github.com/MultiQC/MultiQC/pull/2473))
  - CrosscheckFingerprints: support multiple files, preserve sample order in heatmap ([#2454](https://github.com/MultiQC/MultiQC/pull/2454))
- **qc3C**
  - Fix detecting sample name for relative path ([#2502](https://github.com/MultiQC/MultiQC/pull/2502))
- **QualiMap**
  - BamQC: when trimming long tails, keep at least 20x ([#2431](https://github.com/MultiQC/MultiQC/pull/2431))
- **Samtools**
  - Add support for `markdup` ([#2254](https://github.com/MultiQC/MultiQC/pull/2254))
  - Add violin multiple datasets & samtools flagstat percentage switch ([#2430](https://github.com/MultiQC/MultiQC/pull/2430))
- **Space Ranger**
  - fix for missing `genomic_dna` section ([#2429](https://github.com/MultiQC/MultiQC/pull/2429))
- **xengsort**
  - Fix parsing long files (do no use `content_lines`) ([#2484](https://github.com/MultiQC/MultiQC/pull/2484))

## [MultiQC v1.21](https://github.com/MultiQC/MultiQC/releases/tag/v1.21) - 2024-02-28

### Highlights

#### Box plot

Added a new plot type: box plot. It's useful to visualise a distribution when you have a set of values for each sample.

```py
from multiqc.plots import box
self.add_section(
    ...,
    plot=box.plot(
        {
            "sample 1": [4506, 4326, 3137, 1563, 1730, 3254, 2259, 3670, 2719, ...],
            "sample 2": [2145, 2011, 3368, 2132, 1673, 1993, 6635, 1635, 4984, ...],
            "sample 3": [1560, 1845, 3247, 1701, 2829, 2775, 3179, 1724, 1828, ...],
        },
        pconfig={
            "title": "Iso-Seq: Insert Length",
        },
    )
)
```

<img width="400" src="https://raw.githubusercontent.com/MultiQC/MultiQC/1e2ad5529d547ab9dc04c99274da49ad6a2e556f/docs/images/changelog/v1.21-boxplot.png">

Note the difference with the violin plot: the box plot visualises the distributions of many values within one sample, whereas the violin plot shows the distribution of one metric across many samples.

#### `pyproject.toml`

The `setup.py` file has been superseded by `pyproject.toml` for the build configuration.
Note that now for new modules, an entry point should be added to `pyproject.toml` instead of `setup.py`, e.g.:

```toml
[project.entry-points."multiqc.modules.v1"]
afterqc = "multiqc.modules.afterqc:MultiqcModule"
```

#### Heatmap

The heatmap plot now supports passing a dict as input data, and also supports a `zlab`
parameter to set the label for the z-axis:

```py
from multiqc.plots import heatmap
self.add_section(
    ...,
    plot=heatmap.plot(
        {
            "sample 1": {"sample 2": 0, "sample 3": 1},
            "sample 2": {"sample 1": 0, "sample 3": 0},
            "sample 3": {"sample 1": 1, "sample 2": 0, "sample 3": 1},
        },
        pconfig={
            "title": "Sample comparison",
            "zlab": "Match",
        },
    )
)
```

### MultiQC updates

- New plot type: box plot ([#2358](https://github.com/MultiQC/MultiQC/pull/2358))
- Add "Export to CSV" button for tables ([#2394](https://github.com/MultiQC/MultiQC/pull/2394))
- Replace `setup.py` with `pyproject.toml` ([#2353](https://github.com/MultiQC/MultiQC/pull/2353))
- Heatmap: allow a dict dicts of data ([#2386](https://github.com/MultiQC/MultiQC/pull/2386))
- Heatmap: add `zlab` config parameter. Show `xlab`, `ylab`, `zlab` in tooltip ([#2387](https://github.com/MultiQC/MultiQC/pull/2387))
- Warn if `run_modules` contains a non-existent module ([#2322](https://github.com/MultiQC/MultiQC/pull/2322))
- Catch non-hashable values (dicts, lists) passed as a table cell value ([#2348](https://github.com/MultiQC/MultiQC/pull/2348))
- Always create JSON even when MegaQC upload is disabled ([#2330](https://github.com/MultiQC/MultiQC/pull/2330))
- Use generic font family for Plotly ([#2368](https://github.com/MultiQC/MultiQC/pull/2368))
- Use a padded span with `nowrap` instead of `&nbsp;` before suffixes in table cells ([#2395](https://github.com/MultiQC/MultiQC/pull/2395))
- Refactor: fix unescaped regex strings ([#2384](https://github.com/MultiQC/MultiQC/pull/2384))

Fixes:

- Pin the required Plotly version and add a runtime version check ([#2325](https://github.com/MultiQC/MultiQC/pull/2325))
- Bar plot: preserve the sample order ([#2339](https://github.com/MultiQC/MultiQC/pull/2339))
- Bar plot: fix inner gap in group mode ([#2321](https://github.com/MultiQC/MultiQC/pull/2321))
- Violin: filter `Inf` values ([#2380](https://github.com/MultiQC/MultiQC/pull/2380))
- Table: Fix use of the `no_violin` (ex-`no_beeswarm`) table config flag ([#2376](https://github.com/MultiQC/MultiQC/pull/2376))
- Heatmap: prevent from parsing numerical sample names ([#2349](https://github.com/MultiQC/MultiQC/pull/2349))
- Work around call of `full_figure_for_development` to avoid Kaleido errors ([#2359](https://github.com/MultiQC/MultiQC/pull/2359))
- Auto-generate plot `id` when `pconfig=None` ([#2337](https://github.com/MultiQC/MultiQC/pull/2337))
- Fix: infinite `dmax` or `dmin` fail JSON dump load in JavaScript ([#2354](https://github.com/MultiQC/MultiQC/pull/2354))
- Fix: dump `pconfig` for MegaQC ([#2344](https://github.com/MultiQC/MultiQC/pull/2344))

### New modules

- [**IsoSeq**](https://github.com/PacificBiosciences/IsoSeq)
  - Iso-Seq contains the newest tools to identify transcripts in PacBio single-molecule sequencing data (HiFi reads). `cluster` and `refine` commands are supported.
- [**Space Ranger**](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger)
  - Works with data from 10X Genomics Visium. Processes sequencing reads and images created using
    the 10x Visium platform to generate count matrices with spatial information.
  - New MultiQC module parses Space Ranger quality reports.

### Module updates

- **bcl2fastq**: fix the top undetermined barcodes plot ([#2340](https://github.com/MultiQC/MultiQC/pull/2340))
- **DRAGEN**: add few coverage metrics in general stats ([#2341](https://github.com/MultiQC/MultiQC/pull/2341))
- **DRAGEN**: fix showing the number of found samples ([#2347](https://github.com/MultiQC/MultiQC/pull/2347))
- **DRAGEN**: support `gvcf_metrics` ([#2327](https://github.com/MultiQC/MultiQC/pull/2327))
- **fastp**: fix detection of JSON files ([#2334](https://github.com/MultiQC/MultiQC/pull/2334))
- **HTSeq Count**: robust file reading loop, ignore `.parquet` files ([#2364](https://github.com/MultiQC/MultiQC/pull/2364))
- **Illumina InterOp Statistics**: do not set `'scale': False` as a default ([#2350](https://github.com/MultiQC/MultiQC/pull/2350))
- **mosdepth**: fix regression in showing general stats ([#2346](https://github.com/MultiQC/MultiQC/pull/2346))
- **Picard**: Crosscheck Fingerprints updates ([#2388](https://github.com/MultiQC/MultiQC/pull/2388))
  - add a heatmap for LOD scores besides a table
  - if too many pairs in table, skip those with `Expected` status
  - use the `warn` status for `Inconclusive`
  - add a separate sample-wise table instead of general stats
  - sort tables by status, not by sample name
  - add a column "Best match" and "Best match LOD" in tables
  - hide the LOD Threshold column
- **PURPLE**: support v4.0.1 output without `version` column ([#2366](https://github.com/MultiQC/MultiQC/pull/2366))
- **Samtools**: support new `coverage` command ([#2356](https://github.com/MultiQC/MultiQC/pull/2356))
- **UMI-tools**: support new `extract` command ([#2296](https://github.com/MultiQC/MultiQC/pull/2296))
- **Whatshap**: make robust when a stdout is appended to TSV ([#2361](https://github.com/MultiQC/MultiQC/pull/2361))

## [MultiQC v1.20](https://github.com/MultiQC/MultiQC/releases/tag/v1.20) - 2024-02-12

### Highlights

#### New plotting library

MultiQC v1.20 comes with totally new plotting code for MultiQC reports. This is a huge change to the report output. We've done our best to maintain feature parity with the previous plotting code, but please do let us know if you spot any bugs or changes in behaviour by creating a GitHub issue.

This change comes with many improvements and new features, and paves the way for more in the future. To find out more, read the [associated blog post](https://seqera.io/blog/multiqc-plotly/).

For now, you can revert to the previous plotting code by using the `highcharts` report template (`multiqc --template highcharts`). This will be removed in v1.22.

Note that there are several plotting configuration options which have been removed:

- `click_func`
- `cursor`
- `tt_percentages` (use `tt_suffix: "%"`)
- Bar plot:
  - `use_legend` (automatically hidden if there is only 1 category)
- Line plot:
  - `labelSize`
  - `xDecimals`, `yDecimals` (automatic if all values can be cast to int)
  - `xLabelFormat`, `yLabelFormat` (use `tt_label`)
  - `pointFormat`
- Heatmap:
  - `datalabel_colour`
  - `borderWidth`

#### Moved GitHub and docker repositories

The v1.20 release is also the first release we've had since we moved the MultiQC repositories. Please note that the code is now at [MultiQC/MultiQC](https://github.com/MultiQC/MultiQC) (formerly [ewels/MultiQC](https://github.com/ewels/MultiQC)) and the same for the Docker repository. The GitHub repo should automatically redirect, but it's still good to update any references you may have.

### MultiQC updates

- Support Plotly as a new backend for plots ([#2079](https://github.com/MultiQC/MultiQC/pull/2079))
  - The default template now uses Plotly for all plots
  - Added a new plot type `violin` (replaces `beeswarm`)
  - Moved legacy Highcharts/Matplotlib code under an optional template `highcharts`
    ([#2292](https://github.com/MultiQC/MultiQC/pull/2292))
- Move GitHub repository to `MultiQC` organisation ([#2243](https://github.com/MultiQC/MultiQC/pull/2243))
- Update all GitHub actions to their latest versions ([#2242](https://github.com/MultiQC/MultiQC/pull/2242))
- Update docs to work with Astro 4 ([#2256](https://github.com/MultiQC/MultiQC/pull/2256))
- Remove unused dependency on `future` library ([#2258](https://github.com/MultiQC/MultiQC/pull/2258))
- Fix incorrect scale IDs caught by linting ([#2272](https://github.com/MultiQC/MultiQC/pull/2272))
- Docs: fix missing `v` prefix in docker image tags ([#2273](https://github.com/MultiQC/MultiQC/pull/2273))
- Unicode file reading errors: attempt to skip non-unicode characters ([#2275](https://github.com/MultiQC/MultiQC/pull/2275))
- Heatmap: check if value is numeric when calculating min and max ([#2276](https://github.com/MultiQC/MultiQC/pull/2276))
- Add `filesearch_file_shared` config option, remove unnecessary per-module `shared` flags in search patterns ([#2227](https://github.com/MultiQC/MultiQC/pull/2227))
- Use alternative method to walk directory using pathlib ([#2277](https://github.com/MultiQC/MultiQC/pull/2277))
- Export `config.output_dir` in MegaQC JSON ([#2287](https://github.com/MultiQC/MultiQC/pull/2287))
- Drop support for module tags ([#2278](https://github.com/MultiQC/MultiQC/pull/2278))
- Pin `Pillow` package, wrap add_logo in try-except ([#2312](https://github.com/MultiQC/MultiQC/pull/2312))
- Custom content: support multiple datasets ([#2291](https://github.com/MultiQC/MultiQC/pull/2291))
- Configuration: fix reading config.output_fn_name and --filename ([#2314](https://github.com/MultiQC/MultiQC/pull/2314))

### New modules

- [**Bamdst**](https://https://github.com/shiquan/bamdst) ([#2161](https://github.com/MultiQC/MultiQC/pull/2161))
  - Bamdst is a lightweight tool to stat the depth coverage of target regions of bam file(s).
- [**MetaPhlAn**](https://github.com/biobakery/MetaPhlAn) ([#2262](https://github.com/MultiQC/MultiQC/pull/2262))
  - MetaPhlAn is a computational tool for profiling the composition of microbial communities from metagenomic shotgun sequencing data.
- [**MEGAHIT**](https://github.com/voutcn/megahit) ([#2222](https://github.com/MultiQC/MultiQC/pull/2222))
  - MEGAHIT is an ultra-fast and memory-efficient NGS assembler
- [**Nonpareil**](https://github.com/lmrodriguezr/nonpareil) ([#2215](https://github.com/MultiQC/MultiQC/pull/2215))
  - Estimate metagenomic coverage and sequence diversity.

### Module updates

- **Bcftools**: order variant depths plot categories ([#2289](https://github.com/MultiQC/MultiQC/pull/2289))
- **Bcftools**: add missing `self.ignore_samples` in stats ([#2288](https://github.com/MultiQC/MultiQC/pull/2288))
- **BCL Convert**: add index, project names to sample statistics and calculate mean quality for lane statistics. ([#2261](https://github.com/MultiQC/MultiQC/pull/2261))
- **BCL Convert**: fix duplicated `yield` for 3.9.3+ when the yield is provided explicitly in Quality_Metrics ([#2253](https://github.com/MultiQC/MultiQC/pull/2253))
- **BCL Convert**: handle samples with zero yield ([#2297](https://github.com/MultiQC/MultiQC/pull/2297))
- **Bismark**: fix old link in docs ([#2252](https://github.com/MultiQC/MultiQC/pull/2252))
- **Cutadapt**: support JSON format ([#2281](https://github.com/MultiQC/MultiQC/pull/2281))
- **HiFiasm**: account for lines with no asterisk ([#2268](https://github.com/MultiQC/MultiQC/pull/2268))
- **HUMID**: add cluster statistics ([#2265](https://github.com/MultiQC/MultiQC/pull/2265))
- **mosdepth**: add additional summaries to general stats #2257 ([#2257](https://github.com/MultiQC/MultiQC/pull/2257))
- **Picard**: fix using multiple times in report: do not pass `module.anchor` to `self.find_log_files` ([#2255](https://github.com/MultiQC/MultiQC/pull/2255))
- **QualiMap**: address NBSP as thousands separators ([#2282](https://github.com/MultiQC/MultiQC/pull/2282))
- **Seqera Platform CLI**: updates for v0.9.2 ([#2248](https://github.com/MultiQC/MultiQC/pull/2248))
- **Seqera Platform CLI**: handle failed tasks ([#2286](https://github.com/MultiQC/MultiQC/pull/2286))

## [MultiQC v1.19](https://github.com/MultiQC/MultiQC/releases/tag/v1.19) - 2023-12-18

### MultiQC updates

- Add missing table `id` in DRAGEN modules, and require `id` in plot configs in strict mode ([#2228](https://github.com/MultiQC/MultiQC/pull/2228))
- Config `table_columns_visible` and `table_columns_name`: support flat config and `table_id` as a group ([#2191](https://github.com/MultiQC/MultiQC/pull/2191))
- Add `sort_samples: false` config option for bar graphs ([#2210](https://github.com/MultiQC/MultiQC/pull/2210))
- Upgrade the jQuery tablesorter plugin to v2 ([#1666](https://github.com/MultiQC/MultiQC/pull/1666))
- Refactor pre-Python-3.6 code, prefer f-strings over `.format()` calls ([#2224](https://github.com/MultiQC/MultiQC/pull/2224))
- Allow specifying default sort columns for tables with `defaultsort` ([#1667](https://github.com/MultiQC/MultiQC/pull/1667))
- Create CODE_OF_CONDUCT.md ([#2195](https://github.com/MultiQC/MultiQC/pull/2195))
- Add `.cram` to sample name cleaning defaults ([#2209](https://github.com/MultiQC/MultiQC/pull/2209))

### MultiQC bug fixes

- Re-add `run` into the `multiqc` namespace ([#2202](https://github.com/MultiQC/MultiQC/pull/2202))
- Fix the `"square": True` flag to scatter plot to actually make the plot square ([#2189](https://github.com/MultiQC/MultiQC/pull/2189))
- Fix running with the `--no-report` flag ([#2212](https://github.com/MultiQC/MultiQC/pull/2212))
- Fix guessing custom content plot type: do not assume first row of a bar plot data are sample names ([#2208](https://github.com/MultiQC/MultiQC/pull/2208))
- Fix detection of changed specific module in Changelog CI ([#2234](https://github.com/MultiQC/MultiQC/pull/2234))

### Module updates

- **BCLConvert**: fix mean quality, fix count-per-lane bar plot ([#2197](https://github.com/MultiQC/MultiQC/pull/2197))
- **deepTools**: handle missing data in `plotProfile` ([#2229](https://github.com/MultiQC/MultiQC/pull/2229))
- **Fastp**: search content instead of file name ([#2213](https://github.com/MultiQC/MultiQC/pull/2213))
- **GATK**: square the `BaseRecalibrator` scatter plot ([#2189](https://github.com/MultiQC/MultiQC/pull/2189))
- **HiC-Pro**: add missing search patterns and better handling of missing data ([#2233](https://github.com/MultiQC/MultiQC/pull/2233))
- **Kraken**: fix `UnboundLocalError` ([#2230](https://github.com/MultiQC/MultiQC/pull/2230))
- **Kraken**: fixed column keys in genstats ([#2205](https://github.com/MultiQC/MultiQC/pull/2205))
- **QualiMap**: fix `BamQC` for global-only stats ([#2207](https://github.com/MultiQC/MultiQC/pull/2207))
- **Picard**: add more search patterns for `MarkDuplicates`, including `MarkDuplicatesSpark` ([#2226](https://github.com/MultiQC/MultiQC/pull/2226))
- **Salmon**: add `library_types`, `compatible_fragment_ratio`, `strand_mapping_bias` to the general stats table ([#1485](https://github.com/MultiQC/MultiQC/pull/1485))

## [MultiQC v1.18](https://github.com/MultiQC/MultiQC/releases/tag/v1.18) - 2023-11-17

### Highlights

#### Better configs

As of this release, you can now set all of your config variables via environment variables! (see [docs](https://multiqc.info/docs/getting_started/config/#config-with-environment-variables)).

Better still, YAML config files can now use string interpolation to parse environment variables within strings (see [docs](https://multiqc.info/docs/getting_started/config/#referencing-environment-variables-in-yaml-configs)), eg:

```yaml
report_header_info:
  - Contact E-mail: !ENV "${NAME:info}@${DOMAIN:example.com}"
```

#### Picard refactoring

In this release, there was a significant refactoring of the Picard module.
It has been generalized for better code sharing with other Picard-based software, like Sentieon and Parabricks.
As a result of this, the standalone Sentieon module was removed: Sentieon QC files will be interpreted directly as Picard QC files.

If you were using the Sentieon module in your pipelines, make sure to update any places that reference the module name:

- MultiQC command line (e.g. replace `--module sentieon` with `--module picard`).
- MultiQC configs (e.g. replace `sentieon` with `picard` in options like `run_modules`, `exclude_modules`, `module_order`).
- Downstream code that relies on names of the files in `multiqc_data` or `multiqc_plots` saves (e.g., `multiqc_data/multiqc_sentieon_AlignmentSummaryMetrics.txt` becomes `multiqc_data/multiqc_picard_AlignmentSummaryMetrics.txt`).
- Code that parses data files like `multiqc_data/multiqc_data.json`.
- Custom plugins and templates that rely on HTML anchors (e.g. `#sentieon_aligned_reads` becomes `#picard_AlignmentSummaryMetrics`).
- Also, note that Picard fetches sample names from the commands it finds inside the QC headers (e.g. `# net.sf.picard.analysis.CollectMultipleMetrics INPUT=Szabo_160930_SN583_0215_AC9H20ACXX.bam ...` -> `Szabo_160930_SN583_0215_AC9H20ACXX`), whereas the removed Sentieon module prioritized the QC file names. To revert to the old Sentieon approach, use the [`use_filename_as_sample_name` config flag](https://multiqc.info/docs/getting_started/config/#using-log-filenames-as-sample-names).

### MultiQC updates

- Config can be set with environment variables, including env var interpolation ([#2178](https://github.com/MultiQC/MultiQC/pull/2178))
- Try find config in `~/.config` or `$XDG_CONFIG_HOME` ([#2183](https://github.com/MultiQC/MultiQC/pull/2183))
- Better sample name cleaning with pairs of input filenames ([#2181](https://github.com/MultiQC/MultiQC/pull/2181))
- Software versions: allow any string as a version tag ([#2166](https://github.com/MultiQC/MultiQC/pull/2166))
- Table columns with non-numeric values and now trigger a linting error if `scale` is set ([#2176](https://github.com/MultiQC/MultiQC/pull/2176))
- Stricter config variable typing ([#2178](https://github.com/MultiQC/MultiQC/pull/2178))
- Remove `position:absolute` CSS from table values ([#2169](https://github.com/MultiQC/MultiQC/pull/2169))
- Fix column sorting in exported TSV files from a matplotlib linegraph plot ([#2143](https://github.com/MultiQC/MultiQC/pull/2143))
- Fix custom anchors for kraken ([#2170](https://github.com/MultiQC/MultiQC/pull/2170))
- Fix logging spillover bug ([#2174](https://github.com/MultiQC/MultiQC/pull/2174))

### New Modules

- [**Seqera Platform CLI**](https://github.com/seqeralabs/tower-cli) ([#2151](https://github.com/MultiQC/MultiQC/pull/2151))
  - Seqera Platform CLI reports statistics generated by the Seqera Platform CLI.
- [**Xenome**](https://github.com/data61/gossamer/blob/master/docs/xenome.md) ([#1860](https://github.com/MultiQC/MultiQC/pull/1860))
  - A tool for classifying reads from xenograft sources.
- [**xengsort**](https://gitlab.com/genomeinformatics/xengsort) ([#2168](https://github.com/MultiQC/MultiQC/pull/2168))
  - xengsort is a fast xenograft read sorter based on space-efficient k-mer hashing

### Module updates

- **fastp**: add version parsing ([#2159](https://github.com/MultiQC/MultiQC/pull/2159))
- **fastp**: correctly parse sample name from `--in1`/`--in2` in bash command. Prefer file name if not `fastp.json`; fallback to file name when error ([#2139](https://github.com/MultiQC/MultiQC/pull/2139))
- **Kaiju**: fix `division by zero` error ([#2179](https://github.com/MultiQC/MultiQC/pull/2179))
- **Nanostat**: account for both tab and spaces in `v1.41+` search pattern ([#2155](https://github.com/MultiQC/MultiQC/pull/2155))
- **Pangolin**: update for v4: add QC Note , update tool versions columns ([#2157](https://github.com/MultiQC/MultiQC/pull/2157))
- **Picard**: Generalize to directly support Sentieon and Parabricks outputs ([#2110](https://github.com/MultiQC/MultiQC/pull/2110))
- **Sentieon**: Removed the module in favour of directly supporting parsing by the **Picard** module ([#2110](https://github.com/MultiQC/MultiQC/pull/2110))
  - Note that any code that relies on the module name needs to be updated, e.g. `-m sentieon` will no longer work
  - The exported plot and data files will be now be prefixed as `picard` instead of `sentieon`, etc.
  - Note that the Sentieon module used to fetch the sample names from the file names by default, and now it follows the Picard module's logic, and prioritizes the commands recorded in the logs. To override, use the `use_filename_as_sample_name` config flag

## [MultiQC v1.17](https://github.com/MultiQC/MultiQC/releases/tag/v1.17) - 2023-10-17

### The one with the new logo

Highlights:

- Introducing the new MultiQC logo!
- Adding support for Python 3.12 and dropping support for Python 3.7
- New `--require-logs` to fail if expected tool outputs are not found
- Rename `--lint` to `--strict`
- Modules should now use `ModuleNotFoundError` instead of `UserWarning` when no logs are found
- 2 new modules and updates to 9 modules.

### MultiQC updates

- Add CI action [changelog.yml](.github%2Fworkflows%2Fchangelog.yml) to populate the changelog from PR titles, triggered by a comment `@multiqc-bot changelog` ([#2025](https://github.com/MultiQC/MultiQC/pull/2025), [#2102](https://github.com/MultiQC/MultiQC/pull/2102), [#2115](https://github.com/MultiQC/MultiQC/pull/2115))
- Add GitHub Actions bot workflow to fix code linting from a PR comment ([#2082](https://github.com/MultiQC/MultiQC/pull/2082))
- Use custom exception type instead of `UserWarning` when no samples are found. ([#2049](https://github.com/MultiQC/MultiQC/pull/2049))
- Lint modules for missing `self.add_software_version` ([#2081](https://github.com/MultiQC/MultiQC/pull/2081))
- Strict mode: rename `config.lint` to `config.strict`, crash early on module or template error. Add `MULTIQC_STRICT=1` ([#2101](https://github.com/MultiQC/MultiQC/pull/2101))
- Matplotlib line plots now respect `xLog: True` and `yLog: True` in config ([#1632](https://github.com/MultiQC/MultiQC/pull/1632))
- Fix matplotlib linegraph and bargraph for the case when `xmax` `<` `xmin` in config ([#2124](https://github.com/MultiQC/MultiQC/pull/2124))
- Add `--require-logs` flag to error out if requested modules not used ([#2109](https://github.com/MultiQC/MultiQC/pull/2109))
- Fixes for python 3.12
  - Replace removed `distutils` ([#2113](https://github.com/MultiQC/MultiQC/pull/2113))
  - Bundle lzstring ([#2119](https://github.com/MultiQC/MultiQC/pull/2119))
- Drop Python 3.6 and 3.7 support, add 3.12 ([#2121](https://github.com/MultiQC/MultiQC/pull/2121))
- Just run CI on the oldest + newest supported Python versions ([#2074](https://github.com/MultiQC/MultiQC/pull/2074))
- <img src="./multiqc/templates/default/assets/img/favicon-16x16.png" alt="///" width="10px"/> New logo
- Set name and anchor for the custom content "module" [#2131](https://github.com/MultiQC/MultiQC/pull/2131)
- Fix use of `shutil.copytree` when overriding existing template files in `tmp_dir` ([#2133](https://github.com/MultiQC/MultiQC/pull/2133))

### New Modules

- [**Bracken**](https://ccb.jhu.edu/software/bracken/)
  - A highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample.
- [**Truvari**](https://github.com/ACEnglish/truvari) ([#1751](https://github.com/MultiQC/MultiQC/pull/1751))
  - Truvari is a toolkit for benchmarking, merging, and annotating structural variants

### Module updates

- **Dragen**: make sure all inputs are recorded in multiqc_sources.txt ([#2128](https://github.com/MultiQC/MultiQC/pull/2128))
- **Cellranger**: Count submodule updated to parse Antibody Capture summary ([#2118](https://github.com/MultiQC/MultiQC/pull/2118))
- **fastp**: parse unescaped sample names with white spaces ([#2108](https://github.com/MultiQC/MultiQC/pull/2108))
- **FastQC**: Add top overrepresented sequences table ([#2075](https://github.com/MultiQC/MultiQC/pull/2075))
- **HiCPro**: Fix parsing scientific notation in hicpro-ashic. Thanks @Just-Roma ([#2126](https://github.com/MultiQC/MultiQC/pull/2126))
- **HTSeq Count**: allow counts files with more than 2 columns ([#2129](https://github.com/MultiQC/MultiQC/pull/2129))
- **mosdepth**: fix prioritizing region over global information ([#2106](https://github.com/MultiQC/MultiQC/pull/2106))
- **Picard**: Adapt WgsMetrics to parabricks bammetrics outputs ([#2127](https://github.com/MultiQC/MultiQC/pull/2127))
- **Picard**: MarkDuplicates: Fix parsing mixed strings/numbers, account for missing trailing `0` ([#2083](https://github.com/MultiQC/MultiQC/pull/2083), [#2094](https://github.com/MultiQC/MultiQC/pull/2094))
- **Samtools**: Add MQ0 reads to the Percent Mapped barplot in Stats submodule ([#2123](https://github.com/MultiQC/MultiQC/pull/2123))
- **WhatsHap**: Process truncated input with no ALL chromosome ([#2095](https://github.com/MultiQC/MultiQC/pull/2095))

## [MultiQC v1.16](https://github.com/MultiQC/MultiQC/releases/tag/v1.16) - 2023-09-22

### Highlight: Software versions

New in v1.16 - software version information can now automatically parsed from log output where available,
and added to MultiQC in a standardised manner. It's shown in the MultiQC report next to section
headings and in a dedicated report section, as well as being saved to `multiqc_data`.
Where version information is not available in logs, it can be submitted manually by using a new
special file type with filename pattern `*_mqc_versions.yml`.
There's the option of representing groups of versions, useful for a tool that uses sub-tools,
or pipelines that want to report version numbers per analysis step.

There are a handful of new config scopes to control behaviour:
`software_versions`, `skip_versions_section`, `disable_version_detection`, `versions_table_group_header`.
See the documentation for more ([writing modules](https://multiqc.info/docs/development/modules/#saving-version-information),
[supplying stand-alone](https://multiqc.info/docs/reports/customisation/#listing-software-versions))

Huge thanks to [@pontushojer](https://https://github.com/pontushojer) for the
contribution ([#1927](https://github.com/MultiQC/MultiQC/pull/1927)).
This idea goes way back to [issue #290](https://github.com/MultiQC/MultiQC/issues/290), made in 2016!

### MultiQC updates

- Removed `simplejson` unused dependency ([#1973](https://github.com/MultiQC/MultiQC/pull/1973))
- Give config `custom_plot_config` priority over column-specific settings set by modules
- When exporting plots, make a more clear error message for unsupported FastQC dot plot ([#1976](https://github.com/MultiQC/MultiQC/pull/1976))
- Fixed parsing of `plot_type: "html"` `data` in json custom content
- Replace deprecated `pkg_resources`
- [Fix](https://github.com/MultiQC/MultiQC/issues/2032) the module groups configuration for modules where the namespace is passed explicitly to `general_stats_addcols`. Namespace is now always appended to the module name in the general stats ([2037](https://github.com/MultiQC/MultiQC/pull/2037)).
- Do not call `sys.exit()` in the `multiqc.run()` function, to avoid breaking interactive environments. [#2055](https://github.com/MultiQC/MultiQC/pull/2055)
- Fixed the DOI exports in `multiqc_data` to include more than just the MultiQC paper ([#2058](https://github.com/MultiQC/MultiQC/pull/2058))
- Fix table column color scaling then there are negative numbers ([1869](https://github.com/MultiQC/MultiQC/issues/1869))
- Export plots as static images and data in a ZIP archive. Fixes the [issue](https://github.com/MultiQC/MultiQC/issues/1873) when only 10 plots maximum were downloaded due to the browser limitation.

### New Modules

- [**Bakta**](https://github.com/oschwengers/bakta)
  - Rapid and standardized annotation of bacterial genomes, MAGs & plasmids.
- [**mapDamage**](https://github.com/ginolhac/mapDamage)
  - mapDamage2 is a computational framework written in Python and R, which tracks and quantifies DNA damage patterns among ancient DNA sequencing reads generated by Next-Generation Sequencing platforms.
- [**Sourmash**](https://github.com/sourmash-bio/sourmash)
  - Quickly search, compare, and analyze genomic and metagenomic data sets.

### Module updates

- **BcfTools**
  - Stats: fix parsing multi-sample logs ([#2052](https://github.com/MultiQC/MultiQC/issues/2052))
- **Custom content**
  - Don't convert sample IDs to floats ([#1883](https://github.com/MultiQC/MultiQC/issues/1883))
- **DRAGEN**
  - Make DRAGEN module use `fn_clean_exts` instead of hardcoded file names. Fixes working with arbitrary file names ([#1994](https://github.com/MultiQC/MultiQC/pull/1994))
- **FastQC**:
  - fix `UnicodeDecodeError` when parsing `fastqc_data.txt`: try latin-1 or fail gracefully ([#2024](https://github.com/MultiQC/MultiQC/issues/2024))
- **Kaiju**:
  - Fix `UnboundLocalError` on outputs when Kanju was run with the `-e` flag ([#2023](https://github.com/MultiQC/MultiQC/pull/2023))
- **Kraken**
  - Parametrize top-N through config ([#2060](https://github.com/MultiQC/MultiQC/pull/2060))
  - Fix bug where ranks incorrectly assigned to tabs ([#1766](https://github.com/MultiQC/MultiQC/issues/1766)).
- **Mosdepth**
  - Add X/Y relative coverage plot, analogous to the one in samtools-idxstats ([#1978](https://github.com/MultiQC/MultiQC/issues/1978))
  - Added the `perchrom_fraction_cutoff` option into the config to help avoid clutter in contig-level plots
  - Fix a bug happening when both `region` and `global` coverage histograms for a sample are available (i.e. when mosdepth was run with `--by`, see [mosdepth docs](https://github.com/brentp/mosdepth#usage)). In this case, data was effectively merged. Instead, summarise it separately and add a separate report section for the region-based coverage data.
  - Do not fail when all input samples have no coverage ([#2005](https://github.com/MultiQC/MultiQC/pull/2005)).
- **NanoStat**
  - Support new format ([#1997](https://github.com/MultiQC/MultiQC/pull/1997)).
- **RSeQC**
  - Fix `max() arg is an empty sequence` error ([#1985](https://github.com/MultiQC/MultiQC/issues/1985))
  - Fix division by zero on all-zero input ([#2040](https://github.com/MultiQC/MultiQC/pull/2040))
- **Samtools**
  - Stats: fix "Percent Mapped" plot when samtools was run with read filtering ([#1972](https://github.com/MultiQC/MultiQC/pull/1972))
- **Qualimap**
  - BamQC: Include `% On Target` in General Stats table ([#2019](https://github.com/MultiQC/MultiQC/issues/2019))
- **WhatsHap**
  - Bugfix: ensure that TSV is only split on tab character. Allows sample names with spaces ([#1981](https://github.com/MultiQC/MultiQC/pull/1981))

## [MultiQC v1.15](https://github.com/MultiQC/MultiQC/releases/tag/v1.15) - 2023-08-04

#### Potential breaking change in some edge cases

This release of MultiQC introduces speed improvements to the file search.
One way it does this is by limiting the number of lines loaded by each search pattern.
For the vast majority of users, this should have no effect except faster searches.
However, in some edge cases it may break things. Hypothetically, for example:

- If you concatenate log files from multiple tools
- If you have a custom plugin module that we haven't tested

See the [troubleshooting docs](https://multiqc.info/docs/usage/troubleshooting/#long-log-files)
for more information.

### MultiQC updates

- Refactor file search for performance improvements ([#1904](https://github.com/MultiQC/MultiQC/pull/1904))
- Bump `log_filesize_limit` default (to skip large files in the search) from 10MB to 50MB.
- Table code now tolerates lambda function calls with bad data ([#1739](https://github.com/MultiQC/MultiQC/issues/1739))
- Beeswarm plot now saves data to `multiqc_data`, same as tables ([#1861](https://github.com/MultiQC/MultiQC/issues/1861))
- Don't print DOI in module if it's set to an empty string.
- Optimize parsing of 2D data dictionaries in multiqc.utils.utils_functions.write_data_file() ([#1891](https://github.com/MultiQC/MultiQC/pull/1891))
- Don't sort table headers alphabetically if we don't have an `OrderedDict` - regular dicts are fine in Py3 ([#1866](https://github.com/MultiQC/MultiQC/issues/1866))
- New back-end to preview + deploy the new website when the docs are edited.
- Fixed a _lot_ of broken links in the documentation from the new website change in structure.

### New Modules

- [**Freyja**](https://github.com/andersen-lab/Freyja)
  - Freyja is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset.
- [**Librarian**](https://github.com/DesmondWillowbrook/Librarian)
  - A tool to predict the sequencing library type from the base composition of a supplied FastQ file.

### Module updates

- **Conpair**
  - Bugfix: allow to find and proprely parse the `concordance` output of Conpair, which may output 2 kinds of format for `concordance` depending if it's ran with or without `--outfile` ([#1851](https://github.com/MultiQC/MultiQC/issues/1851))
- **Cell Ranger**
  - Bugfix: avoid multiple `KeyError` exceptions when parsing Cell Ranger 7.x `web_summary.html` ([#1853](https://github.com/MultiQC/MultiQC/issues/1853), [#1871](https://github.com/MultiQC/MultiQC/issues/1871))
- **DRAGEN**
  - Restored functionality to show target BED coverage metrics ([#1844](https://github.com/MultiQC/MultiQC/issues/1844))
  - Update filename pattern in RNA quant metrics ([#1958](https://github.com/MultiQC/MultiQC/pull/1958))
- **filtlong**
  - Handle reports from locales that use `.` as a thousands separator ([#1843](https://github.com/MultiQC/MultiQC/issues/1843))
- **GATK**
  - Adds support for [AnalyzeSaturationMutagenesis submodule](https://gatk.broadinstitute.org/hc/en-us/articles/360037594771-AnalyzeSaturationMutagenesis-BETA-)
- **HUMID**
  - Fix bug that prevent HUMID stats files from being parsed ([#1856](https://github.com/MultiQC/MultiQC/issues/1856))
- **Mosdepth**
  - Fix data not written to `mosdepth_cumcov_dist.txt` and `mosdepth_cov_dist.txt` ([#1868](https://github.com/MultiQC/MultiQC/issues/1868))
  - Update documentation with new file `{prefix}.mosdepth.summary.txt` ([#1868](https://github.com/MultiQC/MultiQC/issues/1868))
  - Fill in missing values for general stats table ([#1868](https://github.com/MultiQC/MultiQC/issues/1868))
  - Include mosdepth/summary file paths in `multiqc_sources.txt` ([#1868](https://github.com/MultiQC/MultiQC/issues/1868))
  - Enable log switch for Coverage per contig plot ([#1868](https://github.com/MultiQC/MultiQC/issues/1868))
  - Fix y-axis scaling for Coverage distribution plot ([#1868](https://github.com/MultiQC/MultiQC/issues/1868))
  - Handle case of intermediate missing coverage x-values in the `*_dist.txt` file causing a distorted Coverage distribution plot ([#1960](https://github.com/MultiQC/MultiQC/issues/1960))
- **Picard**
  - WgsMetrics: Fix wrong column label ([#1888](https://github.com/MultiQC/MultiQC/issues/1888))
  - HsMetrics: Add missing field descriptions ([#1928](https://github.com/MultiQC/MultiQC/pull/1928))
- **Porechop**
  - Don't render bar graphs if no samples had any adapters trimmed ([#1850](https://github.com/MultiQC/MultiQC/issues/1850))
  - Added report section listing samples that had no adapters trimmed
- **RSeQC**
  - Fix `ZeroDivisionError` error for `bam_stat` results when there are 0 reads ([#1735](https://github.com/MultiQC/MultiQC/issues/1735))
- **UMI-tools**
  - Fix bug that broke the module with paired-end data ([#1845](https://github.com/MultiQC/MultiQC/issues/1845))

## [MultiQC v1.14](https://github.com/MultiQC/MultiQC/releases/tag/v1.14) - 2023-01-08

### MultiQC new features

- Rewrote the `Dockerfile` to build multi-arch images (amd64 + arm), run through a non-privileged user and build tools for non precompiled python binaries ([#1541](https://github.com/MultiQC/MultiQC/pull/1541), [#1541](https://github.com/MultiQC/MultiQC/pull/1541))
- Add a new lint test to check that colour scale names are valid ([#1835](https://github.com/MultiQC/MultiQC/pull/1835))
- Update github actions to run tests on a single module if it is the only file affected by the PR ([#915](https://github.com/MultiQC/MultiQC/issues/915))
- Add CI testing for Python 3.10 and 3.11
- Optimize line-graph generation to remove an n^2 loop ([#1668](https://github.com/MultiQC/MultiQC/pull/1668))
- Parsing output file column headers is much faster.

### MultiQC code cleanup

- Remove Python 2-3 compatability `from __future__` imports
- Remove unused `#!/usr/bin/env python` hashbangs from module files
- Add new code formatting tool [isort](https://pycqa.github.io/isort/) to standardise the order and formatting of Python module imports
- Add [Pycln](https://hadialqattan.github.io/pycln/#/) pre-commit hook to remove unused imports

### MultiQC updates

- Bugfix: Make `config.data_format` work again ([#1722](https://github.com/MultiQC/MultiQC/issues/1722))
- Bump minimum version of Jinja2 to `>=3.0.0` ([#1642](https://github.com/MultiQC/MultiQC/issues/1642))
- Disable search progress bar if running with `--quiet` or `--no-ansi` ([#1638](https://github.com/MultiQC/MultiQC/issues/1638))
- Allow path filters without full paths by trying to prefix analysis dir when filtering ([#1308](https://github.com/MultiQC/MultiQC/issues/1308))
- Fix sorting of table columns with text values
- Don't crash if a barplot is given an empty list of categories ([#1540](https://github.com/MultiQC/MultiQC/issues/1540))
- New logos! MultiQC is now developed and maintained at [Seqera Labs](https://seqera.io/). Updated logos and email addresses accordingly.

### New Modules

- [**Anglerfish**](https://github.com/remiolsen/anglerfish)
  - A tool designed to assess pool balancing, contamination and insert sizes of Illumina library dry runs on Oxford Nanopore data.
- [**BBDuk**](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
  - Combines most common data-quality-related trimming, filtering, and masking operations via kmers into a single high-performance tool.
- [**Cell Ranger**](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
  - Works with data from 10X Genomics Chromium. Processes Chromium single cell data to align reads, generate feature-barcode matrices, perform clustering and other secondary analysis, and more.
  - New MultiQC module parses Cell Ranger quality reports from VDJ and count analysis
- [**DIAMOND**](https://github.com/bbuchfink/diamond)
  - A high-throughput program for aligning DNA reads or protein sequences against a protein reference database.
- [**DRAGEN-FastQC**](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)
  - Illumina Bio-IT Platform that uses FPGA for accelerated primary and secondary analysis
  - Finally merged the epic 2.5-year-old pull request, with 3.5k new lines of code.
  - Please report any bugs you find!
- [**Filtlong**](https://github.com/rrwick/Filtlong)
  - A tool for filtering long reads by quality.
- [**GoPeaks**](https://github.com/maxsonBraunLab/gopeaks)
  - GoPeaks is used to call peaks in CUT&TAG/CUT&RUN datasets.
- [**HiFiasm**](https://github.com/chhylp123/hifiasm)
  - A haplotype-resolved assembler for accurate Hifi reads
- [**HUMID**](https://github.com/jfjlaros/dedup)
  - HUMID is a tool to quickly and easily remove duplicate reads from FastQ files, with or without UMIs.
- [**mOTUs**](https://motu-tool.org/)
  - Microbial profiling through marker gene (MG)-based operational taxonomic units (mOTUs)
- [**Nextclade**](https://github.com/nextstrain/nextclade)
  - Tool that assigns clades to SARS-CoV-2 samples
- [**Porechop**](https://github.com/rrwick/Porechop)
  - A tool for finding and removing adapters from Oxford Nanopore reads
- [**PRINSEQ++**](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus)
  - PRINSEQ++ is a C++ of `prinseq-lite.pl` program for filtering, reformating or trimming genomic and metagenomic sequence data.
- [**UMI-tools**](https://umi-tools.readthedocs.io)
  - Work with Unique Molecular Identifiers (UMIs) / Random Molecular Tags (RMTs) and single cell RNA-Seq cell barcodes.

### Module updates

- **Bcftools stats**
  - Bugfix: Do not show empty bcftools stats variant depth plots ([#1777](https://github.com/MultiQC/MultiQC/pull/1777))
  - Bugfix: Avoid exception when `PSC nMissing` column is not present ([#1832](https://github.com/MultiQC/MultiQC/issues/1832))
- **BCL Convert**
  - Handle single-end read data correctly when setting cluster length instead of always assuming paired-end reads ([#1697](https://github.com/MultiQC/MultiQC/issues/1697))
  - Handle different R1 and R2 read-lengths correctly instead of assuming they are the same ([#1774](https://github.com/MultiQC/MultiQC/issues/1774))
  - Handle single-index paired-end data correctly
  - Added a config option to enable the creation of barplots with undetermined barcodes (`create_unknown_barcode_barplots` with `False` as default) ([#1709](https://github.com/MultiQC/MultiQC/pull/1709))
- **BUSCO**
  - Update BUSCO pass/warning/fail scheme to be more clear for users
- **Bustools**
  - Show median reads per barcode statistic
- **Custom content**
  - Create a report even if there's only Custom Content General Stats there
  - Attempt to cooerce line / scatter x-axes into floats so as not to lose labels ([#1242](https://github.com/MultiQC/MultiQC/issues/1242))
  - Multi-sample line-graph TSV files that have no sample name in row 1 column 1 now use row 1 as x-axis labels ([#1242](https://github.com/MultiQC/MultiQC/issues/1242))
- **fastp**
  - Add total read count (after filtering) to general stats table ([#1744](https://github.com/MultiQC/MultiQC/issues/1744))
  - Don't crash for invalid JSON files ([#1652](https://github.com/MultiQC/MultiQC/issues/1652))
- **FastQC**
  - Report median read-length for fastqc in addition to mean ([#1745](https://github.com/MultiQC/MultiQC/pull/1745))
- **Kaiju**
  - Don't crash if we don't have any data for the top-5 barplot ([#1540](https://github.com/MultiQC/MultiQC/issues/1540))
- **Kallisto**
  - Fix `ZeroDivisionError` when a sample has zero reads ([#1746](https://github.com/MultiQC/MultiQC/issues/1746))
- **Kraken**
  - Fix duplicate heatmap to account for missing taxons ([#1779](https://github.com/MultiQC/MultiQC/pull/1779))
  - Make heatmap full width
  - Handle empty kreports gracefully ([#1637](https://github.com/MultiQC/MultiQC/issues/1637))
  - Fix regex error with very large numbers of unclassified reads ([#1639](https://github.com/MultiQC/MultiQC/pull/1639))
- **malt**
  - Fixed division by 0 in malt module ([#1683](https://github.com/MultiQC/MultiQC/issues/1683))
- **miRTop**
  - Avoid `KeyError` - don't assume all fields present in logs ([#1778](https://github.com/MultiQC/MultiQC/issues/1778))
- **Mosdepth**
  - Don't pad the General Stats table with zeros for missing data ([#1810](https://github.com/MultiQC/MultiQC/pull/1810))
- **Picard**
  - HsMetrics: Allow custom columns in General Stats too, with `HsMetrics_genstats_table_cols` and `HsMetrics_genstats_table_cols_hidden`
- **Qualimap**
  - Added additional columns in general stats for BamQC results that displays region on-target stats if region bed has been supplied (hidden by default) ([#1798](https://github.com/MultiQC/MultiQC/pull/1798))
  - Bugfix: Remove General Stats rows for filtered samples ([#1780](https://github.com/MultiQC/MultiQC/issues/1780))
- **RSeQC**
  - Update `geneBody_coverage` to plot normalized coverages using a similar formula to that used by RSeQC itself ([#1792](https://github.com/MultiQC/MultiQC/pull/1792))
- **Sambamba Markdup**
  - Catch zero division in sambamba markdup ([#1654](https://github.com/MultiQC/MultiQC/issues/1654))
- **Samtools**
  - Added additional column for `flagstat` that displays percentage of mapped reads in a bam (hidden by default) ([#1733](https://github.com/MultiQC/MultiQC/issues/1733))
- **VEP**
  - Don't crash with `ValueError` if there are zero variants ([#1681](https://github.com/MultiQC/MultiQC/issues/1681))

## [MultiQC v1.13](https://github.com/MultiQC/MultiQC/releases/tag/v1.13) - 2022-09-08

### MultiQC updates

- Major spruce of the command line help, using the new [rich-click](https://github.com/ewels/rich-click) package
- Drop some of the Python 2k compatability code (eg. custom requirements)
- Improvements for running MultiQC in a Python environment, such as a Jupyter Notebook or script
  - Fixed bug raised when removing logging file handlers between calls that arose when configuring the root logger with dictConfig ([#1643](https://github.com/MultiQC/MultiQC/issues/1643))
- Added new config option `custom_table_header_config` to override any config for any table header
- Fixed edge-case bug in custom content where a `description` that doesn't terminate in `.` gave duplicate section descriptions.
- Tidied the verbose log to remove some very noisy statements and add summaries for skipped files in the search
- Add timezone to time in reports
- Add nix flake support
- Added automatic tweet about new releases
- Breaking: Removed `--cl_config` option. Please use `--cl-config` instead.

### Module updates

- **AdapterRemoval**
  - Finally merge a fix for counts of reads that are discarded/collapsed ([#1647](https://github.com/MultiQC/MultiQC/issues/1647))
- **VEP**
  - Fixed bug when `General Statistics` have a value of `-` ([#1656](https://github.com/MultiQC/MultiQC/pull/1656))
- **Custom content**
  - Only set id for custom content when id not set by metadata ([#1629](https://github.com/MultiQC/MultiQC/issues/1629))
  - Fix bug where module wouldn't run if all content was within a MultiQC config file ([#1686](https://github.com/MultiQC/MultiQC/issues/1686))
  - Fix crash when `info` isn't set ([#1688](https://github.com/MultiQC/MultiQC/issues/1688))
- **Nanostat**
  - Removed HTML escaping of special characters in the log to fix bug in jinja2 v3.10 removing `jinja2.escape()` ([#1659](https://github.com/MultiQC/MultiQC/pull/1659))
  - Fix bug where module would crash if input does not contain quality scores ([#1717](https://github.com/MultiQC/MultiQC/issues/1717))
- **Pangolin**
  - Updated module to handle outputs from Pangolin v4 ([#1660](https://github.com/MultiQC/MultiQC/pull/1660))
- **Somalier**
  - Handle zero mean X depth in _Sex_ plot ([#1670](https://github.com/MultiQC/MultiQC/pull/1670))
- **Fastp**
  - Include low complexity and too long reads in filtering bar chart
- **miRTop**
  - Fix module crashing when `ref_miRNA_sum` is missing in file. ([#1712](https://github.com/MultiQC/MultiQC/issues/1712))
  - Fix module crashing due to zero division ([#1719](https://github.com/MultiQC/MultiQC/issues/1719))
- **FastQC**
  - Fixed error when parsing duplicate ratio when there is `nan` values in the report. ([#1725](https://github.com/MultiQC/MultiQC/pull/1725))

## [MultiQC v1.12](https://github.com/MultiQC/MultiQC/releases/tag/v1.12) - 2022-02-08

### MultiQC - new features

- Added option to customise default plot height in plot config ([#1432](https://github.com/MultiQC/MultiQC/issues/1432))
- Added `--no-report` flag to skip report generation ([#1462](https://github.com/MultiQC/MultiQC/issues/1462))
- Added support for priting tool DOI in report sections ([#1177](https://github.com/MultiQC/MultiQC/issues/1177))
- Added support for `--custom-css-file` / `config.custom_css_files` option to include custom CSS in the final report ([#1573](https://github.com/MultiQC/MultiQC/pull/1573))
- New plot config option `labelSize` to customise font size for axis labels in flat MatPlotLib charts ([#1576](https://github.com/MultiQC/MultiQC/pull/1576))
- Added support for customising table column names ([#1255](https://github.com/MultiQC/MultiQC/issues/1255))

### MultiQC - updates

- MultiQC now skips modules for which no files were found - gives a small performance boost ([#1463](https://github.com/MultiQC/MultiQC/issues/1463))
- Improvements for running MultiQC in a Python environment, such as a Jupyter Notebook or script
  - Fixed logger bugs when calling `multiqc.run` multiple times by removing logging file handlers between calls ([#1141](https://github.com/MultiQC/MultiQC/issues/1141))
  - Init/reset global state between runs ([#1596](https://github.com/MultiQC/MultiQC/pull/1596))
- Added commonly missing functions to several modules ([#1468](https://github.com/MultiQC/MultiQC/issues/1468))
- Wrote new script to check for the above function calls that should be in every module (`.github/workflows/code_checks.py`), runs on GitHub actions CI
- Make table _Conditional Formatting_ work at table level as well as column level. ([#761](https://github.com/MultiQC/MultiQC/issues/761))
- CSS Improvements to make printed reports more attractive / readable ([#1579](https://github.com/MultiQC/MultiQC/pull/1579))
- Fixed a problem with numeric filenames ([#1606](https://github.com/MultiQC/MultiQC/issues/1606))
- Fixed nasty bug where line charts with a categorical x-axis would take categories from the last sample only ([#1568](https://github.com/MultiQC/MultiQC/issues/1568))
- Ignore any files called `multiqc_data.json` ([#1598](https://github.com/MultiQC/MultiQC/issues/1598))
- Check that the config `path_filters` is a list, convert to list if a string is supplied ([#1539](https://github.com/MultiQC/MultiQC/issues/1539))

### New Modules

- [**CheckQC**](https://github.com/Molmed/checkQC)
  - A program designed to check a set of quality criteria against an Illumina runfolder
- [**pbmarkdup**](https://github.com/PacificBiosciences/pbmarkdup)
  - Mark duplicate reads from PacBio sequencing of an amplified library
- [**WhatsHap**](https://whatshap.readthedocs.io)
  - WhatsHap is a software for phasing genomic variants using DNA sequencing reads
- [**SeqWho**](https://daehwankimlab.github.io/seqwho/)
  - Tool to determine a FASTQ(A) sequencing file identity, both source protocol and species of origin.

### Module feature additions

- **BBMap**
  - Added handling for `qchist` output ([#1021](https://github.com/MultiQC/MultiQC/issues/1021))
- **bcftools**
  - Added a plot with samplewise number of sites, Ts/Tv, number of singletons and sequencing depth ([#1087](https://github.com/MultiQC/MultiQC/issues/1087))
- **Mosdepth**
  - Added mean coverage [#1566](https://github.com/MultiQC/MultiQC/issues/1566)
- **NanoStat**
  - Recognize FASTA and FastQ report flavors ([#1547](https://github.com/MultiQC/MultiQC/issues/1547))

### Module updates

- **BBMap**
  - Correctly handle adapter stats files with additional columns ([#1556](https://github.com/MultiQC/MultiQC/issues/1556))
- **BCL Convert**
  - Handle change in output format in v3.9.3 with new `Quality_Metrics.csv` file ([#1563](https://github.com/MultiQC/MultiQC/issues/1563))
- **bowtie**
  - Minor update to handle new log wording in bowtie v1.3.0 ([#1615](https://github.com/MultiQC/MultiQC/issues/1615))
- **CCS**
  - Tolerate compound IDs generated by pbcromwell ccs in the general statistics ([#1486](https://github.com/MultiQC/MultiQC/pull/1486))
  - Fix report parsing. Update test on attributes ids ([#1583](https://github.com/MultiQC/MultiQC/issues/1583))
- **Custom content**
  - Fixed module failing when writing data to file if there is a `/` in the section name ([#1515](https://github.com/MultiQC/MultiQC/issues/1515))
  - Use filename for section header in files with no headers ([#1550](https://github.com/MultiQC/MultiQC/issues/1550))
  - Sort custom content bargraph data by default ([#1412](https://github.com/MultiQC/MultiQC/issues/1412))
  - Always save `custom content` data to file with a name reflecting the section name. ([#1194](https://github.com/MultiQC/MultiQC/issues/1194))
- **DRAGEN**
  - Fixed bug in sample name regular expression ([#1537](https://github.com/MultiQC/MultiQC/pull/1537))
- **Fastp**
  - Fixed % pass filter statistics ([#1574](https://github.com/MultiQC/MultiQC/issues/1574))
- **FastQC**
  - Fixed several bugs occuring when FastQC sections are skipped ([#1488](https://github.com/MultiQC/MultiQC/issues/1488), [#1533](https://github.com/MultiQC/MultiQC/issues/1533))
  - Clarify general statistics table header for length
- **goleft/indexcov**
  - Fix `ZeroDivisionError` if no bins are found ([#1586](https://github.com/MultiQC/MultiQC/issues/1586))
- **HiCPro**
  - Better handling of errors when expected data keys are not found ([#1366](https://github.com/MultiQC/MultiQC/issues/1366))
- **Lima**
  - Move samples that have been renamed using `--replace-names` into the _General Statistics_ table ([#1483](https://github.com/MultiQC/MultiQC/pull/1483))
- **miRTrace**
  - Replace hardcoded RGB colours with Hex to avoid errors with newer versions of matplotlib ([#1263](https://github.com/MultiQC/MultiQC/pull/1263))
- **Mosdepth**
  - Fixed issue [#1568](https://github.com/MultiQC/MultiQC/issues/1568)
  - Fixed a bug when reporting per contig coverage
- **Picard**
  - Update `ExtractIlluminaBarcodes` to recognise more log patterns in newer versions of Picard ([#1611](https://github.com/MultiQC/MultiQC/pull/1611))
- **Qualimap**
  - Fix `ZeroDivisionError` in `QM_RNASeq` and skip genomic origins plot if no aligned reads are found ([#1492](https://github.com/MultiQC/MultiQC/issues/1492))
- **QUAST**
  - Clarify general statistics table header for length
- **RSeQC**
  - Fixed minor bug in new TIN parsing where the sample name was not being correctly cleaned ([#1484](https://github.com/MultiQC/MultiQC/issues/1484))
  - Fixed bug in the `junction_saturation` submodule ([#1582](https://github.com/MultiQC/MultiQC/issues/1582))
  - Fixed bug where empty files caused `tin` submodule to crash ([#1604](https://github.com/MultiQC/MultiQC/issues/1604))
  - Fix bug in `read_distribution` for samples with zero tags ([#1571](https://github.com/MultiQC/MultiQC/issues/1571))
- **Sambamba**
  - Fixed issue with a change in the format of output from `sambamba markdup` 0.8.1 ([#1617](https://github.com/MultiQC/MultiQC/issues/1617))
- **Skewer**
  - Fix `ZeroDivisionError` if no available reads are found ([#1622](https://github.com/MultiQC/MultiQC/issues/1622))
- **Somalier**
  - Plot scaled X depth instead of mean for _Sex_ plot ([#1546](https://github.com/MultiQC/MultiQC/issues/1546))
- **VEP**
  - Handle table cells containing `-` instead of numbers ([#1597](https://github.com/MultiQC/MultiQC/issues/1597))

## [MultiQC v1.11](https://github.com/MultiQC/MultiQC/releases/tag/v1.11) - 2021-07-05

### MultiQC new features

- New interactive slider controls for controlling heatmap colour scales ([#1427](https://github.com/MultiQC/MultiQC/issues/1427))
- Added new `--replace-names` / config `sample_names_replace` option to replace sample names during report generation
- Added `use_filename_as_sample_name` config option / `--fn_as_s_name` command line flag ([#949](https://github.com/MultiQC/MultiQC/issues/949), [#890](https://github.com/MultiQC/MultiQC/issues/890), [#864](https://github.com/MultiQC/MultiQC/issues/864), [#998](https://github.com/MultiQC/MultiQC/issues/998), [#1390](https://github.com/MultiQC/MultiQC/issues/1390))
  - Forces modules to use the log filename for the sample identifier, even if the module usually takes this from the file contents
  - Required a change to the `clean_s_name()` function arguments. All core MultiQC modules updated to reflect this.
  - Should be backwards compatible for custom modules. To adopt new behaviour, supply `f` instead of `f["root"]` as the second argument.
  - See the documenation for details: [Using log filenames as sample names](https://multiqc.info/docs/#using-log-filenames-as-sample-names) and [Custom sample names](https://multiqc.info/docs/#custom-sample-names).

### MultiQC updates

- Make the module crash tracebacks much prettier using `rich`
- Refine the cli log output a little (nicely formatted header line + drop the `[INFO]`)
- Added docs describing tools for downstream analysis of MultiQC outputs.
- Added CI tests for Python 3.9, pinned `networkx` package to `>=2.5.1` ([#1413](https://github.com/MultiQC/MultiQC/issues/1413))
- Added patterns to `config.fn_ignore_paths` to avoid error with parsing installation dir / singularity cache ([#1416](https://github.com/MultiQC/MultiQC/issues/1416))
- Print a log message when flat-image plots are used due to sample size surpassing `plots_flat_numseries` config ([#1254](https://github.com/MultiQC/MultiQC/issues/1254))
- Fix the `mqc_colours` util function to lighten colours even when passing categorical or single-length lists.
- Bugfix for Custom Content, using YAML configuration (eg. section headers) for images should now work

### New Modules

- [**BCL Convert**](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html)
  - Tool that converts / demultiplexes Illumina Binary Base Call (BCL) files to FASTQ files
- [**Bustools**](https://bustools.github.io/)
  - Tools for working with BUS files
- [**ccs**](https://github.com/PacificBiosciences/ccs)
  - Generate highly accurate single-molecule consensus reads from PacBio data
- [**GffCompare**](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
  - GffCompare can annotate and estimate accuracy of one or more GFF files compared with a reference annotation.
- [**Lima**](https://github.com/PacificBiosciences/barcoding)
  - The PacBio Barcode Demultiplexer
- [**NanoStat**](https://github.com/wdecoster/nanostat)
  - Calculate various statistics from a long read sequencing datasets
- [**ODGI**](https://github.com/pangenome/odgi)
  - Optimized dynamic genome/graph implementation
- [**Pangolin**](https://github.com/cov-lineages/pangolin)
  - Added MultiQC support for Pangolin, the tool that determines SARS-CoV-2 lineages
- [**Sambamba Markdup**](https://lomereiter.github.io/sambamba/docs/sambamba-markdup.html)
  - Added MultiQC module to add duplicate rate calculated by Sambamba Markdup.
- [**Snippy**](https://github.com/tseemann/snippy)
  - Rapid haploid variant calling and core genome alignment.
- [**VEP**](https://www.ensembl.org/info/docs/tools/vep/index.html)
  - Added MultiQC module to add summary statistics of Ensembl VEP annotations.
  - Handle error from missing variants in VEP stats file. ([#1446](https://github.com/MultiQC/MultiQC/issues/1446))

### Module feature additions

- **Cutadapt**
  - Added support for linked adapters [#1329](https://github.com/MultiQC/MultiQC/issues/1329)]
  - Parse whether trimming was 5' or 3' for _Lengths of Trimmed Sequences_ plot where possible
- **Mosdepth**
  - Include or exclude contigs based on patterns for coverage-per-contig plots
- **Picard**
  - Add support for `CollectIlluminaBasecallingMetrics`, `CollectIlluminaLaneMetrics`, `ExtractIlluminaBarcodes` and `MarkIlluminaAdapters` ([#1336](https://github.com/MultiQC/MultiQC/pull/1336))
  - New `insertsize_xmax` configuration option to limit the plotted maximum insert size for `InsertSizeMetrics`
- **Qualimap**
  - Added new percentage coverage plot in `QM_RNASeq` ([#1258](https://github.com/MultiQC/MultiQC/issues/1258))
- **RSeQC**
  - Added a long-requested submodule to support showing the [**TIN**](http://rseqc.sourceforge.net/#tin-py) (Transcript Integrity Number) ([#737](https://github.com/MultiQC/MultiQC/issues/737))

### Module updates

- **biscuit**
  - Duplicate Rate and Cytosine Retention tables are now bargraphs.
  - Refactor code to only calculate alignment statistics once.
  - Fixed bug where cytosine retentions values would not be properly read if in scientific notation.
- **bcl2fastq**
  - Added sample name cleaning so that prepending directories with the `-d` flag works properly.
- **Cutadapt**
  - Plot filtered reads even when no filtering category is found ([#1328](https://github.com/MultiQC/MultiQC/issues/1328))
  - Don't take the last command line string for the sample name if it looks like a command-line flag ([#949](https://github.com/MultiQC/MultiQC/issues/949))
- **Dragen**
  - Handled MultiQC crashing when run on single-end output from Dragen ([#1374](https://github.com/MultiQC/MultiQC/issues/1374))
- **fastp**
  - Handle a `ZeroDivisionError` if there are zero reads ([#1444](https://github.com/MultiQC/MultiQC/issues/1444))
- **FastQC**
  - Added check for if `overrepresented_sequences` is missing from reports ([#1281](https://github.com/MultiQC/MultiQC/issues/1444))
- **Flexbar**
  - Fixed bug where reports with 0 reads would crash MultiQC ([#1407](https://github.com/MultiQC/MultiQC/issues/1407))
- **Kraken**
  - Handle a `ZeroDivisionError` if there are zero reads ([#1440](https://github.com/MultiQC/MultiQC/issues/1440))
  - Updated search patterns to handle edge case ([#1428](https://github.com/MultiQC/MultiQC/issues/1428))
- **Mosdepth**
  - Show barplot instead of line graph for coverage-per-contig plot if there is only one contig.
- **Picard**
  - `RnaSeqMetrics` - fix assignment barplot labels to say bases instead of reads ([#1408](https://github.com/MultiQC/MultiQC/issues/1408))
  - `CrosscheckFingerprints` - fix bug where LOD threshold was not detected when invoked with "new" picard cli style. fixed formatting bug ([#1414](https://github.com/MultiQC/MultiQC/issues/1414))
  - Made checker for comma as decimal separator in `HsMetrics` more robust ([#1296](https://github.com/MultiQC/MultiQC/issues/1296))
- **qc3C**
  - Updated module to not fail on older field names.
- **Qualimap**
  - Fixed wrong units in tool tip label ([#1258](https://github.com/MultiQC/MultiQC/issues/1258))
- **QUAST**
  - Fixed typo causing wrong number of contigs being displayed ([#1442](https://github.com/MultiQC/MultiQC/issues/1442))
- **Sentieon**
  - Handled `ZeroDivisionError` when input files have zero reads ([#1420](https://github.com/MultiQC/MultiQC/issues/1420))
- **RSEM**
  - Handled `ZeroDivisionError` when input files have zero reads ([#1040](https://github.com/MultiQC/MultiQC/issues/1040))
- **RSeQC**
  - Fixed double counting of some categories in `read_distribution` bar graph. ([#1457](https://github.com/MultiQC/MultiQC/issues/1457))

## [MultiQC v1.10.1](https://github.com/MultiQC/MultiQC/releases/tag/v1.10.1) - 2021-04-01

### MultiQC updates

- Dropped the `Skipping search pattern` log message from a warning to debug
- Moved directory prepending with `-d` back to before sample name cleaning (as it was before v1.7) ([#1264](https://github.com/MultiQC/MultiQC/issues/1264))
- If linegraph plot data goes above `ymax`, only _discard_ the data if the line doesn't come back again ([#1257](https://github.com/MultiQC/MultiQC/issues/1257))
- Allow scientific notation numbers in colour scheme generation
  - Fixed bug with very small minimum numbers that only revelead itself after a bugfix done in the v1.10 release
- Allow `top_modules` to be specified as empty dicts ([#1274](https://github.com/MultiQC/MultiQC/issues/1274))
- Require at least `rich` version `9.4.0` to avoid `SpinnerColumn` `AttributeError` ([#1393](https://github.com/MultiQC/MultiQC/issues/1393))
- Properly ignore `.snakemake` folders as intended ([#1395](https://github.com/MultiQC/MultiQC/issues/1395))

#### Module updates

- **bcftools**
  - Fixed bug where `QUAL` value `.` would crash MultiQC ([#1400](https://github.com/MultiQC/MultiQC/issues/1400))
- **bowtie2**
  - Fix bug where HiSAT2 paired-end bar plots were missing unaligned reads ([#1230](https://github.com/MultiQC/MultiQC/issues/1230))
- **Deeptools**
  - Handle `plotProfile` data where no upstream / downstream regions have been calculated around genes ([#1317](https://github.com/MultiQC/MultiQC/issues/1317))
  - Fix `IndexError` caused by mysterious `-1` in code.. ([#1275](https://github.com/MultiQC/MultiQC/issues/1275))
- **FastQC**
  - Replace `NaN` with `0` in the _Per Base Sequence Content_ plot to avoid crashing the plot ([#1246](https://github.com/MultiQC/MultiQC/issues/1246))
- **Picard**
  - Fixed bug in `ValidateSamFile` module where additional whitespace at the end of the file would cause MultiQC to crash ([#1397](https://github.com/MultiQC/MultiQC/issues/1397))
- **Somalier**
  - Fixed bug where using sample name cleaning in a config would trigger a `KeyError` ([#1234](https://github.com/MultiQC/MultiQC/issues/1234))

## [MultiQC v1.10](https://github.com/MultiQC/MultiQC/releases/tag/v1.10) - 2021-03-08

### Update for developers: Code linting

This is a big change for MultiQC developers. I have added automated code formatting and code linting
(style checks) to MultiQC. This helps to keep the MultiQC code base consistent despite having many
contributors and helps me to review pull-requests without having to consider whitespace.

Specifically, MultiQC now uses three main tools:

- [Black](https://github.com/psf/black) - Python Code
- [Prettier](https://prettier.io/) - Everything else (almost)
- [markdownlint-cli](https://github.com/igorshubovych/markdownlint-cli) - Stricter markdown rules

**All developers must run these tools when submitting changes via Pull-Requests!**
Automated CI tests now run with GitHub actions to check that all files pass the above tests.
If any files do not, that test will fail giving a red :x: next to the pull request.

For further information, please see the [documentation](https://multiqc.info/docs/#coding-with-multiqc).

### MultiQC updates

#### New MultiQC Features

- `--sample-filters` now also accepts `show_re` and `hide_re` in addition to `show` and `hide`. The `_re` options use regex, while the "normal" options use globbing.
- MultiQC config files now work with `.yml` file extension as well as `.yaml`
  - `.yaml` will take preference if both found.
- Section comments can now also be added for _General Statistics_
  - `section_comments: { general_stats: "My comment" }`
- New table header config option `bgcols` allows background colours for table cells with categorical data.
- New table header config options `cond_formatting_rules` and `cond_formatting_colours`
  - Comparable functionality to user config options `table_cond_formatting_rules` and `table_cond_formatting_colours`,
    allowes module developers to format table cell values as labels.
- New CI test looks for git merge markers in files
- Beautiful new [progress bar](https://rich.readthedocs.io/en/stable/progress.html) from the amazing [willmcgugan/rich](https://github.com/willmcgugan/rich) package.
- Added a bunch of new default sample name trimming suffixes ([see `8ac5c7b`](https://github.com/MultiQC/MultiQC/commit/8ac5c7b6e4ea6003ca2c9b681953ab3f22c5dd66))
- Added `timeout-minutes: 10` to the CI test workflow to check that changes aren't negatively affecting run time too much.
- New table header option `bars_zero_centrepoint` to treat `0` as zero width bars and plot bar length based on absolute values

#### New Modules

- [**EigenStratDatabaseTools**](https://github.com/TCLamnidis/EigenStratDatabaseTools)
  - Added MultiQC module to report SNP coverages from `eigenstrat_snp_coverage.py` in the general stats table.
- [**HOPS**](https://www.github.com/rhubler/HOPS)
  - Post-alignment ancient DNA analysis tool for MALT
- [**JCVI**](https://github.com/tanghaibao/jcvi)
  - Computes statistics on genome annotation.
- [**ngsderive**](https://github.com/stjudecloud/ngsderive)
  - Forensic analysis tool useful in backwards computing information from next-generation sequencing data.
- [**OptiType**](https://github.com/FRED-2/OptiType)
  - Precision HLA typing from next-generation sequencing data
- [**PURPLE**](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator)
  - A purity, ploidy and copy number estimator for whole genome tumor data
- [**Pychopper**](https://github.com/nanoporetech/pychopper)
  - Identify, orient and trim full length Nanopore cDNA reads
- [**qc3C**](https://github.com/cerebis/qc3C)
  - Reference-free QC of Hi-C sequencing data
- [**Sentieon**](https://www.sentieon.com/products/)
  - Submodules added to catch Picard-based QC metrics files

#### Module updates

- **DRAGEN**
  - Fix issue where missing out fields could crash the module ([#1223](https://github.com/MultiQC/MultiQC/issues/1223))
  - Added support for whole-exome / targetted data ([#1290](https://github.com/MultiQC/MultiQC/issues/1290))
- **featureCounts**
  - Add support for output from [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) ([#1022](https://github.com/MultiQC/MultiQC/issues/1022))
- **fgbio**
  - Fix `ErrorRateByReadPosition` to calculate `ymax` not just on the overall `error_rate`, but also specific base errors (ex. `a_to_c_error_rate`, `a_to_g_error_rate`, ...). ([#1215](https://github.com/MultiQC/MultiQC/pull/1251))
  - Fix `ErrorRateByReadPosition` plotted line names to no longer concatenate multiple read identifiers and no longer have off-by-one read numbering (ex. `Sample1_R2_R3` -> `Sample1_R2`) ([#[1304](https://github.com/MultiQC/MultiQC/pull/1304))
- **Fastp**
  - Fixed description for duplication rate (pre-filtering, not post) ([#[1313](https://github.com/MultiQC/MultiQC/pull/1313))
- **GATK**
  - Add support for the creation of a "Reported vs Empirical Quality" graph to the Base Recalibration module.
- **hap.py**
  - Updated module to plot both SNP and INDEL stats ([#1241](https://github.com/MultiQC/MultiQC/issues/1241))
- **indexcov**
  - Fixed bug when making the PED file plots ([#1265](https://github.com/MultiQC/MultiQC/issues/1265))
- **interop**
  - Added the `% Occupied` metric to `Read Metrics per Lane` table which is reported for NovaSeq and iSeq platforms.
- **Kaiju**
  - Fixed bug affecting inputs with taxa levels other than Phylum ([#1217](https://github.com/MultiQC/MultiQC/issues/1217))
  - Rework barplot, add top 5 taxons ([#1219](https://github.com/MultiQC/MultiQC/issues/1219))
- **Kraken**
  - Fix `ZeroDivisionError` ([#1276](https://github.com/MultiQC/MultiQC/issues/1276))
  - Add distinct minimizer heatmap for KrakenUniq style duplication information ([#1333](https://github.com/MultiQC/MultiQC/pull/1380))
- **MALT**
  - Fix y-axis labelling in bargraphs
- **MACS2**
  - Add number of peaks to the _General Statistics_ table.
- **mosdepth**
  - Enable prepending of directory to sample names
  - Display contig names in _Coverage per contig_ plot tooltip
- **Picard**
  - Fix `HsMetrics` bait percentage columns ([#1212](https://github.com/MultiQC/MultiQC/issues/1212))
  - Fix `ConvertSequencingArtifactToOxoG` files not being found ([#1310](https://github.com/MultiQC/MultiQC/issues/1310))
  - Make `WgsMetrics` histogram smoothed if more than 1000 data points (avoids huge plots that crash the browser)
  - Multiple new config options for `WgsMetrics` to customise coverage histogram and speed up MultiQC with very high coverage files.
  - Add additional datasets to Picard Alignment Summary ([#1293](https://github.com/MultiQC/MultiQC/issues/1293))
  - Add support for `CrosscheckFingerprints` ([#1327](https://github.com/MultiQC/MultiQC/issues/1327))
- **PycoQC**
  - Log10 x-axis for _Read Length_ plot ([#1214](https://github.com/MultiQC/MultiQC/issues/1214))
- **Rockhopper**
  - Fix issue with parsing genome names in Rockhopper summary files ([#1333](https://github.com/MultiQC/MultiQC/issues/1333))
  - Fix issue properly parsing multiple samples within a single Rockhopper summary file
- **Salmon**
  - Only try to generate a plot for fragment length if the data was found.
- **verifyBamID**
  - Fix `CHIP` value detection ([#1316](https://github.com/MultiQC/MultiQC/pull/1316)).

#### New Custom Content features

- General Stats custom content now gives a log message
- If `id` is not set in `JSON` or `YAML` files, it defaults to the sample name instead of just `custom_content`
- Data from `JSON` or `YAML` now has `data` keys (sample names) run through the `clean_s_name()` function to apply sample name cleanup
- Fixed minor bug which caused custom content YAML files with a string `data` type to not be parsed

#### Bug Fixes

- Disable preservation of timestamps / modes when copying temp report files, to help issues with network shares ([#1333](https://github.com/MultiQC/MultiQC/issues/1333))
- Fixed MatPlotLib warning: `FixedFormatter should only be used together with FixedLocator`
- Fixed long-standing min/max bug with shared minimum values for table columns using `shared_key`
- Made table colour schemes work with negative numbers (don't strip `-` from values when making scheme)

## [MultiQC v1.9](https://github.com/MultiQC/MultiQC/releases/tag/v1.9) - 2020-05-30

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

- Now using [GitHub Actions](https://github.com/features/actions) for all CI testing
  - Dropped Travis and AppVeyor, everything is now just on GitHub
  - Still testing on both Linux and Windows, with multiple versions of Python
  - CI tests should now run automatically for anyone who forks the MultiQC repository
- Linting with `--lint` now checks line graphs as well as bar graphs
- New `gathered` template with no tool name sections ([#1119](https://github.com/MultiQC/MultiQC/issues/1119))
- Added `--sample-filters` option to add _show_/_hide_ buttons at the top of the report ([#1125](https://github.com/MultiQC/MultiQC/issues/1125))
  - Buttons control the report toolbox Show/Hide tool, filtering your samples
  - Allows reports to be pre-configured based on a supplied list of sample names at report-generation time.
- Line graphs can now have `Log10` buttons (same functionality as bar graphs)
- Importing and running `multiqc` in a script is now a _little_ Better
  - `multiqc.run` now returns the `report` and `config` as well as the exit code. This means that you can explore the MultiQC run time a little in the Python environment.
  - Much more refactoring is needed to make MultiQC as useful in Python scripts as it could be. Watch this space.
- If a custom module `anchor` is set using `module_order`, it's now used a bit more:
  - Prefixed to module _section_ IDs
  - Appended to files saved in `multiqc_data`
  - Should help to prevent duplicates requiring `-1` suffixes when running a module multiple times
- New heatmap plot config options `xcats_samples` and `ycats_samples`
  - If set to `False`, the report toolbox options (_highlight_, _rename_, _show/hide_) do not affect that axis.
  - Means that the _Show only matching samples_ report toolbox option works on FastQC Status Checks, for example ([#1172](https://github.com/MultiQC/MultiQC/issues/1172))
- Report header time and analysis paths can now be hidden
  - New config options `show_analysis_paths` and `show_analysis_time` ([#1113](https://github.com/MultiQC/MultiQC/issues/1113))
- New search pattern key `skip: true` to skip specific searches when modules look for a lot of different files (eg. Picard).
- New `--profile-runtime` command line option (`config.profile_runtime`) to give analysis of how long the report takes to be generated
  - Plots of the file search results and durations are added to the end of the MultiQC report as a special module called _Run Time_
  - A summary of the time taken for the major stages of MultiQC execution are printed to the command line log.
- New table config option `only_defined_headers`
  - Defaults to `true`, set to `false` to also show any data columns that are not defined as headers
  - Useful as allows table-wide defaults to be set with column-specific overrides
- New `module` key allowed for `config.extra_fn_clean_exts` and `config.fn_clean_exts`
  - Means you can limit the action of a sample name cleaning pattern to specific MultiQC modules ([#905](https://github.com/MultiQC/MultiQC/issues/905))

#### New Custom Content features

- Improve support for HTML files - now just end your HTML filename with `_mqc.html`
  - Native handling of HTML snippets as files, no MultiQC config or YAML file required.
  - Also with embedded custom content configuration at the start of the file as a HTML comment.
- Add ability to group custom-content files into report sections
  - Use the new `parent_id`, `parent_name` and `parent_description` config keys to group content together like a regular module ([#1008](https://github.com/MultiQC/MultiQC/issues/1008))
- Custom Content files can now be configured using `custom_data`, without giving search patterns or data
  - Allows you to set descriptions and nicer titles for images and other 'blunt' data types in reports ([#1026](https://github.com/MultiQC/MultiQC/issues/1026))
  - Allows configuration of custom content separately from files themselves (`tsv`, `csv`, `txt` formats) ([#1205](https://github.com/MultiQC/MultiQC/issues/1205))

#### New Modules

- [**DRAGEN**](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)
  - Illumina Bio-IT Platform that uses FPGA for secondary NGS analysis
- [**iVar**](https://github.com/andersen-lab/ivar)
  - Added support for iVar: a computational package that contains functions broadly useful for viral amplicon-based sequencing.
- [**Kaiju**](http://kaiju.binf.ku.dk/)
  - Fast and sensitive taxonomic classification for metagenomics
- [**Kraken**](https://ccb.jhu.edu/software/kraken2/)
  - K-mer matching tool for taxonomic classification. Module plots bargraph of counts for top-5 hits across each taxa rank. General stats summary.
- [**MALT**](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/malt/)
  - Megan Alignment Tool: Metagenomics alignment tool.
- [**miRTop**](https://github.com/miRTop/mirtop)
  - Command line tool to annotate miRNAs with a standard mirna/isomir naming (mirGFF3)
  - Module started by [@oneillkza](https://github.com/oneillkza/) and completed by [@FlorianThibord](https://github.com/FlorianThibord/)
- [**MultiVCFAnalyzer**](https://github.com/alexherbig/multivcfanalyzer)
  - Combining multiple VCF files into one coherent report and format for downstream analysis.
- **Picard** - new submodules for `QualityByCycleMetrics`, `QualityScoreDistributionMetrics` & `QualityYieldMetrics`
  - See [#1116](https://github.com/MultiQC/MultiQC/issues/1114)
- [**Rockhopper**](https://cs.wellesley.edu/~btjaden/Rockhopper/)
  - RNA-seq tool for bacteria, includes bar plot showing where features map.
- [**Sickle**](https://github.com/najoshi/sickle)
  - A windowed adaptive trimming tool for FASTQ files using quality
- [**Somalier**](https://github.com/brentp/somalier)
  - Relatedness checking and QC for BAM/CRAM/VCF for cancer, DNA, BS-Seq, exome, etc.
- [**VarScan2**](https://github.com/dkoboldt/varscan)
  - Variant calling and somatic mutation/CNV detection for next-generation sequencing data

#### Module updates

- **BISCUIT**
  - Major rewrite to work with new BISCUIT QC script (BISCUIT `v0.3.16+`)
    - This change breaks backwards-compatability with previous BISCUIT versions. If you are unable to upgrade BISCUIT, please use MultiQC v1.8.
  - Fixed error when missing data in log files ([#1101](https://github.com/MultiQC/MultiQC/issues/1101))
- **bcl2fastq**
  - Samples with multiple library preps (i.e barcodes) will now be handled correctly ([#1094](https://github.com/MultiQC/MultiQC/issues/1094))
- **BUSCO**
  - Updated log search pattern to match new format in v4 with auto-lineage detection option ([#1163](https://github.com/MultiQC/MultiQC/issues/1163))
- **Cutadapt**
  - New bar plot showing the proportion of reads filtered out for different criteria (eg. _too short_, _too many Ns_) ([#1198](https://github.com/MultiQC/MultiQC/issues/1198))
- **DamageProfiler**
  - Removes redundant typo in init name. This makes referring to the module's column consistent with other modules when customising general stats table.
- **DeDup**
  - Updates plots to make compatible with 0.12.6
  - Fixes reporting errors - barplot total represents _mapped_ reads, not total reads in BAM file
  - New: Adds 'Post-DeDup Mapped Reads' column to general stats table.
- **FastQC**
  - Fixed tooltip text in _Sequence Duplication Levels_ plot ([#1092](https://github.com/MultiQC/MultiQC/issues/1092))
  - Handle edge-case where a FastQC report was for an empty file with 0 reads ([#1129](https://github.com/MultiQC/MultiQC/issues/1129))
- **FastQ Screen**
  - Don't skip plotting `% No Hits` even if it's `0%` ([#1126](https://github.com/MultiQC/MultiQC/issues/1126))
  - Refactor parsing code. Avoids error with `-0.00 %Unmapped` ([#1126](https://github.com/MultiQC/MultiQC/issues/1126))
  - New plot for _Bisulfite Reads_, if data is present
  - Categories in main plot are now sorted by the total read count and hidden if 0 across all samples
- **fgbio**
  - New: Plot error rate by read position from `ErrorRateByReadPosition`
  - GroupReadsByUmi plot can now be toggled to show relative percents ([#1147](https://github.com/MultiQC/MultiQC/pull/1147))
- **FLASh**
  - Logs not reporting innie and outine uncombined pairs now plot combined pairs instead ([#1173](https://github.com/MultiQC/MultiQC/issues/1173))
- **GATK**
  - Made parsing for VariantEval more tolerant, so that it will work with output from the tool when run in different modes ([#1158](https://github.com/MultiQC/MultiQC/issues/1158))
- **MTNucRatioCalculator**
  - Fixed misleading value suffix in general stats table
- **Picard MarkDuplicates**
  - **Major change** - previously, if multiple libraries (read-groups) were found then only the first would be used and all others ignored. Now, values from all libraries are merged and `PERCENT_DUPLICATION` and `ESTIMATED_LIBRARY_SIZE` are recalculated. Libraries can be kept as separate samples with a new MultiQC configuration option - `picard_config: markdups_merge_multiple_libraries: False`
  - **Major change** - Updated `MarkDuplicates` bar plot to double the read-pair counts, so that the numbers stack correctly. ([#1142](https://github.com/MultiQC/MultiQC/issues/1142))
- **Picard HsMetrics**
  - Updated large table to use columns specified in the MultiQC config. See [docs](https://multiqc.info/docs/#hsmetrics). ([#831](https://github.com/MultiQC/MultiQC/issues/831))
- **Picard WgsMetrics**
  - Updated parsing code to recognise new java class string ([#1114](https://github.com/MultiQC/MultiQC/issues/1114))
- **QualiMap**
  - Fixed QualiMap mean coverage calculation [#1082](https://github.com/MultiQC/MultiQC/issues/1082), [#1077](https://github.com/MultiQC/MultiQC/issues/1082)
- **RSeqC**
  - Support added for output from `geneBodyCoverage2.py` script ([#844](https://github.com/MultiQC/MultiQC/issues/844))
  - Single sample view in the _"Junction saturation"_ plot now works with the toolbox properly _(rename, hide, highlight)_ ([#1133](https://github.com/MultiQC/MultiQC/issues/1133))
- **RNASeQC2**
  - Updated to handle the parsing metric files from the [newer rewrite of RNA-SeqQC](https://github.com/broadinstitute/rnaseqc).
- **Samblaster**
  - Improved parsing to handle variable whitespace ([#1176](https://github.com/MultiQC/MultiQC/issues/1176))
- **Samtools**
  - Removes hardcoding of general stats column names. This allows column names to indicate when a module has been run twice ([https://github.com/MultiQC/MultiQC/issues/1076](https://github.com/MultiQC/MultiQC/issues/1076)).
  - Added an observed over expected read count plot for `idxstats` ([#1118](https://github.com/MultiQC/MultiQC/issues/1118))
  - Added additional (by default hidden) column for `flagstat` that displays number total number of reads in a bam
- **sortmerna**
  - Fix the bug for the latest sortmerna version 4.2.0 ([#1121](https://github.com/MultiQC/MultiQC/issues/1121))
- **sexdeterrmine**
  - Added a scatter plot of relative X- vs Y-coverage to the generated report.
- **VerifyBAMID**
  - Allow files with column header `FREEMIX(alpha)` ([#1112](https://github.com/MultiQC/MultiQC/issues/1112))

#### Bug Fixes

- Added a new test to check that modules work correctly with `--ignore-samples`. A lot of them didn't:
  - `Mosdepth`, `conpair`, `Qualimap BamQC`, `RNA-SeQC`, `GATK BaseRecalibrator`, `SNPsplit`, `SeqyClean`, `Jellyfish`, `hap.py`, `HOMER`, `BBMap`, `DeepTools`, `HiCExplorer`, `pycoQC`, `interop`
  - These modules have now all been fixed and `--ignore-samples` should work as you expect for whatever data you have.
- Removed use of `shutil.copy` to avoid problems with working on multiple filesystems ([#1130](https://github.com/MultiQC/MultiQC/issues/1130))
- Made folder naming behaviour of `multiqc_plots` consistent with `multiqc_data`
  - Incremental numeric suffixes now added if folder already exists
  - Plots folder properly renamed if using `-n`/`--filename`
- Heatmap plotting function is now compatible with MultiQC toolbox `hide` and `highlight` ([#1136](https://github.com/MultiQC/MultiQC/issues/1136))
- Plot config `logswitch_active` now works as advertised
- When running MultiQC modules several times, multiple data files are now created instead of overwriting one another ([#1175](https://github.com/MultiQC/MultiQC/issues/1175))
- Fixed minor bug where tables could report negative numbers of columns in their header text
- Fixed bug where numeric custom content sample names could trigger a `TypeError` ([#1091](https://github.com/MultiQC/MultiQC/issues/1091))
- Fixed custom content bug HTML data in a config file would trigger a `ValueError` ([#1071](https://github.com/MultiQC/MultiQC/issues/1071))
- Replaced deprecated 'warn()' with 'warning()' of the logging module
- Custom content now supports `section_extra` config key to add custom HTML after description.
- Barplots with `ymax` set now ignore this when you click the _Percentages_ tab.

## [MultiQC v1.8](https://github.com/MultiQC/MultiQC/releases/tag/v1.8) - 2019-11-20

#### New Modules

- [**fgbio**](http://fulcrumgenomics.github.io/fgbio/)
  - Process family size count hist data from `GroupReadsByUmi`
- [**biobambam2**](https://github.com/gt1/biobambam2)
  - Added submodule for `bamsormadup` tool
  - Totally cheating - it uses Picard MarkDuplicates but with a custom search pattern and naming
- [**SeqyClean**](https://github.com/ibest/seqyclean)
  - Adds analysis for seqyclean files
- [**mtnucratio**](https://github.com/apeltzer/MTNucRatioCalculator)
  - Added little helper tool to compute mt to nuclear ratios for NGS data.
- [**mosdepth**](https://github.com/brentp/mosdepth)
  - fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing
- [**SexDetErrmine**](https://github.com/TCLamnidis/Sex.DetERRmine)
  - Relative coverage and error rate of X and Y chromosomes
- [**SNPsplit**](https://github.com/FelixKrueger/SNPsplit)
  - Allele-specific alignment sorting

#### Module updates

- **bcl2fastq**
  - Added handling of demultiplexing of more than 2 reads
  - Allow bcl2fastq to parse undetermined barcode information in situations when lane indexes do not start at 1
- **BBMap**
  - Support for scafstats output marked as not yet implemented in docs
- **DeDup**
  - Added handling clusterfactor and JSON logfiles
- **damageprofiler**
  - Added writing metrics to data output file.
- **DeepTools**
  - Fixed Python3 bug with int() conversion ([#1057](https://github.com/MultiQC/MultiQC/issues/1057))
  - Handle varied TES boundary labels in plotProfile ([#1011](https://github.com/MultiQC/MultiQC/issues/1011))
  - Fixed bug that prevented running on only plotProfile files when no other deepTools files found.
- **fastp**
  - Fix faulty column handling for the _after filtering_ Q30 rate ([#936](https://github.com/MultiQC/MultiQC/issues/936))
- **FastQC**
  - When including a FastQC section multiple times in one report, the Per Base Sequence Content heatmaps now behave as you would expect.
  - Added heatmap showing FastQC status checks for every section report across all samples
  - Made sequence content individual plots work after samples have been renamed ([#777](https://github.com/MultiQC/MultiQC/issues/777))
  - Highlighting samples from status - respect chosen highlight colour in the toolbox ([#742](https://github.com/MultiQC/MultiQC/issues/742))
- **FastQ Screen**
  - When including a FastQ Screen section multiple times in one report, the plots now behave as you would expect.
  - Fixed MultiQC linting errors
- **fgbio**
  - Support the new output format of `ErrorRateByReadPosition` first introduced in version `1.3.0`, as well as the old output format.
- **GATK**
  - Refactored BaseRecalibrator code to be more consistent with MultiQC Python style
  - Handle zero count errors in BaseRecalibrator
- **HiC Explorer**
  - Fixed bug where module tries to parse `QC_table.txt`, a new log file in hicexplorer v2.2.
  - Updated the format of the report to fits the changes which have been applied to the QC report of hicexplorer v3.3
  - Updated code to save parsed results to `multiqc_data`
- **HTSeq**
  - Fixed bug where module would crash if a sample had zero reads ([#1006](https://github.com/MultiQC/MultiQC/issues/1006))
- **LongRanger**
  - Added support for the LongRanger Align pipeline.
- **miRTrace**
  - Fixed bug where a sample in some plots was missed. ([#932](https://github.com/MultiQC/MultiQC/issues/932))
- **Peddy**
  - Fixed bug where sample name cleaning could lead to error. ([#1024](https://github.com/MultiQC/MultiQC/issues/1024))
  - All plots (including _Het Check_ and _Sex Check_) now hidden if no data
- **Picard**
  - Modified `OxoGMetrics` so that it will find files created with GATK `CollectMultipleMetrics` and `ConvertSequencingArtifactToOxoG`.
- **QoRTs**
  - Fixed bug where `--dirs` broke certain input files. ([#821](https://github.com/MultiQC/MultiQC/issues/821))
- **Qualimap**
  - Added in mean coverage computation for general statistics report
  - Creates now tables of collected data in `multiqc_data`
- **RNA-SeQC**
  - Updated broken URL link
- **RSeQC**
  - Fixed bug where Junction Saturation plot when clicking a single sample was mislabelling the lines.
  - When including a RSeQC section multiple times in one report, clicking Junction Saturation plot now behaves as you would expect.
  - Fixed bug where exported data in `multiqc_rseqc_read_distribution.txt` files had incorrect values for `_kb` fields ([#1017](https://github.com/MultiQC/MultiQC/issues/1017))
- **Samtools**
  - Utilize in-built `read_count_multiplier` functionality to plot `flagstat` results more nicely
- **SnpEff**
  - Increased the default summary csv file-size limit from 1MB to 5MB.
- **Stacks**
  - Fixed bug where multi-population sum stats are parsed correctly ([#906](https://github.com/MultiQC/MultiQC/issues/906))
- **TopHat**
  - Fixed bug where TopHat would try to run with files from Bowtie2 or HiSAT2 and crash
- **VCFTools**
  - Fixed a bug where `tstv_by_qual.py` produced invalid json from infinity-values.
- **snpEff**
  - Added plot of effects

#### New MultiQC Features

- Added some installation docs for windows
- Added some docs about using MultiQC in bioinformatics pipelines
- Rewrote Docker image
  - New base image `czentye/matplotlib-minimal` reduces image size from ~200MB to ~80MB
  - Proper installation method ensures latest version of the code
  - New entrypoint allows easier command-line usage
- Support opening MultiQC on websites with CSP `script-src 'self'` with some sha256 exceptions
  - Plot data is no longer intertwined with javascript code so hashes stay the same
- Made `config.report_section_order` work for module sub-sections as well as just modules.
- New config options `exclude_modules` and `run_modules` to complement `-e` and `-m` cli flags.
- Command line output is now coloured by default :rainbow: (use `--no-ansi` to turn this off)
- Better launch comparability due to code refactoring by [@KerstenBreuer](https://github.com/KerstenBreuer) and [@ewels](https://github.com/ewels)
  - Windows support for base `multiqc` command
  - Support for running as a python module: `python -m multiqc .`
  - Support for running within a script: `import multiqc` and `multiqc.run('/path/to/files')`
- Config option `custom_plot_config` now works for bargraph category configs as well ([#1044](https://github.com/MultiQC/MultiQC/issues/1044))
- Config `table_columns_visible` can now be given a module namespace and it will hide all columns from that module ([#541](https://github.com/MultiQC/MultiQC/issues/541))

#### Bug Fixes

- MultiQC now ignores all `.md5` files
- Use `SafeLoader` for PyYaml load calls, avoiding recent warning messages.
- Hide `multiqc_config_example.yaml` in the `test` directory to stop people from using it without modification.
- Fixed matplotlib background colour issue (@epakarin - [#886](https://github.com/MultiQC/MultiQC/issues))
- Table rows that are empty due to hidden columns are now properly hidden on page load ([#835](https://github.com/MultiQC/MultiQC/issues/835))
- Sample name cleaning: All sample names are now truncated to their basename, without a path.
  - This includes for `regex` and `replace` (before was only the default `truncate`).
  - Only affects modules that take sample names from file contents, such as cutadapt.
  - See [#897](https://github.com/MultiQC/MultiQC/issues/897) for discussion.

## [MultiQC v1.7](https://github.com/MultiQC/MultiQC/releases/tag/v1.7) - 2018-12-21

#### New Modules

- [**BISCUIT**](https://github.com/zwdzwd/biscuit)
  - BISuilfite-seq CUI Toolkit
  - Module written by [@zwdzwd](https://github.com/zwdzwd/)
- [**DamageProfiler**](https://github.com/Integrative-Transcriptomics/DamageProfiler)
  - A tool to determine ancient DNA misincorporation rates.
  - Module written by [@apeltzer](https://github.com/apeltzer/)
- [**FLASh**](https://ccb.jhu.edu/software/FLASH/)
  - FLASH (Fast Length Adjustment of SHort reads)
  - Module written by [@pooranis](https://github.com/pooranis/)
- [**MinIONQC**](https://github.com/roblanf/minion_qc)
  - QC of reads from ONT long-read sequencing
  - Module written by [@ManavalanG](https://github.com/ManavalanG)
- [**phantompeakqualtools**](https://www.encodeproject.org/software/phantompeakqualtools)
  - A tool for informative enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data.
  - Module written by [@chuan-wang](https://github.com/chuan-wang/)
- [**Stacks**](http://catchenlab.life.illinois.edu/stacks/)
  - A software for analyzing restriction enzyme-based data (e.g. RAD-seq). Support for Stacks >= 2.1 only.
  - Module written by [@remiolsen](https://github.com/remiolsen/)

#### Module updates

- **AdapterRemoval**
  - Handle error when zero bases are trimmed. See [#838](https://github.com/MultiQC/MultiQC/issues/838).
- **Bcl2fastq**
  - New plot showing the top twenty of undetermined barcodes by lane.
  - Informations for R1/R2 are now separated in the General Statistics table.
  - SampleID is concatenate with SampleName because in Chromium experiments several sample have the same SampleName.
- **deepTools**
  - New PCA plots from the `plotPCA` function (written by [@chuan-wang](https://github.com/chuan-wang/))
  - New fragment size distribution plots from `bamPEFragmentSize --outRawFragmentLengths` (written by [@chuan-wang](https://github.com/chuan-wang/))
  - New correlation heatmaps from the `plotCorrelation` function (written by [@chuan-wang](https://github.com/chuan-wang/))
  - New sequence distribution profiles around genes, from the `plotProfile` function (written by [@chuan-wang](https://github.com/chuan-wang/))
  - Reordered sections
- **Fastp**
  - Fixed bug in parsing of empty histogram data. See [#845](https://github.com/MultiQC/MultiQC/issues/845).
- **FastQC**
  - Refactored _Per Base Sequence Content_ plots to show original underlying data, instead of calculating it from the page contents. Now shows original FastQC base-ranges and fixes 100% GC bug in final few pixels. See [#812](https://github.com/MultiQC/MultiQC/issues/812).
  - When including a FastQC section multiple times in one report, the summary progress bars now behave as you would expect.
- **FastQ Screen**
  - Don't hide genomes in the simple plot, even if they have zero unique hits. See [#829](https://github.com/MultiQC/MultiQC/issues/829).
- **InterOp**
  - Fixed bug where read counts and base pair yields were not displaying in tables correctly.
  - Number formatting for these fields can now be customised in the same way as with other modules, as described [in the docs](http://multiqc.info/docs/#number-base-multiplier)
- **Picard**
  - InsertSizeMetrics: You can now configure to what degree the insert size plot should be smoothed.
  - CollectRnaSeqMetrics: Add warning about missing rRNA annotation.
  - CollectRnaSeqMetrics: Add chart for counts/percentage of reads mapped to the correct strand.
  - Now parses VariantCallingMetrics reports. (Similar to GATK module's VariantEval.)
- **phantompeakqualtools**
  - Properly clean sample names
- **Trimmomatic**
  - Updated Trimmomatic module documentation to be more helpful
  - New option to use filenames instead of relying on the command line used. See [#864](https://github.com/MultiQC/MultiQC/issues/864).

#### New MultiQC Features

- Embed your custom images with a new Custom Content feature! Just add `_mqc` to the end of the filename for `.png`, `.jpg` or `.jpeg` files.
- Documentation for Custom Content reordered to make it a little more sane
- You can now add or override any config parameter for any MultiQC plot! See [the documentation](http://multiqc.info/docs/#customising-plots) for more info.
- Allow `table_columns_placement` config to work with table IDs as well as column namespaces. See [#841](https://github.com/MultiQC/MultiQC/issues/841).
- Improved visual spacing between grouped bar plots

#### Bug Fixes

- Custom content no longer clobbers `col1_header` table configs
- The option `--file-list` that refers to a text file with file paths to analyse will no longer ignore directory paths
- [Sample name directory prefixes](https://multiqc.info/docs/#sample-names-prefixed-with-directories) are now added _after_ cleanup.
- If a module is run multiple times in one report, it's CSS and JS files will only be included once (`default` template)

## [MultiQC v1.6](https://github.com/MultiQC/MultiQC/releases/tag/v1.6) - 2018-08-04

Some of these updates are thanks to the efforts of people who attended the [NASPM](https://twitter.com/NordicGenomics) 2018 MultiQC hackathon session. Thanks to everyone who attended!

#### New Modules

- [**fastp**](https://github.com/OpenGene/fastp)
  - An ultra-fast all-in-one FASTQ preprocessor (QC, adapters, trimming, filtering, splitting...)
  - Module started by [@florianduclot](https://github.com/florianduclot/) and completed by [@ewels](https://github.com/ewels/)
- [**hap.py**](https://github.com/Illumina/hap.py)
  - Hap.py is a set of programs based on htslib to benchmark variant calls against gold standard truth datasets
  - Module written by [@tsnowlan](https://github.com/tsnowlan/)
- [**Long Ranger**](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger)
  - Works with data from the 10X Genomics Chromium. Performs sample demultiplexing, barcode processing, alignment, quality control, variant calling, phasing, and structural variant calling.
  - Module written by [@remiolsen](https://github.com/remiolsen/)
- [**miRTrace**](https://github.com/friedlanderlab/mirtrace)
  - A quality control software for small RNA sequencing data.
  - Module written by [@chuan-wang](https://github.com/chuan-wang/)

#### Module updates

- **BCFtools**
  - New plot showing SNP statistics versus quality of call from bcftools stats ([@MaxUlysse](https://github.com/MaxUlysse) and [@Rotholandus](https://github.com/Rotholandus))
- **BBMap**
  - Support added for BBDuk kmer-based adapter/contaminant filtering summary stats ([@boulund](https://github.com/boulund)
- **FastQC**
  - New read count plot, split into unique and duplicate reads if possible.
  - Help text added for all sections, mostly copied from the excellent FastQC help.
  - Sequence duplication plot rescaled
- **FastQ Screen**
  - Samples in large-sample-number plot are now sorted alphabetically ([@hassanfa](https://github.com/hassanfa)
- **MACS2**
  - Output is now more tolerant of missing data (no plot if no data)
- **Peddy**
  - Background samples now shown in ancestry PCA plot ([@roryk](https://github.com/roryk))
  - New plot showing sex checks versus het ratios, supporting unknowns ([@oyvinev](https://github.com/oyvinev))
- **Picard**
  - New submodule to handle `ValidateSamFile` reports ([@cpavanrun](https://github.com/cpavanrun))
  - WGSMetrics now add the mean and standard-deviation coverage to the general stats table (hidden) ([@cpavanrun](https://github.com/cpavanrun))
- **Preseq**
  - New config option to plot preseq plots with unique old coverage on the y axis instead of read count
  - Code refactoring by [@vladsaveliev](https://github.com/vladsaveliev)
- **QUAST**
  - Null values (`-`) in reports now handled properly. Bargraphs always shown despite varying thresholds. ([@vladsaveliev](https://github.com/vladsaveliev))
- **RNA-SeQC**
  - Don't create the report section for Gene Body Coverage if no data is given
- **Samtools**
  - Fixed edge case bug where MultiQC could crash if a sample had zero count coverage with idxstats.
  - Adds % proper pairs to general stats table
- **Skewer**
  - Read length plot rescaled
- **Tophat**
  - Fixed bug where some samples could be given a blank sample name ([@lparsons](https://github.com/lparsons))
- **VerifyBamID**
  - Change column header help text for contamination to match percentage output ([@chapmanb](https://github.com/chapmanb))

#### New MultiQC Features

- New config option `remove_sections` to skip specific report sections from modules
- Add `path_filters_exclude` to exclude certain files when running modules multiple times. You could previously only include certain files.
- New `exclude_*` keys for file search patterns
  - Have a subset of patterns to exclude otherwise detected files with, by filename or contents
- Command line options all now use mid-word hyphens (not a mix of hyphens and underscores)
  - Old underscore terms still maintained for backwards compatibility
- Flag `--view-tags` now works without requiring an "analysis directory".
- Removed Python dependency for `enum34` ([@boulund](https://github.com/boulund))
- Columns can be added to `General Stats` table for custom content/module.
- New `--ignore-symlinks` flag which will ignore symlinked directories and files.
- New `--no-megaqc-upload` flag which disables automatically uploading data to MegaQC

#### Bug Fixes

- Fix path filters for `top_modules/module_order` configuration only selecting if _all_ globs match. It now filters searches that match _any_ glob.
- Empty sample names from cleaning are now no longer allowed
- Stop prepend_dirs set in the config from getting clobbered by an unpassed CLI option ([@tsnowlan](https://github.com/tsnowlan))
- Modules running multiple times now have multiple sets of columns in the General Statistics table again, instead of overwriting one another.
- Prevent tables from clobbering sorted row orders.
- Fix linegraph and scatter plots data conversion (sporadically the incorrect `ymax` was used to drop data points) ([@cpavanrun](https://github.com/cpavanrun))
- Adjusted behavior of ceiling and floor axis limits
- Adjusted multiple file search patterns to make them more specific
  - Prevents the wrong module from accidentally slurping up output from a different tool. By [@cpavanrun](https://github.com/cpavanrun) (see [PR #727](https://github.com/MultiQC/MultiQC/pull/727))
- Fixed broken report bar plots when `-p`/`--export-plots` was specified (see issue [#801](https://github.com/MultiQC/MultiQC/issues/801))

## [MultiQC v1.5](https://github.com/MultiQC/MultiQC/releases/tag/v1.5) - 2018-03-15

#### New Modules

- [**HiCPro**](https://github.com/nservant/HiC-Pro) - New module!
  - HiCPro: Quality controls and processing of Hi-C
  - Module written by [@nservant](https://github.com/nservant),
- [**DeDup**](http://www.github.com/apeltzer/DeDup) - New module!
  - DeDup: Improved Duplicate Removal for merged/collapsed reads in ancient DNA analysis
  - Module written by [@apeltzer](https://github.com/apeltzer),
- [**Clip&Merge**](http://github.com/apeltzer/ClipAndMerge) - New module!
  - Clip&Merge: Adapter clipping and read merging for ancient DNA analysis
  - Module written by [@apeltzer](https://github.com/apeltzer),

#### Module updates

- **bcl2fastq**
  - Catch `ZeroDivisionError` exceptions when there are 0 reads ([@aledj2](https://github.com/aledj2))
  - Add parsing of `TrimmedBases` and new General Stats column for % bases trimmed ([@matthdsm](https://github.com/matthdsm)).
- **BUSCO**
  - Fixed configuration bug that made all sample names become `'short'`
- **Custom Content**
  - Parsed tables now exported to `multiqc_data` files
- **Cutadapt**
  - Refactor parsing code to collect all length trimming plots
- **FastQC**
  - Fixed starting y-axis label for GC-content lineplot being incorrect.
- **HiCExplorer**
  - Updated to work with v2.0 release.
- **Homer**
  - Made parsing of `tagInfo.txt` file more resilient to variations in file format so that it works with new versions of Homer.
  - Kept order of chromosomes in coverage plot consistent.
- **Peddy**
  - Switch `Sex error` logic to `Correct sex` for better highlighting ([@aledj2](https://github.com/aledj2))
- **Picard**
  - Updated module and search patterns to recognise new output format from Picard version >= 2.16 and GATK output.
- **Qualimap BamQC**
  - Fixed bug where start of _Genome Fraction_ could have a step if target is 100% covered.
- **RNA-SeQC**
  - Added rRNA alignment stats to summary table [@Rolandde](https://github.com/Rolandde)
- **RSeqC**
  - Fixed read distribution plot by adding category for `other_intergenic` (thanks to [@moxgreen](https://github.com/moxgreen))
  - Fixed a dodgy plot title (Read GC content)
- **Supernova**
  - Added support for Supernova 2.0 reports. Fixed a TypeError bug when using txt reports only. Also a bug when parsing empty histogram files.

#### New MultiQC Features

- Invalid choices for `--module` or `--exclude` now list the available modules alphabetically.
- Linting now checks for presence in `config.module_order` and tags.

#### Bug Fixes

- Excluding modules now works in combination with using module tags.
- Fixed edge-case bug where certain combinations of `output_fn_name` and `data_dir_name` could trigger a crash
- Conditional formatting - values are now longer double-labelled
- Made config option `extra_series` work in scatter plots the same way that it works for line plots
- Locked the `matplotlib` version to `v2.1.0` and below
  - Due to [two](https://github.com/matplotlib/matplotlib/issues/10476) [bugs](https://github.com/matplotlib/matplotlib/issues/10784) that appeared in `v2.2.0` - will remove this constraint when there's a new release that works again.

## [MultiQC v1.4](https://github.com/MultiQC/MultiQC/releases/tag/v1.4) - 2018-01-11

A slightly earlier-than-expected release due to a new problem with dependency packages that is breaking MultiQC installations since 2018-01-11.

#### New Modules

- [**Sargasso**](http://statbio.github.io/Sargasso/)
  - Parses output from Sargasso - a tool to separate mixed-species RNA-seq reads according to their species of origin
  - Module written by [@hxin](https://github.com/hxin/)
- [**VerifyBAMID**](https://genome.sph.umich.edu/wiki/VerifyBamID)
  - Parses output from VerifyBAMID - a tool to detect contamination in BAM files.
  - Adds the `CHIPMIX` and `FREEMIX` columns to the general statistics table.
  - Module written by [@aledj2](https://github.com/aledj2/)

#### Module updates

- **MACS2**
  - Updated to work with output from older versions of MACS2 by [@avilella](https://github.com/avilella/)
- **Peddy**
  - Add het check plot to suggest potential contamination by [@aledj2](https://github.com/aledj2)
- **Picard**
  - Picard HsMetrics `HS_PENALTY` plot now has correct axis labels
  - InsertSizeMetrics switches commas for points if it can't convert floats. Should help some european users.
- **QoRTs**
  - Added support for new style of output generated in the v1.3.0 release
- **Qualimap**
  - New `Error rate` column in General Statistics table, added by [@Cashalow](https://github.com/Cashalow/)
    - Hidden by default - customise your MultiQC config to always show this column (see [docs](http://multiqc.info/docs/#hiding-columns))
- **QUAST**
  - New option to customise the default display of contig count and length (eg. `bp` instead of `Mbp`).
  - See [documentation](http://multiqc.info/docs/#quast). Written by [@ewels](https://github.com/ewels/) and [@Cashalow](https://github.com/Cashalow/)
- **RSeQC**
  - Removed normalisation in Junction Saturation plot. Now raw counts instead of % of total junctions.

#### New MultiQC Features

- Conditional formatting / highlighting of cell contents in tables
  - If you want to make values that match a criteria stand out more, you can now write custom rules and formatting instructions for tables.
  - For instructions, see [the documentation](http://multiqc.info/docs/#conditional-formatting)
- New `--lint` option which is strict about best-practices for writing new modules
  - Useful when writing new modules and code as it throws warnings
  - Currently only implemented for bar plots and a few other places. More linting coming soon...
- If MultiQC breaks and shows am error message, it now reports the filename of the last log it found
  - Hopefully this will help with debugging / finding dodgy input data

#### Bug Fixes

- Addressed new dependency error with conflicting package requirements
  - There was a conflict between the `networkx`, `colormath` and `spectra` releases.
  - I previously forced certain software versions to get around this, but `spectra` has now updated with the unfortunate effect of introducing a new dependency clash that halts installation.
- Fixed newly introduced bug where Custom Content MultiQC config file search patterns had been broken
- Updated pandoc command used in `--pdf` to work with new releases of Pandoc
- Made config `table_columns_visible` module name key matching case insensitive to make less frustrating

## [MultiQC v1.3](https://github.com/MultiQC/MultiQC/releases/tag/v1.3) - 2017-11-03

#### Breaking changes - custom search patterns

Only for users with custom search patterns for the `bowtie` or `star`: you will
need to update your config files - the `bowtie` search key is now `bowtie1`,
`star_genecounts` is now `star/genecounts`.

For users with custom modules - search patterns _must_ now conform to the search
pattern naming convention: `modulename` or `modulename/anything` (the search pattern
string beginning with the name of your module, anything you like after the first `/`).

#### New Modules

- [**10X Supernova**](https://support.10xgenomics.com/de-novo-assembly/software/overview/welcome)
  - Parses statistics from the _de-novo_ Supernova software.
  - Module written by [@remiolsen](https://github.com/remiolsen/)
- [**BBMap**](https://sourceforge.net/projects/bbmap/)
  - Plot metrics from a number of BBMap tools, a suite of DNA/RNA mapping tools and utilities
  - Module written by [@boulund](https://github.com/boulund/) and [@epruesse](https://github.com/epruesse/)
- [**deepTools**](https://github.com/fidelram/deepTools) - new module!
  - Parse text output from `bamPEFragmentSize`, `estimateReadFiltering`, `plotCoverage`, `plotEnrichment`, and `plotFingerprint`
  - Module written by [@dpryan79](https://github.com/dpryan79/)
- [**Homer Tag Directory**](http://homer.ucsd.edu/homer/ngs/tagDir.html) - new submodule!
  - Module written by [@rdali](https://github.com/rdali/)
- [**illumina InterOp**](http://illumina.github.io/interop/index.html)
  - Module to parse metrics from illumina sequencing runs and demultiplexing, generated by the InterOp package
  - Module written by [@matthdsm](https://github.com/matthdsm/)
- [**RSEM**](https://deweylab.github.io/RSEM/) - new module!
  - Parse `.cnt` file comming from rsem-calculate-expression and plot read repartitions (Unalignable, Unique, Multi ...)
  - Module written by [@noirot](https://github.com/noirot/)
- [**HiCExplorer**](https://github.com/maxplanck-ie/HiCExplorer)
  - New module to parse the log files of `hicBuildMatrix`.
  - Module written by [@joachimwolff](https://github.com/joachimwolff/)

#### Module updates

- **AfterQC**
  - Handle new output format where JSON summary key changed names.
- **bcl2fastq**
  - Clusters per sample plot now has tab where counts are categoried by lane.
- **GATK**
  - New submodule to handle Base Recalibrator stats, written by [@winni2k](https://github.com/winni2k/)
- **HiSAT2**
  - Fixed bug where plot title was incorrect if both SE and PE bargraphs were in one report
- **Picard HsMetrics**
  - Parsing code can now handle commas for decimal places
- **Preseq**
  - Updated odd file-search pattern that limited input files to 500kb
- **QoRTs**
  - Added new plots, new helptext and updated the module to produce a lot more output.
- **Qualimap BamQC**
  - Fixed edge-case bug where the refactored coverage plot code could raise an error from the `range` call.
- Documentation and link fixes for Slamdunk, GATK, bcl2fastq, Adapter Removal, FastQC and main docs
  - Many of these spotted and fixed by [@juliangehring](https://github.com/juliangehring/)
- Went through all modules and standardised plot titles
  - All plots should now have a title with the format _Module name: Plot name_

#### New MultiQC Features

- New MultiQC docker image
  - Ready to use docker image now available at <https://hub.docker.com/r/MultiQC/MultiQC/> (200 MB)
  - Uses automated builds - pull `:latest` to get the development version, future releases will have stable tags.
  - Written by [@MaxUlysse](https://github.com/MaxUlysse/)
- New `module_order` config options allow modules to be run multiple times
  - Filters mean that a module can be run twice with different sets of files (eg. before and after trimming)
  - Custom module config parameters can be passed to module for each run
- File search refactored to only search for running modules
  - Makes search much faster when running with lots of files and limited modules
  - For example, if using `-m star` to only use the STAR module, all other file searches now skipped
- File search now warns if an unrecognised search type is given
- MultiQC now saves nearly all parsed data to a structured output file by default
  - See `multiqc_data/multiqc_data.json`
  - This can be turned off by setting `config.data_dump_file: false`
- Verbose logging when no log files found standardised. Less duplication in code and logs easier to read!
- New documentation section describing how to use MultiQC with Galaxy
- Using `shared_key: 'read_counts'` in table header configs now applies relevant defaults

#### Bug Fixes

- Installation problem caused by changes in upstream dependencies solved by stricter installation requirements
- Minor `default_dev` directory creation bug squashed
- Don't prepend the directory separator (`|`) to sample names with `-d` when there are no subdirs
- `yPlotLines` now works even if you don't set `width`

## [MultiQC v1.2](https://github.com/MultiQC/MultiQC/releases/tag/v1.2) - 2017-08-16

#### CodeFest 2017 Contributions

We had a fantastic group effort on MultiQC at the [2017 BOSC CodeFest](https://www.open-bio.org/wiki/Codefest_2017).
Many thanks to those involved!

#### New Modules

- [**AfterQC**](https://github.com/OpenGene/AfterQC) - New module!
  - Added parsing of the _AfterQC_ json file data, with a plot of filtered reads.
  - Work by [@raonyguimaraes](https://github.com/raonyguimaraes)
- [**bcl2fastq**](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
  - bcl2fastq can be used to both demultiplex data and convert BCL files to FASTQ file formats for downstream analysis
  - New module parses JSON output from recent versions and summarises some key statistics from the demultiplexing process.
  - Work by [@iimog](https://github.com/iimog) (with a little help from [@tbooth](https://github.com/tbooth) and [@ewels](https://github.com/ewels))
- [**leeHom**](https://github.com/grenaud/leeHom)
  - leeHom is a program for the Bayesian reconstruction of ancient DNA
- [**VCFTools**](https://vcftools.github.io)
  - Added initial support for VCFTools `relatedness2`
  - Added support for VCFTools `TsTv-by-count` `TsTv-by-qual` `TsTv-summary`
  - Module written by [@mwhamgenomics](https://github.com/mwhamgenomics)

#### Module updates

- **FastQ Screen**
  - Gracefully handle missing data from very old FastQ Screen versions.
- **RNA-SeQC**
  - Add new transcript-associated reads plot.
- **Picard**
  - New submodule to handle output from `TargetedPcrMetrics`
- **Prokka**
  - Added parsing of the `# CRISPR arrays` data from Prokka when available ([@asetGem](https://github.com/asetGem))
- **Qualimap**
  - Some code refactoring to radically improve performance and run times, especially with high coverage datasets.
  - Fixed bug where _Cumulative coverage genome fraction_ plot could be truncated.

#### New MultiQC Features

- New module help text
  - Lots of additional help text was written to make MultiQC report plots easier to interpret.
  - Updated modules:
    - Bowtie
    - Bowtie 2
    - Prokka
    - Qualimap
    - SnpEff
  - Elite team of help-writers:
    - [@tabwalsh](https://github.com/tabwalsh)
    - [@ddesvillechabrol](https://github.com/tabwalsh)
    - [@asetGem](https://github.com/asetGem)
- New config option `section_comments` allows you to add custom comments above specific sections in the report
- New `--tags` and `--view_tags` command line options
  - Modules can now be given tags (keywords) and filtered by those. So running `--tags RNA` will only run MultiQC modules related to RNA analysis.
  - Work by [@Hammarn](https://github.com/Hammarn)
- Back-end configuration options to specify the order of table columns
  - Modules and user configs can set priorities for columns to customise where they are displayed
  - Work by [@tbooth](https://github.com/tbooth)
- Added framework for proper unit testing
  - Previous start on unit tests tidied up, new blank template and tests for the `clean_sample_name` functionality.
  - Added to Travis and Appveyor for continuous integration testing.
  - Work by [@tbooth](https://github.com/tbooth)
- Bug fixes and refactoring of report configuration saving / loading
  - Discovered and fixed a bug where a report config could only be loaded once
  - Work by [@DennisSchwartz](https://github.com/DennisSchwartz)
- Table column row headers (sample names) can now be numeric-only.
  - Work by [@iimog](https://github.com/iimog)
- Improved sample name cleaning functionality
  - Added option `regex_keep` to clean filenames by _keeping_ the matching part of a pattern
  - Work by [@robinandeer](https://github.com/robinandeer)
- Handle error when invalid regexes are given in reports
  - Now have a nice toast error warning you and the invalid regexes are highlighted
  - Previously this just crashed the whole report without any warning
  - Work by [@robinandeer](https://github.com/robinandeer)
- Command line option `--dirs-depth` now sets `-d` to `True` (so now works even if `-d` isn't also specified).
- New config option `config.data_dump_file` to export as much data as possible to `multiqc_data/multiqc_data.json`
- New code to send exported JSON data to a a web server
  - This is in preparation for the upcoming MegaQC project. Stay tuned!

#### Bug Fixes

- Specifying multiple config files with `-c`/`--config` now works as expected
  - Previously this would only read the last specified
- Fixed table rendering bug that affected Chrome v60 and IE7-11
  - Table cell background bars weren't showing up. Updated CSS to get around this rendering error.
- HTML ID cleanup now properly cleans strings so that they work with jQuery as expected.
- Made bar graph sample highlighting work properly again
- Config `custom_logo` paths can now be relative to the config file (or absolute as before)
- Report doesn't keep annoyingly telling you that toolbox changes haven't been applied
  - Now uses more subtle _toasts_ and only when you close the toolbox (not every click).
- Switching report toolbox options to regex mode now enables the _Apply_ button as it should.
- Sorting table columns with certain suffixes (eg. `13X`) no works properly (numerically)
- Fixed minor bug in line plot data smoothing (now works with unsorted keys)

---

## [MultiQC v1.1](https://github.com/MultiQC/MultiQC/releases/tag/v1.1) - 2017-07-18

#### New Modules

- [**BioBloom Tools**](https://github.com/bcgsc/biobloom)
  - Create Bloom filters for a given reference and then to categorize sequences
- [**Conpair**](https://github.com/nygenome/Conpair)
  - Concordance and contamination estimator for tumor–normal pairs
- [**Disambiguate**](https://github.com/AstraZeneca-NGS/disambiguate)
  - Bargraph displaying the percentage of reads aligning to two different reference genomes.
- [**Flexbar**](https://github.com/seqan/flexbar)
  - Flexbar is a tool for flexible barcode and adapter removal.
- [**HISAT2**](https://ccb.jhu.edu/software/hisat2/)
  - New module for the HISAT2 aligner.
  - Made possible by updates to HISAT2 logging by @infphilo (requires `--new-summary` HISAT2 flag).
- [**HOMER**](http://homer.ucsd.edu/homer/)
  - Support for summary statistics from the `findPeaks` tool.
- [**Jellyfish**](http://www.cbcb.umd.edu/software/jellyfish/)
  - Histograms to estimate library complexity and coverage from k-mer content.
  - Module written by @vezzi
- [**MACS2**](https://github.com/taoliu/MACS)
  - Summary of redundant rate from MACS2 peak calling.
- [**QoRTs**](http://hartleys.github.io/QoRTs/)
  - QoRTs is toolkit for analysis, QC and data management of RNA-Seq datasets.
- [**THetA2**](http://compbio.cs.brown.edu/projects/theta/)
  - THeTA2 _(Tumor Heterogeneity Analysis)_ estimates tumour purity and clonal / subclonal copy number.

#### Module updates

- **BCFtools**
  - Option to collapse complementary changes in substitutions plot, useful for non-strand specific experiments (thanks to @vladsaveliev)
- **Bismark**
  - M-Bias plots no longer show read 2 for single-end data.
- **Custom Content**
  - New option to print raw HTML content to the report.
- **FastQ Screen**
  - Fixed edge-case bug where many-sample plot broke if total number of reads was less than the subsample number.
  - Fixed incorrect logic of config option `fastqscreen_simpleplot` (thanks to @daler)
  - Organisms now alphabetically sorted in fancy plot so that order is nonrandom (thanks to @daler)
  - Fixed bug where `%No Hits` was missed in logs from recent versions of FastQ Screen.
- **HTSeq Counts**
  - Fixed but so that module still works when `--additional-attr` is specified in v0.8 HTSeq above (thanks to @nalcala)
- **Picard**
  - CollectInsertSize: Fixed bug that could make the General Statistics _Median Insert Size_ value incorrect.
  - Fixed error in sample name regex that left trailing `]` characters and was generally broken (thanks to @jyh1 for spotting this)
- **Preseq**
  - Improved plots display (thanks to @vladsaveliev)
- **Qualimap**
  - Only calculate bases over target coverage for values in General Statistics. Should give a speed increase for very high coverage datasets.
- **QUAST**
  - Module is now compatible with runs from [MetaQUAST](http://quast.sourceforge.net/metaquast) (thanks to @vladsaveliev)
- **RSeQC**
  - Changed default order of sections
  - Added config option to reorder and hide module report sections

#### New MultiQC features

- If a report already exists, execution is no longer halted.
  - `_1` is appended to the filename, iterating if this also exists.
  - `-f`/`--force` still overwrites existing reports as before
  - Feature written by [@Hammarn](https://github.com/Hammarn)
- New ability to run modules multiple times in a single report
  - Each run can be given different configuration options, including filters for input files
  - For example, have FastQC after trimming as well as FastQC before trimming.
  - See the relevant [documentation](http://multiqc.info/docs/#order-of-modules) for more instructions.
- New option to customise the order of report _sections_
  - This is in addition / alternative to changing the order of module execution
  - Allows one module to have sections in multiple places (eg. Custom Content)
- Tables have new column options `floor`, `ceiling` and `minRange`.
- Reports show warning if JavaScript is disabled
- Config option `custom_logo` now works with file paths relative to config file directory and cwd.

#### Bug Fixes

- Table headers now sort columns again after scrolling the table
- Fixed buggy table header tooltips
- Base `clean_s_name` function now strips excess whitespace.
- Line graphs don't smooth lines if not needed (number of points < maximum number allowed)
- PDF output now respects custom output directory.

---

## [MultiQC v1.0](https://github.com/MultiQC/MultiQC/releases/tag/v1.0) - 2017-05-17

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
    fn: "fastqc_data.txt"
  zip:
    fn: "*_fastqc.zip"
```

to this:

```yaml
fastqc/data:
  fn: "fastqc_data.txt"
fastqc/zip:
  fn: "*_fastqc.zip"
```

See the [documentation](http://multiqc.info/docs/#step-1-find-log-files) for instructions on how to write the new file search syntax.

See [`search_patterns.yaml`](multiqc/config/search_patterns.yaml) for the new module search keys
and more examples.

#### 2. You have custom plugins / modules / external code

To see what changes need to applied to your custom plugin code, please see the [MultiQC docs](http://multiqc.info/docs/#v1.0-updates).

#### New Modules

- [**Adapter Removal**](https://github.com/mikkelschubert/adapterremoval)
  - AdapterRemoval v2 - rapid adapter trimming, identification, and read merging
- [**BUSCO**](http://busco.ezlab.org/)
  - New module for the `BUSCO v2` tool, used for assessing genome assembly and annotation completeness.
- [**Cluster Flow**](http://clusterflow.io)
  - Cluster Flow is a workflow tool for bioinformatics pipelines. The new module parses executed tool commands.
- [**RNA-SeQC**](http://archive.broadinstitute.org/cancer/cga/rna-seqc)
  - New module to parse output from RNA-SeQC, a java program which computes a series
    of quality control metrics for RNA-seq data.
- [**goleft indexcov**](https://github.com/brentp/goleft/tree/master/indexcov)
  - [goleft indexcov](https://github.com/brentp/goleft/tree/master/indexcov) uses the PED and ROC
    data files to create diagnostic plots of coverage per sample, helping to identify sample gender and coverage issues.
  - Thanks to @chapmanb and @brentp
- [**SortMeRNA**](http://bioinfo.lifl.fr/RNA/sortmerna/)
  - New module for `SortMeRNA`, commonly used for removing rRNA contamination from datasets.
  - Written by @bschiffthaler

#### Module updates

- **Bcftools**
  - Fixed bug with display of indels when only one sample
- **Cutadapt**
  - Now takes the filename if the sample name is `-` (stdin). Thanks to @tdido
- **FastQC**
  - Data for the Sequence content plot can now be downloaded from reports as a JSON file.
- **FastQ Screen**
  - Rewritten plotting method for high sample numbers plot (~ > 20 samples)
  - Now shows counts for single-species hits and bins all multi-species hits
  - Allows plot to show proper percentage view for each sample, much easier to interpret.
- **HTSeq**
  - Fix bug where header lines caused module to crash
- **Picard**
  - New `RrbsSummaryMetrics` Submodule!
  - New `WgsMetrics` Submodule!
  - `CollectGcBiasMetrics` module now prints summary statistics to `multiqc_data` if found. Thanks to @ahvigil
- **Preseq**
  - Now trims the x axis to the point that meets 90% of `min(unique molecules)`.
    Hopefully prevents ridiculous x axes without sacrificing too much useful information.
  - Allows to show estimated depth of coverage instead of less informative molecule counts
    (see [details](http://multiqc.info/docs/#preseq)).
  - Plots dots with externally calculated real read counts (see [details](http://multiqc.info/docs/#preseq)).
- **Qualimap**
  - RNASeq Transcript Profile now has correct axis units. Thanks to @roryk
  - BamQC module now doesn't crash if reports don't have genome gc distributions
- **RSeQC**
  - Fixed Python3 error in Junction Saturation code
  - Fixed JS error for Junction Saturation that made the single-sample combined plot only show _All Junctions_

#### Core MultiQC updates

- Change in module structure and import statements (see [details](http://multiqc.info/docs/#v1.0-updates)).
- Module file search has been rewritten (see above changes to configs)
  - Significant improvement in search speed (test dataset runs in approximately half the time)
  - More options for modules to find their logs, eg. filename and contents matching regexes (see the [docs](http://multiqc.info/docs/#step-1-find-log-files))
- Report plot data is now compressed, significantly reducing report filesizes.
- New `--ignore-samples` option to skip samples based on parsed sample name
  - Alternative to filtering by input filename, which doesn't always work
  - Also can use config vars `sample_names_ignore` (glob patterns) and `sample_names_ignore_re` (regex patterns).
- New `--sample-names` command line option to give file with alternative sample names
  - Allows one-click batch renaming in reports
- New `--cl_config` option to supply MultiQC config YAML directly on the command line.
- New config option to change numeric multiplier in General Stats
  - For example, if reports have few reads, can show `Thousands of Reads` instead of `Millions of Reads`
  - Set config options `read_count_multiplier`, `read_count_prefix` and `read_count_desc`
- Config options `decimalPoint_format` and `thousandsSep_format` now apply to tables as well as plots
  - By default, thosands will now be separated with a space and `.` used for decimal places.
- Tables now have a maximum-height by default and scroll within this.
  - Speeds up report rendering in the web browser and makes report less stupidly long with lots of samples
  - Button beneath table toggles full length if you want a zoomed-out view
  - Refactored and removed previous code to make the table header "float"
  - Set `config.collapse_tables` to `False` to disable table maximum-heights
- Bar graphs and heatmaps can now be zoomed in on
  - Interactive plots sometimes hide labels due to lack of space. These can now be zoomed in on to see specific samples in more detail.
- Report plots now load sequentially instead of all at once
  - Prevents the browser from locking up when large reports load
- Report plot and section HTML IDs are now sanitised and checked for duplicates
- New template available (called _sections_) which has faster loading
  - Only shows results from one module at a time
  - Makes big reports load in the browser much more quickly, but requires more clicking
  - Try it out by specifying `-t sections`
- Module sections tidied and refactored
  - New helper function `self.add_section()`
  - Sections hidden in nav if no title (no more need for the hacky `self.intro +=`)
  - Content broken into `description`, `help` and `plot`, with automatic formatting
  - Empty module sections are now skipped in reports. No need to check if a plot function returns `None`!
  - Changes should be backwards-compatible
- Report plot data export code refactored
  - Now doesn't export hidden samples (uses HighCharts [export-csv](https://github.com/highcharts/export-csv) plugin)
- Handle error when `git` isn't installed on the system.
- Refactored colouring of table cells
  - Was previously done in the browser using [chroma.js](http://gka.github.io/chroma.js/)
  - Now done at report generation time using the [spectra](https://pypi.python.org/pypi/spectra) package
  - Should helpfully speed up report rendering time in the web browser, especially for large reports
- Docs updates (thanks to @varemo)
- Previously hidden log file `.multiqc.log` renamed to `multiqc.log` in `multiqc_data`
- Added option to load MultiQC config file from a path specified in the environment variable `MULTIQC_CONFIG_PATH`
- New table configuration options
  - `sortRows: False` prevents table rows from being sorted alphabetically
  - `col1_header` allows the default first column header to be changed from "Sample Name"
- Tables no longer show _Configure Columns_ and _Plot_ buttons if they only have a single column
- Custom content updates
  - New `custom_content`/`order` config option to specify order of Custom Content sections
  - Tables now use the header for the first column instead of always having `Sample Name`
  - JSON + YAML tables now remember order of table columns
  - Many minor bugfixes
- Line graphs and scatter graphs axis limits
  - If limits are specified, data exceeding this is no longer saved in report
  - Visually identical, but can make report file sizes considerable smaller in some cases
- Creating multiple plots without a config dict now works (previously just gave grey boxes in report)
- All changes are now tested on a Windows system, using [AppVeyor](https://ci.appveyor.com/project/MultiQC/MultiQC/)
- Fixed rare error where some reports could get empty General Statistics tables when no data present.
- Fixed minor bug where config option `force: true` didn't work. Now you don't have to always specify `-f`!

---

## [MultiQC v0.9](https://github.com/MultiQC/MultiQC/releases/tag/v0.9) - 2016-12-21

A major new feature is released in v0.9 - support for _custom content_. This means
that MultiQC can now easily include output from custom scripts within reports without
the need for a new module or plugin. For more information, please see the
[MultiQC documentation](http://multiqc.info/docs/#custom-content).

#### New Modules

- [**HTSeq**](http://www-huber.embl.de/HTSeq/doc/count.html)
  - New module for the `htseq-count` tool, often used in RNA-seq analysis.
- [**Prokka**](http://www.vicbioinformatics.com/software.prokka.shtml)
  - Prokka is a software tool for the rapid annotation of prokaryotic genomes.
- [**Slamdunk**](http://t-neumann.github.io/slamdunk/)
  - Slamdunk is a software tool to analyze SLAMSeq data.
- [**Peddy**](https://github.com/brentp/peddy)
  - Peddy calculates genotype :: pedigree correspondence checks, ancestry checks and sex checks using VCF files.

#### Module updates

- **Cutadapt**
  - Fixed bug in General Stats table number for old versions of cutadapt (pre v1.7)
  - Added support for _really_ old cutadapt logs (eg. v.1.2)
- **FastQC**
  - New plot showing total overrepresented sequence percentages.
  - New option to parse a file containing a theoretical GC curve to display in the background.
    - Human & Mouse Genome / Transcriptome curves bundled, or make your own using
      [fastqcTheoreticalGC](https://github.com/mikelove/fastqcTheoreticalGC). See the
      [MultiQC docs](http://multiqc.info/docs/#fastqc) for more information.
- **featureCounts**
  - Added parsing checks and catch failures for when non-featureCounts files are picked up by accident
- **GATK**
  - Fixed logger error in VariantEval module.
- **Picard**
  - Fixed missing sample overwriting bug in `RnaSeqMetrics`
  - New feature to customise coverage shown from `HsMetrics` in General Statistics table
    see the [docs](http://multiqc.info/docs/#picard) for info).
  - Fixed compatibility problem with output from `CollectMultipleMetrics` for `CollectAlignmentSummaryMetrics`
- **Preseq**
  - Module now recognises output from `c_curve` mode.
- **RSeQC**
  - Made the gene body coverage plot show the percentage view by default
  - Made gene body coverage properly handle sample names
- **Samtools**
  - New module to show duplicate stats from `rmdup` logs
  - Fixed a couple of niggles in the idxstats plot
- **SnpEff**
  - Fixed swapped axis labels in the Variant Quality plot
- **STAR**
  - Fixed crash when there are 0 unmapped reads.
  - Sample name now taken from the directory name if no file prefix found.
- **Qualimap BamQC**
  - Add a line for pre-calculated reference genome GC content
  - Plot cumulative coverage for values above 50x, align with the coverage histogram.
  - New ability to customise coverage thresholds shown in General Statistics table
    (see the [docs](http://multiqc.info/docs/#qualimap) for info).

#### Core MultiQC updates

- Support for _custom content_ (see top of release notes).
- New ninja report tool: make scatter plots of any two table columns!
- Plot data now saved in `multiqc_data` when 'flat' image plots are created
  - Allows you easily re-plot the data (eg. in Excel) for further downstream investigation
- Added _'Apply'_ button to Highlight / Rename / Hide.
  - These tools can become slow with large reports. This means that you can enter several
    things without having to wait for the report to replot each change.
- Report heatmaps can now be sorted by highlight
- New config options `decimalPoint_format` and `thousandsSep_format`
  - Allows you to change the default `1 234.56` number formatting for plots.
- New config option `top_modules` allows you to specify modules that should come at the top of the report
- Fixed bar plot bug where missing categories could shift data between samples
- Report title now printed in the side navigation
- Missing plot IDs added for easier plot exporting
- Stopped giving warnings about skipping directories (now a debug message)
- Added warnings in report about missing functionality for flat plots (exporting and toolbox)
- Export button has contextual text for images / data
- Fixed a bug where user config files were loaded twice
- Fixed bug where module order was random if `--module` or `--exclude` was used.
- Refactored code so that the order of modules can be changed in the user config
- Beefed up code + docs in scatter plots back end and multiple bar plots.
- Fixed a few back end nasties for Tables
  - Shared-key columns are no longer forced to share colour schemes
  - Fixed bug in lambda modified values when format string breaks
  - Supplying just data with no header information now works as advertised
- Improvements to back end code for bar plots
  - New `tt_decimals` and `tt_suffix` options for bar plots
  - Bar plots now support `yCeiling`, `yFloor` and `yMinRange`, as with line plots.
  - New option `hide_zero_cats:False` to force legends to be shown even when all data is 0
- General Stats _Showing x of y_ columns count is fixed on page load.
- Big code whitespace cleanup

---

## [MultiQC v0.8](https://github.com/MultiQC/MultiQC/releases/tag/v0.8) - 2016-09-26

#### New Modules

- [**GATK**](https://software.broadinstitute.org/gatk/)
  - Added support for VariantEval reports, only parsing a little of the information
    in there so far, but it's a start.
  - Module originally written by @robinandeer at the [OBF Codefest](https://www.open-bio.org/wiki/Codefest_2016),
    finished off by @ewels
- [**Bcftools**](https://samtools.github.io/bcftools/)
- [**QUAST**](http://quast.bioinf.spbau.ru/)
  - QUAST is a tool for assessing de novo assemblies against reference genomes.

#### Module updates

- **Bismark** now supports reports from `bam2nuc`, giving Cytosine coverage in General Stats.
- **Bowtie1**
  - Updated to try to find bowtie command before log, handle multiple logs in one file. Same as bowtie2.
- **FastQC**
  - Sample pass/warn/fail lists now display properly even with large numbers of samples
  - Sequence content heatmap display is better with many samples
- **Kallisto**
  - Now supports logs from SE data.
- **Picard**
  - `BaseDistributionByCycle` - new submodule! Written by @mlusignan
  - `RnaSeqMetrics` - new submodule! This one by @ewels ;)
  - `AlignmentSummaryMetrics` - another new submodule!
  - Fixed truncated files crash bug for Python 3 _(#306)_
- **Qualimap RNASeqQC**
  - Fixed parsing bug affecting counts in _Genomic Origin_ plot.
  - Module now works with European style thousand separators (`1.234,56` instead of `1,234.56`)
- **RSeQC**
  - `infer_experiment` - new submodule! Written by @Hammarn
- **Samtools**
  - `stats` submodule now has separate bar graph showing alignment scores
  - `flagstat` - new submodule! Written by @HLWiencko
  - `idxstats` - new submodule! This one by @ewels again

#### Core MultiQC updates

- New `--export`/`-p` option to generate static images plot in `multiqc_plots` (`.png`, `.svg` and `.pdf`)
  - Configurable with `export_plots`, `plots_dir_name` and `export_plot_formats` config options
  - `--flat` option no longer saves plots in `multiqc_data/multiqc_plots`
- New `--comment`/`-b` flag to add a comment to the top of reports.
- New `--dirs-depth`/`-dd` flag to specify how many directories to prepend with `--dirs`/`-d`
  - Specifying a postive number will take that many directories from the end of the path
  - A negative number will take directories from the start of the path.
- Directory paths now appended before cleaning, so `fn_clean_exts` will now affect these names.
- New `custom_logo` attributes to add your own logo to reports.
- New `report_header_info` config option to add arbitrary information to the top of reports.
- New `--pdf` option to create a PDF report
  - Depends on [Pandoc](http://pandoc.org) being installed and is in a beta-stage currently.
  - Note that specifying this will make MultiQC use the `simple` template, giving a HTML report with
    much reduced functionality.
- New `fn_clean_sample_names` config option to turn off sample name cleaning
  - This will print the full filename for samples. Less pretty reports and rows
    on the General Statistics table won't line up, but can prevent overwriting.
- Table header defaults can now be set easily
- General Statistics table now hidden if empty.
- Some new defaults in the sample name cleaning
- Updated the `simple` template.
  - Now has no toolbox or nav, no JavaScript and is better suited for printing / PDFs.
  - New `config.simple_output` config flag so code knows when we're trying to avoid JS.
- Fixed some bugs with config settings (eg. template) being overwritten.
- NFS log file deletion bug fixed by @brainstorm (#265)
- Fixed bug in `--ignore` behaviour with directory names.
- Fixed nasty bug in beeswarm dot plots where sample names were mixed up (#278)
- Beeswarm header text is now more informative (sample count with more info on a tooltip)
- Beeswarm plots now work when reports have > 1000 samples
- Fixed some buggy behaviour in saving / loading report highlighting + renaming configs (#354)

Many thanks to those at the [OpenBio Codefest 2016](https://www.open-bio.org/wiki/Codefest_2016)
who worked on MultiQC projects.

---

## [MultiQC v0.7](https://github.com/MultiQC/MultiQC/releases/tag/v0.7) - 2016-07-04

#### Module updates

- [**Kallisto**](https://pachterlab.github.io/kallisto/) - new module!
- **Picard**
  - Code refactored to make maintenance and additions easier.
  - Big update to `HsMetrics` parsing - more results shown in report, new plots (by @lpantano)
  - Updated `InsertSizeMetrics` to understand logs generated by `CollectMultipleMetrics` (#215)
  - Newlines in picard output. Fixed by @dakl
- **Samtools**
  - Code refactored
  - Rewrote the `samtools stats` code to display more stats in report with a beeswarm plot.
- **Qualimap**
  - Rewritten to use latest methods and fix bugs.
  - Added _Percentage Aligned_ column to general stats for `BamQC` module.
  - Extra table thresholds added by @avilella (hidden by default)
- **General Statistics**
  - Some tweaks to the display defaults (FastQC, Bismark, Qualimap, SnpEff)
  - Now possible to skip the General Statistics section of the report with `--exclude general_stats`
- **Cutadapt** module updated to recognise logs from old versions of cutadapt (<= v1.6)
- **Trimmomatic**
  - Now handles `,` decimal places in percentage values.
  - Can cope with line breaks in log files (see issue #212)
- **FastQC** refactored
  - Now skips zip files if the sample name has already been found. Speeds up MultiQC execution.
  - Code cleaned up. Parsing and data-structures standardised.
  - New popovers on Pass / Warn / Fail status bars showing sample names. Fast highlighting and hiding.
  - New column in General Stats (hidden by default) showing percentage of FastQC modules that failed.
- **SnpEff**
  - Search pattern now more generic, should match reports from others.
  - _Counts by Effect_ plot removed (had hundreds of categories, was fairly unusable).
  - `KeyError` bug fixed.
- **Samblaster** now gets sample name from `ID` instead of `SM` (@dakl)
- **Bowtie 2**
  - Now parses overall alignment rate as intended.
  - Now depends on even less log contents to work with more inputs.
- **MethylQA** now handles variable spacing in logs
- **featureCounts** now splits columns on tabs instead of whitespace, can handle filenames with spaces

#### Core MultiQC updates

- **Galaxy**: MultiQC now available in Galax! Work by @devengineson / @yvanlebras / @cmonjeau
  - See it in the [Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/view/engineson/multiqc/)
- **Heatmap**: New plot type!
- **Scatter Plot**: New plot type!
- **Download raw data** behind plots in reports! Available in the Export toolbox.
  - Choose from tab-separated, comma-separated and the complete JSON.
- **Table columns can be hidden** on page load (shown through _Configure Columns_)
  - Defaults are configurable using the `table_columns_visible` config option.
- **Beeswarm plot**: Added missing rename / highlight / hiding functionality.
- New `-l` / `--file-list` option: specify a file containing a **list of files** to search.
- **Updated HighCharts** to v4.2.5. Added option to export to JPEG.
- Can now **cancel execution** with a single `ctrl+c` rather than having to button mash
- More granular control of **skipping files** during scan (filename, dirname, path matching)
  - Fixed `--exclude` so that it works with directories as well as files
- **New _Clear_ button** in toolbox to bulk remove highlighting / renaming / hiding filters.
- Improved documentation about behaviour for large sample numbers.
- Handle YAML parsing errors for the config file more gracefully
- Removed empty columns from tables again
- Fixed bug in changing module search patterns, reported by @lweasel
- Added timeout parameter to version check to prevent hang on systems with long defaults
- Fixed table display bug in Firefox
- Fixed bug related to order in which config files are loaded
- Fixed bug that broke the _"Show only"_ toolbox feature with multiple names.
- Numerous other small bugs.

---

## [MultiQC v0.6](https://github.com/MultiQC/MultiQC/releases/tag/v0.6) - 2016-04-29

#### Module updates

- New [Salmon](http://combine-lab.github.io/salmon/) module.
- New [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) module.
- New [Bamtools stats](https://github.com/pezmaster31/bamtools) module.
- New beeswarm plot type. General Stats table replaced with this when many samples in report.
- New RSeQC module: Actually a suite of 8 new modules supporting various outputs from RSeQC
- Rewrote bowtie2 module: Now better at parsing logs and tries to scrape input from wrapper logs.
- Made cutadapt show counts by default instead of obs/exp
- Added percentage view to Picard insert size plot

#### Core MultiQC updates

- Dynamic plots now update their labels properly when changing datasets and to percentages
- Config files now loaded from working directory if present
- Started new docs describing how each module works
- Refactored featureCounts module. Now handles summaries describing multiple samples.
- Stopped using so many hidden files. `.multiqc.log` now called `multiqc.log`
- New `-c`/`--config` command line option to specify a MultiQC configuration file
- Can now load run-specific config files called `multiqc_config.yaml` in working directory
- Large code refactoring - moved plotting code out of `BaseModule` and into new `multiqc.plots` submodules
- Generalised code used to generate the General Stats table so that it can be used by modules
- Removed interactive report tour, replaced with a link to a youtube tutorial
- Made it possible to permanently hide the blue welcome message for all future reports
- New option to smooth data for line plots. Avoids mega-huge plots. Applied to SnpEff, RSeQC, Picard.

Bugfixes:

- Qualimap handles infinity symbol (thanks @chapmanb )
- Made SnpEff less fussy about required fields for making plots
- UTF-8 file paths handled properly in Py2.7+
- Extending two config variables wasn't working. Now fixed.
- Dragging the height bar of plots now works again.
- Plots now properly change y axis limits and labels when changing datasets
- Flat plots now have correct path in `default_dev` template

---

## [MultiQC v0.5](https://github.com/MultiQC/MultiQC/releases/tag/v0.5) - 2016-03-29

#### Module updates

- New [Skewer](https://github.com/relipmoc/skewer) module, written by @dakl
- New [Samblaster](https://github.com/GregoryFaust/samblaster) module, written by @dakl
- New [Samtools stats](http://www.htslib.org/) module, written by @lpantano
- New [HiCUP](http://www.bioinformatics.babraham.ac.uk/projects/hicup/) module
- New [SnpEff](http://snpeff.sourceforge.net/) module
- New [methylQA](http://methylqa.sourceforge.net/) module

#### Core MultiQC updates

- New "Flat" image plots, rendered at run time with MatPlotLib
  - By default, will use image plots if > 50 samples (set in config as `plots_flat_numseries`)
  - Means that _very_ large numbers of samples can be viewed in reports. _eg._ single cell data.
  - Templates can now specify their own plotting functions
  - Use `--flat` and `--interactive` to override this behaviour
- MultiQC added to `bioconda` (with help from @dakl)
- New plugin hook: `config_loaded`
- Plugins can now add new command line options (thanks to @robinandeer)
- Changed default data directory name from `multiqc_report_data` to `multiqc_data`
- Removed support for depreciated MultiQC_OSXApp
- Updated logging so that a verbose `multiqc_data/.multiqc.log` file is always written
- Now logs more stuff in verbose mode - command used, user configs and so on.
- Added a call to multiqc.info to check for new versions. Disable with config `no_version_check`
- Removed general stats manual row sorting.
- Made filename matching use glob unix style filename match patterns
- Everything (including the data directory) is now created in a temporary directory and moved when MultiQC is complete.
- A handful of performance updates for large analysis directories

---

## [MultiQC v0.4](https://github.com/MultiQC/MultiQC/releases/tag/v0.4) - 2016-02-16

- New `multiqc_sources.txt` which identifies the paths used to collect all report data for each sample
- Export parsed data as tab-delimited text, `JSON` or `YAML` using the new `-k`/`--data-format` command line option
- Updated HighCharts from `v4.2.2` to `v4.2.3`, fixes tooltip hover bug.
- Nicer export button. Now tied to the export toolbox, hopefully more intuitive.
- FastQC: Per base sequence content heatmap can now be clicked to show line graph for single sample
- FastQC: No longer show adapter contamination datasets with <= 0.1% contamination.
- Picard: Added support for `CollectOxoGMetrics` reports.
- Changed command line option `--name` to `--filename`
- `--name` also used for filename if `--filename` not specified.
- Hide samples toolbox now has switch to _show only_ matching samples
- New regex help box with examples added to report
- New button to copy general stats table to the clipboard
- General Stats table 'floating' header now sorts properly when scrolling
- Bugfix: MultiQC default_dev template now copies module assets properly
- Bufgix: General Stats table floating header now resizes properly when page width changes

---

## [MultiQC v0.3.2](https://github.com/MultiQC/MultiQC/releases/tag/v0.3.2) - 2016-02-08

- All modules now load their log file search parameters from a config
  file, allowing you to overwrite them using your user config file
  - This is useful if your analysis pipeline renames program outputs
- New Picard (sub)modules - Insert Size, GC Bias & HsMetrics
- New Qualimap (sub)module - RNA-Seq QC
- Made Picard MarkDups show percent by default instead of counts
- Added M-Bias plot to Bismark
- New option to stream report HTML to `stdout`
- Files can now be specified as well as directories
- New options to specify whether the parsed data directory should be created
  - command line flags: `--data` / `--no-data`
  - config option name: `make_data_dir`
- Fixed bug with incorrect path to installation dir config YAML file
- New toolbox drawer for bulk-exporting graph images
- Report side navigation can now be hidden to maximise horizontal space
- Mobile styling improved for narrow screen
- More vibrant colours in the general stats table
- General stats table numbers now left aligned
- Settings now saved and loaded to named localstorage locations
  - Simplified interface - no longer global / single report saving
  - Removed static file config. Solves JS error, no-one was doing this
    since we have standalone reports anyway.
- Added support for Python 3.5
- Fixed bug with module specific CSS / JS includes in some templates
- Made the 'ignore files' config use unix style file pattern matching
- Fixed some bugs in the FastQ Screen module
- Fixed some bugs in the FastQC module
- Fixed occasional general stats table bug
- Table sorting on sample names now works after renaming
- Bismark module restructure
  - Each report type now handled independently (alignment / dedup / meth extraction)
  - M-Bias plot now shows R1 and R2
- FastQC GC content plot now has option for counts or percentages
  - Allows comparison between samples with very different read counts
- Bugfix for reports javascript
  - Caused by updated to remotely loaded HighCharts export script
  - Export script now bundled with multiqc, so does not depend on internet connection
  - Other JS errors fixed in this work
- Bugfix for older FastQC reports - handle old style sequence dup data
- Bugfix for varying Tophat alignment report formats
- Bugfix for Qualimap RNA Seq reports with paired end data

---

## [MultiQC v0.3.1](https://github.com/MultiQC/MultiQC/releases/tag/v0.3.1) - 2015-11-04

- Hotfix patch to fix broken FastQC module (wasn't finding `.zip` files properly)
- General Stats table colours now flat. Should improve browser speed.
- Empty rows now hidden if appear due to column removal in general stats
- FastQC Kmer plot removed until we have something better to show.

---

## [MultiQC v0.3](https://github.com/MultiQC/MultiQC/releases/tag/v0.3) - 2015-11-04

- Lots of lovely new documentation!
- Child templates - easily customise specific parts of the default report template
- Plugin hooks - allow other tools to execute custom code during MultiQC execution
- New Preseq module
- New design for general statistics table (snazzy new background bars)
- Further development of toolbox
  - New button to clear all filters
  - Warnings when samples are hidden, plus empty plots and table cols are hidden
  - Active toolbar tab buttons are highlighted
- Lots of refactoring by @moonso to please the Pythonic gods
  - Switched to click instead of argparse to handle command line arguments
  - Code generally conforms to best practices better now.
- Now able to supply multiple directories to search for reports
- Logging output improved (now controlled by `-q` and `-v` for quiet and verbose)
- More HTML output dealt with by the base module, less left to the modules
  - Module introduction text
  - General statistics table now much easier to add to (new helper functions)
- Images, CSS and Javascript now included in HTML, meaning that there is a single
  report file to make sharing easier
- More accessible scrolling in the report - styled scrollbars and 'to top' button.
- Modules and templates now use setuptools entry points, facilitating plugins
  by other packages. Allows niche extensions whilst keeping the core codebase clean.
- The general stats table now has a sticky header row when scrolling, thanks to
  some new javascript wizardry...
- General stats columns can have a _shared key_ which allows common colour schemes
  and data ranges. For instance, all columns describing a read count will now share
  their scale across modules.
- General stats columns can be hidden and reordered with a new modal window.
- Plotting code refactored, reports with many samples (>50 by default) don't
  automatically render to avoid freezing the browser.
- Plots with highlighted and renamed samples now honour this when exporting to
  different file types.

---

## [MultiQC v0.2](https://github.com/MultiQC/MultiQC/releases/tag/v0.2) - 2015-09-18

- Code restructuring for nearly all modules. Common base module
  functions now handle many more functions (plots, config, file import)
  - See the [contributing notes](https://github.com/MultiQC/MultiQC/blob/master/CONTRIBUTING.md)
    for instructions on how to use these new helpers to make your own module
- New report toolbox - sample highlighting, renaming, hiding
  - Config is autosaved by default, can also export to a file for sharing
  - Interactive tour to help users find their way around
- New Tophat, Bowtie 2 and QualiMap modules
  - Thanks to @guillermo-carrasco for the QualiMap module
- Bowtie module now works
- New command line parameter `-d` prefixes sample names with the directory that
  they were found in. Allows duplicate filenames without being overwritten.
- Introduction walkthrough helps show what can be done in the report
- Now compatible with both Python 2 and Python 3
- Software version number now printed on command line properly, and in reports.
- Bugfix: FastQC doesn't break when only one report found
- Bugfix: FastQC seq content heatmap highlighting
- Many, many small bugfixes

---

## [MultiQC v0.1](https://github.com/MultiQC/MultiQC/releases/tag/v0.1) - 2015-09-01

- The first public release of MultiQC, after a month of development. Basic
  structure in place and modules for FastQC, FastQ Screen, Cutadapt, Bismark,
  STAR, Bowtie, Subread featureCounts and Picard MarkDuplicates. Approaching
  stability, though still under fairly heavy development.
