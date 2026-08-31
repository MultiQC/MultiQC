# MultiQC ECharts Template

A child template of `default` that renders report plots with
[Apache ECharts](https://echarts.apache.org/) instead of Plotly. Both plotting engines are
fully supported and selectable; the `default` template and Plotly remain the default for
everyone who doesn't ask for `echarts`.

## Usage

```bash
multiqc . --template echarts
```

Or in your MultiQC configuration file:

```yaml
template: echarts
```

All 7 plot types (bar, line, scatter, box, violin, heatmap and the general stats
table/violin) render with ECharts under this template, including the interactive toolbox
(highlight, rename, hide), dataset/percentage/log10 switches, in-report PNG/SVG export, and
flat/exported image plots.

## How selection works

Selecting the template sets an internal config value, `config.plotting_engine`, to
`"echarts"` (the default is `"plotly"`). This happens automatically when you pass
`--template echarts`; there is no separate `--plotting-engine` flag. A related setting,
`config.echarts_canvas_threshold` (default `3000`), controls how many marks (points, bars,
heatmap cells) a single plot needs before ECharts switches from its SVG renderer to its
canvas renderer for performance. You can override either in a MultiQC config file:

```yaml
plotting_engine: "echarts"
echarts_canvas_threshold: 5000
```

## Static (flat) plot export

Server-side rendering of ECharts plots (used for `--flat`, `--export-plots`, and any plot
that falls back to a static image) needs two extra Python packages: `mini-racer` (an
embedded V8 engine that runs the ECharts JS bundle server-side) and `resvg-py` (rasterises
the resulting SVG to PNG). Install them with:

```bash
pip install 'multiqc[echarts]'
```

If these packages are missing and a flat/exported ECharts plot is requested, MultiQC raises
a loud error pointing at this install command rather than silently falling back to
something else.

## Fonts

Interactive (in-browser) ECharts plots use the CSS font stack from the page, same as any
other web content. Static export is different: `resvg` (the SVG-to-PNG rasteriser) needs a
single, exact font family name, not a comma-separated CSS stack, because ECharts writes
`font-family` straight into the exported SVG's `style` attribute without escaping inner
quotes. A stack like `Inter, "Helvetica Neue", Arial` produces malformed SVG that `resvg`
refuses to parse. For this reason:

- Static/flat exports always use a single bundled font (DejaVu Sans), regardless of any
  `plot_font_family` setting, so exported images render identically on every machine.
- If you set `config.plot_font_family` to customise interactive plots under this template,
  use a single family name (for example `"Courier New"`), not a stack. A multi-family CSS
  stack is not supported by the echarts template and will not behave as expected.

## Known gaps compared to the Plotly (default) template

The ECharts backend covers the vast majority of Plotly behaviour, but a handful of
differences are accepted rather than fixed, either because closing them would need
disproportionate effort for a minor cosmetic difference, or because they touch things
outside the template (module code, the notebook API). These are documented here rather than
silently glossed over:

- **PDF export is not supported.** `--export-plots` with `pdf` in
  `config.export_plot_formats` logs a warning ("PDF export is not supported by the ECharts
  engine, skipping") and skips that format; PNG and SVG export both work. `--pdf` forces the
  `simple` template (which always uses Plotly), so it is unaffected either way.
- **The notebook API (`Plot.show()` / `Plot.save()`) is always Plotly.** These methods build
  and return a Plotly `go.Figure` object regardless of `config.plotting_engine`, so scripts
  or notebooks that call them will always see a Plotly figure, even when
  `--template echarts` is used elsewhere in the run.
- **Toolbox highlight does not recolor sample axis tick labels.** Plotly colors the sample
  names on the axis itself when a highlight pattern matches; ECharts ignores the
  `tickvals`/`ticktext` styling Plotly uses for this, so axis labels stay their normal color.
  Dimming the non-highlighted bars/points/lines remains the primary highlight signal and
  works the same as under the default template.
- **Toolbox rename does not relabel line plot series names**, matching a pre-existing
  quirk of the default template's line plot code (not something introduced by this
  template).
- **Heatmap colour scales are an approximation of custom colour stops.** ECharts'
  `visualMap` interpolates colours evenly between its stops; MultiQC's default 11-stop
  colour scale is even, so it renders identically, but a module-supplied `colstops` list
  with uneven stop positions loses those exact positions (colours are still correct, spacing
  between them is not pixel-identical to the Plotly version).
- **Scatter plots don't show a "N samples hidden" banner** when the toolbox hides samples
  (the bar plot does show this banner). Hidden points are still dropped correctly; only the
  banner is missing.
- **Per-dataset overrides of `x_bands`/`y_bands`/`x_lines`/`y_lines`** (set via
  `pconfig["data_labels"][i]["x_bands"]`, etc.) are not read; only the top-level `pconfig`
  values are used. This is the same class of gap as the general "multi-dataset override"
  pattern and affects a small number of modules.
- **`tickformat` on axes is not converted** (the related `hoverformat`, `ticksuffix` and
  numeric tooltip precision are handled). No shipped module currently relies on
  `tickformat`, so this is a low-risk gap.
- **Heatmaps have no click+drag box-zoom.** Bar, box, line, and scatter plots all support
  Plotly-style rubber-band zoom (drag to zoom, double-click to reset; see
  `converter.zoom_option`). Heatmaps don't: ECharts' toolbox `dataZoomSelect` brush never
  reliably engages for a `heatmap` series in this build, confirmed with agent-browser
  across canvas/SVG renderer and square/non-square layouts, so rather than ship a control
  that silently does nothing, heatmap keeps its pre-existing no-zoom state. Also disabled
  outright, matching Plotly (its axes aren't real data coordinates): violin plots.

None of the above are silent: PDF export logs a warning, and the notebook limitation is
documented here since there is nowhere sensible to log it at the point a user hits it.

## Rebuilding the JavaScript bundles

This template ships **two** separate build outputs, and both need to be rebuilt (and
committed) after changing their sources:

1. **The vendored ECharts library itself** (a tree-shaken, minified build of ECharts 6.1.0
   containing only the chart types and components MultiQC needs): built from
   `echarts-entry.js` by `build-echarts.mjs` (an esbuild script) into
   `assets/js/packages/echarts-6.1.0.custom.min.js`. Rebuild it with:

   ```bash
   cd multiqc/templates/echarts
   npm install
   npm run build:echarts
   ```

   Only rerun this if you change which ECharts chart types or components are imported in
   `echarts-entry.js`, or bump the pinned `echarts` version in `package.json`.

2. **The template's own application JavaScript and CSS** (the `Plot` classes, toolbox
   export, table forks, and Sass), bundled with Vite into `compiled/`: rebuild it with:

   ```bash
   cd multiqc/templates/echarts
   npm install
   npm run build
   ```

   Run this after any change under `src/js/` or a change to the parent `default` template's
   Sass that this template imports.

Both `compiled/` and the vendored ECharts bundle are committed to the repository, the same
way the `default` and `disco` templates commit their build output, so that installing
MultiQC from a wheel does not require a Node.js toolchain.

For live development, use `multiqc ... --template echarts --development` together with
`npm run watch` in this directory: this links the source JS/CSS files instead of inlining
them, so a browser refresh (no rebuild) picks up your changes.
