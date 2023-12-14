---
name: Nonpareil
url: https://github.com/lmrodriguezr/nonpareil
description: >
  Estimate metagenomic coverage and sequence diversity.
---

Nonpareil uses the redundancy of the reads in a metagenomic dataset to estimate
the average coverage and predict the ammount of sequences that will be required
to achieve "nearly complete coverage", defined as ≥95% or ≥99% average coverage.

Since Nonpareil main output has no model information, it is necessary to run its
auxiliary `R` plot functions and save the `curves` object as a `JSON` file. For
more info on how to generate this file, check `nonpareil` [snakemake plot wrapper](https://snakemake-wrappers.readthedocs.io/en/latest/wrappers/nonpareil/plot.html#code). Briefly, call function `export_curve()` on object `curves`:

```r
y <- export_set(curves)
write(y, "output.json")
```

### Module config options

The module builds a line graph for each model with a tab panel to switch between models.
The plot will contain a line for the observed data and a line of the inferred model data.
You can hide one of those lines with the config options below:

```yaml
nonpareil:
  plot_observed: true # default: true
  plot_model: false # default: true
```

You can also add extra dashed lines to represent the dispersion of the model:

```yaml
nonpareil:
  plot_dispersion: sd # default: none
```

The available values for `plot_dispersion`:

- `sd`: one standard deviation around the mean,
- `ci95`: confidence interval at 95%,
- `ci90`: confidence interval at 90%,
- `ci50`: confidence interval at 50%,
- `iq`: inter-quartile range.
