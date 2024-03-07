---
name: Nonpareil
url: https://github.com/lmrodriguezr/nonpareil
description: Estimate metagenomic coverage and sequence diversity.
---

Nonpareil uses the redundancy of the reads in a metagenomic dataset to estimate
the average coverage and predict the ammount of sequences that will be required
to achieve "nearly complete coverage", defined as ≥95% or ≥99% average coverage.

Since Nonpareil main output has no model information, it is necessary to run its
auxiliary `R` plot functions and save the `curves` object as a `JSON` file. Briefly,
call function `export_curve()` on object `curves` (for an example, see [snakemake wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/nonpareil/plot.html#code)):

```r
base::library("jsonlite")
base::message("Exporting model as JSON")

export_curve <- function(object){
  # Extract variables
  n <- names(attributes(object))[c(1:12,21:29)]
  x <- sapply(n, function(v) attr(object,v))
  names(x) <- n
  # Extract vectors
  n <- names(attributes(object))[13:20]
  y <- lapply(n, function(v) attr(object,v))
  names(y) <- n
  curve_json <- c(x, y)

  # Add model
  if (object$has.model) {
    # https://github.com/lmrodriguezr/nonpareil/blob/162f1697ab1a21128e1857dd87fa93011e30c1ba/utils/Nonpareil/R/Nonpareil.R#L330-L332
    x_min <- 1e3
    x_max <- signif(tail(attr(object,"x.adj"), n=1)*10, 1)
    x.model <- exp(seq(log(x_min), log(x_max), length.out=1e3))
    y.model <- predict(object, lr=x.model)
    curve_json <- append(curve_json, list(x.model=x.model))
    curve_json <- append(curve_json, list(y.model=y.model))
  }

  base::print(curve_json)
  curve_json
}

export_set <- function(object){
  y <- lapply(object$np.curves, "export_curve")
  names(y) <- sapply(object$np.curves, function(n) n$label)
  jsonlite::prettify(toJSON(y, auto_unbox=TRUE))
}

y <- export_set(curves)
write(y, "output.json")
```

### Module config options

The module plots a line graph for each sample, with a tab panel to switch between only observed data, only models,
or both combined (model with a dashed line). It will use the colors specified in the JSON file by `nonpareil` and,
if some is missing use one from a MultiQC colour scheme (default: Paired) that can be defined with:

```yaml
nonpareil:
  plot_colours: Paired
```
