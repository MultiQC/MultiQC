---
name: GLIMPSE
urls: ["https://odelaneau.github.io/GLIMPSE/"]
summary: >
  Low-coverage whole genome sequencing imputation
---

The program `GLIMPSE2` is based on the GLIMPSE model and designed for reference panels containing
hundreds of thousands of reference samples, with a special focus on rare variants.
The concordance rates values are displayed in a scatter plot, with the option to switch between
the different concordance metrics.

The supported files are generated from the `GLIMPSE2_concordance` command. The following files are supported:

- `*.error.spl.txt`
- `*.error.grp.txt`

### File search patterns

```yaml
glimpse/err_grp:
  fn: "*.error.grp.txt.gz"
  num_lines: 1
glimpse/err_spl:
  fn: "*.error.spl.txt.gz"
  num_lines: 1
```