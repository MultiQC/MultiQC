---
name: fgbio
urls: ["http://fulcrumgenomics.github.io/fgbio/"]
summary: >
  Processing and evaluating data containing UMIs
---

The module currently supports tool the following outputs:

- [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)
- [ErrorRateByReadPosition](http://fulcrumgenomics.github.io/fgbio/tools/latest/ErrorRateByReadPosition.html)

### File search patterns

```yaml
fgbio/errorratebyreadposition:
  contents: "read_number\tposition\tbases_total\terrors\terror_rate\ta_to_c_error_rate\t\
    a_to_g_error_rate\ta_to_t_error_rate\tc_to_a_error_rate\tc_to_g_error_rate\tc_to_t_error_rate"
  num_lines: 3
fgbio/groupreadsbyumi:
  contents: fraction_gt_or_eq_family_size
  num_lines: 3
```
