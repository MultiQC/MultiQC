---
name: Snippy
urls: ["https://github.com/tseemann/snippy"]
summary: >
  Rapid haploid variant calling and core genome alignment
---

The following commands are implemented:

- `snippy`
  - Variant type descriptive statistics.
  - Parses summary `prefix.txt` file that is generated.
- `snippy-core`
  - Core genome alignment descriptive statistics.
  - Parses summary `prefix.txt` file that is generated.

### File search patterns

```yaml
snippy/snippy:
  contents: snippy
  num_lines: 20
snippy/snippy-core:
  contents_re: ID\tLENGTH\tALIGNED\tUNALIGNED\tVARIANT\tHET\tMASKED\tLOWCOV
  num_lines: 1
```
