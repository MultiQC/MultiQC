---
name: DeDup
urls: ["http://www.github.com/apeltzer/DeDup"]
summary: >
  Improved Duplicate Removal for merged/collapsed reads in ancient DNA analysis
---

By default, tables show read counts in thousands.
To customise this, you can set the following MultiQC config variables:

```yaml
ancient_read_count_prefix: "K"
ancient_read_count_desc: "thousands"
ancient_read_count_multiplier: 0.001
```

### File search patterns

```yaml
dedup:
  fn: "*dedup.json"
```
