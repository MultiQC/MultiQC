---
name: Kaiju
urls: ["http://kaiju.binf.ku.dk/"]
summary: >
  Taxonomic classification for metagenomics
---

The module parses output generated by kaiju2table, e.g:

```bash
kaiju -i R1.fq.gz -j R2.fq.gz -o output_kaiju.txt
kaiju2table -t nodes.dmp -n names.dmp -r species -o kaiju2table_species.txt output_kaiju.txt
kaiju2table -t nodes.dmp -n names.dmp -r phylum -o kaiju2table_phylum.txt output_kaiju.txt
```

### File search patterns

```yaml
kaiju:
  contents_re: file\tpercent\treads\ttaxon_id\ttaxon_name
  num_lines: 1
```