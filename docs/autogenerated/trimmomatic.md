---
name: Trimmomatic
urls: ["http://www.usadellab.org/cms/?page=trimmomatic"]
summary: >
  Read trimming tool for Illumina NGS data
---

The module parses the stderr output, that can be captured by directing it to a file e.g.:

```sh
trimmomatic command 2> trim_out.log
```

By default, the module generates the sample names based on the input FastQ file names in
the command line used by Trimmomatic. If you prefer, you can tell the module to use
the filenames as sample names instead. To do so, use the following config option:

```yaml
trimmomatic:
  s_name_filenames: true
```

### File search patterns

```yaml
trimmomatic:
  contents: Trimmomatic
```
