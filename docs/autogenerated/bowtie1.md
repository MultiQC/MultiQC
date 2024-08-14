---
name: Bowtie 1
urls: ['http://bowtie-bio.sourceforge.net/']
summary: >
  Ultrafast, memory-efficient short read aligner
---

### File search patterns

```yaml
bowtie1:
  contents: '# reads processed:'
  exclude_fn:
  - bowtie.left_kept_reads.log
  - bowtie.left_kept_reads.m2g_um.log
  - bowtie.left_kept_reads.m2g_um_seg1.log
  - bowtie.left_kept_reads.m2g_um_seg2.log
  - bowtie.right_kept_reads.log
  - bowtie.right_kept_reads.m2g_um.log
  - bowtie.right_kept_reads.m2g_um_seg1.log
  - bowtie.right_kept_reads.m2g_um_seg2.log
  shared: true
```