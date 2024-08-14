---
name: BBTools
urls: ["http://jgi.doe.gov/data-and-tools/bbtools/"]
summary: >
  Pre-processing, assembly, alignment, and statistics tools for DNA/RNA sequencing reads
---

The module produces summary statistics from the
[BBMap](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/) suite of tools.
The module can summarise data from the following BBMap output files
(descriptions from command line help output):

- `stats`
  - BBDuk filtering statistics.
- `covstats` _(not yet implemented)_
  - Per-scaffold coverage info.
- `rpkm` _(not yet implemented)_
  - Per-scaffold RPKM/FPKM counts.
- `covhist`
  - Histogram of # occurrences of each depth level.
- `basecov` _(not yet implemented)_
  - Coverage per base location.
- `bincov` _(not yet implemented)_
  - Print binned coverage per location (one line per X bases).
- `scafstats` _(not yet implemented)_
  - Statistics on how many reads mapped to which scaffold.
- `refstats`
  - Statistics on how many reads mapped to which reference file; only for BBSplit.
- `bhist`
  - Base composition histogram by position.
- `qhist`
  - Quality histogram by position.
- `qchist`
  - Count of bases with each quality value.
- `aqhist`
  - Histogram of average read quality.
- `bqhist`
  - Quality histogram designed for box plots.
- `lhist`
  - Read length histogram.
- `gchist`
  - Read GC content histogram.
- `indelhist`
  - Indel length histogram.
- `mhist`
  - Histogram of match, sub, del, and ins rates by read location.
- `statsfile` _(not yet implemented)_
  - Mapping statistics are printed here.

Additional information on the BBMap tools is available on
[SeqAnswers](http://seqanswers.com/forums/showthread.php?t=41057).

### File search patterns

```yaml
bbmap/aqhist:
  contents: "#Quality\tcount1\tfraction1\tcount2\tfraction2"
  num_lines: 1
bbmap/bhist:
  contents: "#Pos\tA\tC\tG\tT\tN"
  num_lines: 1
bbmap/bincov:
  contents: "#RefName\tCov\tPos\tRunningPos"
  num_lines: 3
bbmap/bqhist:
  contents: "#BaseNum\tcount_1\tmin_1\tmax_1\tmean_1\tQ1_1\tmed_1\tQ3_1\tLW_1\tRW_1\t\
    count_2\tmin_2\tmax_2\tmean_2\tQ1_2\tmed_2\tQ3_2\tLW_2\tRW_2"
  num_lines: 1
bbmap/covhist:
  contents: "#Coverage\tnumBases"
  num_lines: 1
bbmap/covstats:
  contents: "#ID\tAvg_fold"
  num_lines: 1
bbmap/ehist:
  contents: "#Errors\tCount"
  num_lines: 1
bbmap/gchist:
  contents: "#GC\tCount"
  num_lines: 5
bbmap/idhist:
  contents: "#Mean_reads"
  num_lines: 1
bbmap/ihist:
  contents: "#InsertSize\tCount"
  num_lines: 6
bbmap/indelhist:
  contents: "#Length\tDeletions\tInsertions"
  num_lines: 1
bbmap/lhist:
  contents: "#Length\tCount"
  num_lines: 1
bbmap/mhist:
  contents: "#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\tMatch2\tSub2\tDel2\t\
    Ins2\tN2\tOther2"
  num_lines: 1
bbmap/qahist:
  contents: "#Deviation"
  num_lines: 1
bbmap/qchist:
  contents_re: "#Quality\tcount1\tfraction1$"
  num_lines: 1
bbmap/qhist:
  contents: "#BaseNum\tRead1_linear\tRead1_log\tRead1_measured\tRead2_linear\tRead2_log\t\
    Read2_measured"
  num_lines: 1
bbmap/rpkm:
  contents: "#File\t"
  num_lines: 1
bbmap/stats:
  contents: "#Name\tReads\tReadsPct"
  num_lines: 4
bbmap/statsfile:
  contents: "Reads Used:"
  num_lines: 1
bbmap/statsfile_machine:
  contents: Reads Used=
  num_lines: 1
```
