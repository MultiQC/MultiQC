---
Name: BBMap
URL: http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/
Description: >
	BBMap is a suite of pre-processing, assembly, alignment, and statistics
	tools for DNA/RNA sequencing reads.
---

The BBMap module produces summary statistics from the
[BBMap](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/) suite of tools.
The module can summarise data from the following BBMap output files
(descriptions from `bbmap.sh` help output):

* `covstats` _(not yet implemented)_
    * Per-scaffold coverage info.
* `rpkm` _(not yet implemented)_
    * Per-scaffold RPKM/FPKM counts.
* `covhist`
    * Histogram of # occurrences of each depth level.
* `basecov` _(not yet implemented)_
    * Coverage per base location.
* `bincov` _(not yet implemented)_
    * Print binned coverage per location (one line per X bases).
* `scafstats`
    * Statistics on how many reads mapped to which scaffold.
* `refstats`
    * Statistics on how many reads mapped to which reference file; only for BBSplit.
* `bhist`
    * Base composition histogram by position.
* `qhist`
    * Quality histogram by position.
* `qchist`
    * Count of bases with each quality value.
* `aqhist`
    * Histogram of average read quality.
* `bqhist`
    * Quality histogram designed for box plots.
* `lhist`
    * Read length histogram.
* `gchist`
    * Read GC content histogram.
* `indelhist`
    * Indel length histogram.
* `mhist`
    * Histogram of match, sub, del, and ins rates by read location.
* `statsfile` _(not yet implemented)_
    * Mapping statistics are printed here.

Additional information on the BBMap tools is available on
[SeqAnswers](http://seqanswers.com/forums/showthread.php?t=41057).
