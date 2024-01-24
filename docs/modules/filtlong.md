---
name: Filtlong
url: https://github.com/rrwick/Filtlong
description: Filtlong is a tool for filtering long reads by quality.
---

[Filtlong](https://github.com/rrwick/Filtlong) is a tool for filtering long reads by quality.
It can take a set of long reads and produce a smaller, better subset. It uses both read length (longer is better) and read identity (higher is better) when choosing which reads pass the filter.

The module takes summary statistics of number of long reads filtered and displays them in the General Stats table.

### Bases Kept

Sometimes, the Filtlong log message contains this:

```
Filtering long reads
  target: 123456789 bp
  reads already fall below target after filtering

Outputting passed long reads
```

In these cases we cannot say for sure how many bases were kept. As such, this field is left blank.
If you have a better solution, please suggest in an issue or pull request.
