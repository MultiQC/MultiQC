---
title: Contributing
description: Guides for how to contribute to the MultiQC code base
---

# Contributing

## Changelog

`CHANGELOG.md` file is populated semi-automatically. To generate an initial state, we use a script:

```sh
python scripts/print_changelog.py
```

It automatically generates the changelog from the merged pull-requests assigned to a milestone (e.g. `v1.26`),
using titles as changelog entries, and tags (labels) to categorize entries in sections (e.g. `New modules`, `Module fixes`, `Infrastructure`, etc.). We run that script before creating a release.

For that reason, **your job is to ensure that your pull-request has a clean thoughtful title**. The title must
summarize the changes, and be written in a good English without typos and errors, start with a capital letter,
and have single spaces between words. Examples of good titles:

- "Prepend plot IDs with `self.anchor` to assure custom anchor is applied"
- "New module: Percolator (semi-supervised learning framework for peptide identification)"
- "GATK BQSR: support Sentieon QualCal output"

Add labels to the PR to automatically categorize them in the changelog:

- `module: new` - "New modules"
- `module: enhancement`, `module: change` - "Module updates"
- `bug: module` - "Module fixes"
- `bug: core` - "Fixes"
- `core: back end`, `core: front end` - "Feature updates and improvements"
- `core: infrastructure` - "Infrastructure and packaging"
- `core: refactoring` - "Refactoring and typing"
- `documentation` - "Chores"

## Docs - Admonitions

Admonitions, sometimes known as call-outs, can be used to highlight relevant information in the docs so that it stands out of the main flow of text.

### Notes

```md
:::note
He had half a mind just to keep on `falling`.
:::
```

:::note
He had half a mind just to keep on `falling`.
:::

### Info

```md
:::info
His face froze for a second or two and then began to do that terribly slow crashing `trick` that Arctic ice floes do so spectacularly in the spring.
:::
```

:::info
His face froze for a second or two and then began to do that terribly slow crashing `trick` that Arctic ice floes do so spectacularly in the spring.
:::

### Tip

```md
:::tip
Her remark would have commanded greater attention had it been generally realized that human beings were only the third most intelligent life form present on the planet Earth.
:::
```

:::tip
Her remark would have commanded greater attention had it been generally realized that human beings were only the third most intelligent life form present on the planet Earth.
:::

### Success

```md
:::success
“I’m afraid you cannot leave,' said Zarniwoop, 'you are entwined in the Improbability Field. You cannot escape.'
:::
```

:::success
“I’m afraid you cannot leave,' said Zarniwoop, 'you are entwined in the Improbability Field. You cannot escape.'
:::

### Warnings

```md
:::warning
He smiled the smile that Zaphod had wanted to hit and this time `Zaphod` hit it.
:::
```

:::warning
He smiled the smile that Zaphod had wanted to hit and this time `Zaphod` hit it.
:::

### Danger

```md
:::danger
One of the troublesome circumstances was the Plural nature of this Galactic Sector, where the possible `continually` interfered with the probable.
:::
```

:::danger
One of the troublesome circumstances was the Plural nature of this Galactic Sector, where the possible `continually` interfered with the probable.
:::
