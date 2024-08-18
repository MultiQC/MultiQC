---
title: Contributing
description: Guides for how to contribute to the MultiQC code base
---

# Contributing

## Changelog

Almost all changes deserve an entry in the `CHANGELOG.md` file, so that people
know what updates are present between versions.

Whilst you can do this yourself by manually editing the file, we prefer to automate
the process by using our friendly MultiQC bot, just before merging.
By doing the changelog entry at the last minute we reduce the risk of having to
solve changelog merge conflicts.

The MultiQC changelog bot works by using the pull-request title.
**Your job is to ensure that your pull-request follows one of the following 3 conventions:**

- `New module: XYZ` - adding a new module named `XYZ`
- `XYZ: Change something in this existing module` - updating module `XYZ`
- `Some other change` - anything else, e.g. core MultiQC changes
- `Typo in docs [skip changelog]` - a change so minor that we don't want to log it at all

The MultiQC bot will automatically build a proper changelog entry based on this title
and (for new modules / module changes) the meta-information in the `MultiqcModule` class.

When a pull request is opened, a GitHub Action script is triggered, that inspects
the PR, updates the changelog and commits the update back to your PR. If your
pull request is not worth a change log entry (i.e. a minor documentation update),
you can append `[skip changelog]` to the PR title.

The action can also be triggered manually by adding the following comment on an open
pull request:

```md
@multiqc-bot changelog
```

It will replace the automatically added changelog entry if the pull request title
was updated after the initial commit. And finally, if you forgot to initially append
`[skip changelog]` and you did it after the initial commit, triggering the bot with a
comment will assure that the changelog line is removed.

## Admonitions

Admonitions, sometimes known as call-outs, can be used to highlight relevant information in the docs so that it stands out of the main flow of text.

The MultiQC website uses [remark-directive](https://github.com/remarkjs/remark-directive) to add a set of custom styled admonition directives to markdown syntax, shown below.

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

### Custom titles

```md
:::note{title="Don't Panic"}
Arthur Dent was grappling with his consciousness the way one grapples with a lost bar of soap in the bath. He lay, panting heavily in the wet air, and tried feeling bits of himself to see where he might be hurt. Wherever he touched himself, he encountered a pain. After a short while he worked out that this was because it was his hand that was hurting. Arthur nodded intelligently. Today was one of those bad days.
:::
```

:::note{title="Don't Panic"}
Arthur Dent was grappling with his consciousness the way one grapples with a lost bar of soap in the bath. He lay, panting heavily in the wet air, and tried feeling bits of himself to see where he might be hurt. Wherever he touched himself, he encountered a pain. After a short while he worked out that this was because it was his hand that was hurting. Arthur nodded intelligently. Today was one of those bad days.
:::
