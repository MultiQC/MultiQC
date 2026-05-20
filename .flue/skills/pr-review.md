# pr-review

Review a MultiQC pull request and produce a structured review.

## Inputs

The workflow passes a JSON payload containing:

- `number` — PR number
- `title` — PR title
- `body` — PR description
- `author` — author login
- `diff` — unified diff of the PR (truncated if very large)
- `changed_files` — list of changed file paths

## Task

Read `../context.md` for repo conventions. Then review the diff and produce:

- A one-paragraph **summary** of what the PR changes.
- A list of **findings**, each tagged by category (`bug`, `style`, `perf`,
  `security`, `test-coverage`, `docs`) with severity (`low`, `medium`, `high`)
  and a short explanation tied to a file/line where possible.
- **Suggestions** — concrete, actionable next steps.
- An overall **verdict** — `approve`, `comment`, or `request-changes`.

Be constructive. Reference MultiQC conventions from `context.md` when
applicable (e.g. `ModuleNoSamplesFound`, `add_software_version`,
`write_data_file` ordering).

## Schema

```ts
v.object({
  summary: v.string(),
  findings: v.array(v.object({
    category: v.picklist(['bug','style','perf','security','test-coverage','docs']),
    severity: v.picklist(['low','medium','high']),
    file: v.optional(v.string()),
    note: v.string(),
  })),
  suggestions: v.array(v.string()),
  verdict: v.picklist(['approve','comment','request-changes']),
})
```
