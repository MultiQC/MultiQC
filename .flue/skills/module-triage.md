# module-triage

Triage a MultiQC module request. These are issues asking for support for a new
bioinformatics tool (typically labelled `module: new`).

## Inputs

The workflow passes a JSON payload containing:

- `number` — issue number
- `title` — issue title
- `body` — issue body
- `author` — submitter login
- `labels` — current label names

## Task

Read `../context.md` for the MultiQC architecture. Then analyze the request
against these dimensions:

1. **Tool identification** — what tool is being requested? Common name, primary
   purpose, output formats.
2. **Existing coverage** — is this tool already supported, or covered by a
   similar/sibling module?
3. **Adoption signal** — does the request indicate widespread use (e.g. cited
   pipeline, popular field, GitHub stars mentioned)?
4. **Implementation complexity** — based on the request, how hard would a
   module be? Consider: structured output (JSON/TSV) vs free-form, multiple
   sample handling, existing parsers in language ecosystems.
5. **Request completeness** — does the issue include example output files,
   links to docs, sample data? Or is it underspecified?

## Output

Suggest a priority label (`module: prio-low`, `module: prio-medium`,
`module: prio-high`) reflecting strength of adoption signal × implementation
feasibility. Suggest a concrete next step the maintainer or community can take.

## Schema

```ts
v.object({
  tool: v.object({
    name: v.string(),
    purpose: v.string(),
    output_formats: v.array(v.string()),
  }),
  existing_coverage: v.string(),
  adoption_signal: v.picklist(['weak','moderate','strong']),
  implementation_complexity: v.picklist(['low','medium','high']),
  request_completeness: v.picklist(['incomplete','partial','complete']),
  suggested_priority_label: v.picklist(['module: prio-low','module: prio-medium','module: prio-high']),
  summary: v.string(),
  next_step: v.string(),
})
```
