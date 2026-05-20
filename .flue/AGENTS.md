# MultiQC Agents

Background agents that run on GitHub events. Inspired by the
[software-factory.dev](https://www.software-factory.dev/) model: humans steer
intent through specs and labels, agents execute on events without `@mention`
invocation.

## Design rules

- **One trigger, one agent.** Each agent listens to one well-defined GitHub
  event class (issue opened, PR opened, etc.).
- **Structured output, no freeform.** Every agent declares a `valibot` schema
  so downstream workflow steps can act on fields, not parse prose.
- **Specs live in markdown.** Prompts and repo context belong in this
  directory's `.md` files. TypeScript files are thin glue: load context, pass
  payload to the model, return typed data.
- **Workflows handle GitHub I/O.** Agents do not call the GitHub API directly.
  The calling workflow fetches inputs (issue body, PR diff) and posts outputs
  (comments, labels) using `gh`. This keeps agents portable across hosts.
- **Pinned model + pinned deps.** No floating versions in `package.json` or
  `init({ model })`. CI runs are reproducible.

## Agents

| Agent           | Trigger                                                        | Output                                                    |
| --------------- | -------------------------------------------------------------- | --------------------------------------------------------- |
| `triage`        | `issues: opened`                                               | category, suggested labels, priority, next steps          |
| `pr-review`     | `pull_request: opened`, `/review` comment                      | summary, findings by category, suggestions                |
| `module-triage` | `issues: opened`/`labeled` with `module: new`, manual dispatch | adoption signal, complexity, completeness, priority label |

## Agent chains

Agents trigger each other via GitHub labels, with no extra dispatch surface;
every hop is visible in the GitHub UI.

```
issues:opened ─► triage ─[category=module-request]─► label "module: new"
                                                            │
                                                            ▼
                                  issues:labeled ─► module-triage ─► label "module: prio-{low,medium,high}"
```

Rules for adding a chain:

- Use **labels** as the transport, not `repository_dispatch` or
  `workflow_dispatch`. Labels are debuggable and resumable.
- The upstream agent must apply only labels from a **curated allow-list**, not
  raw model output, to avoid creating arbitrary labels.
- The downstream agent gates on the label via `if:` so other label changes
  don't fire it accidentally.
- Avoid loops: downstream agents must not apply labels that the upstream agent
  listens for.

## Adding a new agent

1. Write a spec in `skills/<name>.md` describing what the agent does and the
   exact output schema.
2. Create `agents/<name>.ts` (import shared context, call `init`/`session.prompt`
   with the schema).
3. Create `.github/workflows/flue-<name>.yml` that fetches inputs via `gh`,
   invokes `flue run`, and posts outputs.
4. Update the table above.

## Shared files

- `context.md` — MultiQC repo context loaded by every agent
- `skills/<name>.md` — per-agent prompt and schema spec
