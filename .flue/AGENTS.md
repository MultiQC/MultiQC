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

| Agent | Trigger | Output |
|---|---|---|
| `triage` | `issues: opened` | category, suggested labels, priority, next steps |
| `pr-review` | `pull_request: opened`, `/review` comment | summary, findings by category, suggestions |

## Adding a new agent

1. Write a spec in `skills/<name>.md` describing what the agent does and the
   exact output schema.
2. Create `agents/<name>.ts` — import shared context, call `init`/`session.prompt`
   with the schema.
3. Create `.github/workflows/flue-<name>.yml` that fetches inputs via `gh`,
   invokes `flue run`, and posts outputs.
4. Update the table above.

## Shared files

- `context.md` — MultiQC repo context loaded by every agent
- `skills/<name>.md` — per-agent prompt and schema spec
