---
name: triaging-module-requests
description: |
  Triage MultiQC `module: new` GitHub issues: calculate 0-100 priority scores, apply priority labels, post analysis comments with score breakdowns, and give contributors actionable feedback to improve their request. Use when a new `module: new` issue is opened, when a user comments `@claude analyze-module` on a request, during weekly bulk triage, or when manually re-evaluating a module request.
---

# Triage MultiQC Module Requests

Invoked by `.github/workflows/module-requests.yml` on new `module: new`
issues, on `@claude analyze-module` comments, and on a weekly schedule
(Mondays 9 AM UTC). Also runs on manual workflow dispatch.

## Workflow

1. **Pick a mode** from context:
   - `analyze-single` — one issue (new issue or on-demand request)
   - `triage-all` — every open `module: new` issue (weekly batch)
   - `dry-run` — calculate and print, but make no GitHub changes
2. **Fetch issue data** with `gh issue view` / `gh issue list`. See
   [github-actions.md](github-actions.md) for the exact commands.
3. **Fetch tool metrics** via `scripts/fetch-tool-metrics.js` (GitHub stars,
   PyPI/Conda downloads, last commit date).
4. **Calculate the priority score** using the rubric in
   [scoring-criteria.md](scoring-criteria.md).
5. **Apply the priority label** (see Priority bands below). Remove any
   existing priority labels first; see [github-actions.md](github-actions.md).
6. **Post the analysis comment** using
   [analysis-templates.md](analysis-templates.md). Always show the score
   breakdown, what's good, and concrete improvements (with `+X points`
   tags where they apply).

## Priority bands

| Score | Band      | Label                 |
| ----- | --------- | --------------------- |
| ≥70   | 🔴 High   | `module: prio-high`   |
| 40–69 | 🟡 Medium | `module: prio-medium` |
| 20–39 | 🟢 Low    | `module: prio-low`    |
| <20   | ⚪ Hold   | `module: prio-hold`   |

The five score categories (full rubric in [scoring-criteria.md](scoring-criteria.md)):

- Tool Popularity (25) — GitHub stars + maintenance bonus
- Package Downloads (15) — PyPI / Conda / Bioconda monthly
- Community Engagement (35) — reactions, comments, duplicates
- Request Quality (20) — completed fields + example files
- Technical Feasibility (15) — output format, metric clarity, parsing

## Feedback principles

Be specific (point to exact fields), be encouraging (acknowledge strengths
first), be consistent (apply the rubric uniformly), be transparent (show
the calculation). Cache API results when batch-processing to stay under
rate limits.
