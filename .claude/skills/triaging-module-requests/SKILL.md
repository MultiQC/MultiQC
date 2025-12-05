# Triage MultiQC Module Requests

## Description

Analyze and prioritize MultiQC module requests through automated triage. Calculate priority scores, assign labels, update project boards, and provide actionable feedback to contributors.

**Use this skill when:**

- A new module request issue is opened (labeled `module: new`)
- User comments `@claude analyze-module` on a module request
- Weekly bulk triage is scheduled
- Manual analysis is requested via workflow dispatch

## Operation Modes

### analyze-single

Analyze one specific module request issue. Requires issue number.

**When to use:** New issues, on-demand analysis requests

### triage-all

Analyze all open module requests labeled `module: new`.

**When to use:** Weekly batch processing, cleanup operations

### dry-run

Perform analysis without making any GitHub changes (no labels, comments, or board updates).

**When to use:** Testing, validation, debugging

## Quick Start

1. **Determine the mode** from the workflow context or user request
2. **Fetch issue data** using GitHub CLI (`gh issue view` or `gh issue list`)
3. **Calculate priority score** using [scoring-criteria.md](scoring-criteria.md)
4. **Perform GitHub actions** following [github-actions.md](github-actions.md)
5. **Post analysis** using templates from [analysis-templates.md](analysis-templates.md)

## Priority Score Overview

Score is 0-100 based on five weighted categories:

- **Tool Popularity** (25 pts): GitHub metrics
- **Package Downloads** (15 pts): PyPI/Conda/Bioconda downloads
- **Community Engagement** (35 pts): Reactions, comments, duplicates
- **Request Quality** (20 pts): Completeness, example files
- **Technical Feasibility** (15 pts): Output parseability, metrics clarity

**Priority Bands:**

- 🔴 **High** (≥70): `priority: high` label
- 🟡 **Medium** (40-69): `priority: medium` label
- 🟢 **Low** (20-39): `priority: low` label
- ⚪ **Hold** (<20): `needs-triage` label only

See [scoring-criteria.md](scoring-criteria.md) for detailed rubric.

## Workflow Integration

This skill is invoked by `.github/workflows/module-requests.yml` which:

- Triggers on new issues with `module: new` label
- Responds to `@claude analyze-module` comments
- Runs weekly bulk triage (Mondays at 9 AM UTC)
- Supports manual workflow dispatch

## GitHub Operations

Key operations (see [github-actions.md](github-actions.md) for details):

- Fetch issue metadata and body content
- Extract tool information (name, URL, description)
- Add/update priority labels
- Post analysis comments
- Update project board positions (if configured)

## Tool Metrics Collection

Use `scripts/fetch-tool-metrics.js` for reliable API calls:

- GitHub stars, forks, last commit date
- PyPI download statistics
- Bioconda package data
- Repository activity metrics

## Analysis Output

Generate clear, actionable feedback using templates from [analysis-templates.md](analysis-templates.md):

- Current priority score with breakdown
- Specific improvement recommendations
- Comparison to similar requests
- Next steps for increasing priority

## Error Handling

- **Missing tool URL**: Assign low score, request homepage in feedback
- **Private/deleted repository**: Note in analysis, use partial scoring
- **API rate limits**: Implement exponential backoff, cache results
- **Invalid issue format**: Log warning, assign to "Needs Analysis" column

## Best Practices

1. **Be specific**: Point to exact fields that need improvement
2. **Be encouraging**: Frame feedback positively, emphasize what's good
3. **Be consistent**: Apply scoring rubric uniformly across all requests
4. **Be transparent**: Show score calculations in analysis comments
5. **Respect rate limits**: Cache API results, batch operations

## Files in This Skill

- `SKILL.md` (this file): High-level overview and workflow
- `scoring-criteria.md`: Detailed scoring rubric with examples
- `github-actions.md`: GitHub API operations and CLI commands
- `analysis-templates.md`: Comment templates and feedback patterns
- `scripts/fetch-tool-metrics.js`: Tool metrics collection script

## Related Documentation

- [Module Triage System Guide](../../docs/module-triage-system.md)
- [Project Board Setup](../../docs/module-triage-project-setup.md)
- [Workflow Configuration](../../../.github/workflows/module-requests.yml)
