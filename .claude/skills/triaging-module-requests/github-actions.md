# GitHub Actions Guide for Module Triage

This document describes GitHub CLI operations and API interactions for the module triage system.

## Prerequisites

Ensure GitHub CLI (`gh`) is authenticated with appropriate permissions:

- `issues: write` - Add labels, post comments
- `projects: write` - Update project board
- `contents: read` - Access repository data

## Fetching Issue Data

### Get Single Issue

```bash
# Fetch full issue metadata
gh issue view ISSUE_NUMBER --json title,body,labels,reactions,comments,author,createdAt,updatedAt

# Example output:
{
  "title": "New module for FastQC",
  "body": "...",
  "labels": [{"name": "module: new"}],
  "reactions": {"+1": 5, "heart": 2},
  "comments": [...],
  "author": {"login": "username"},
  "createdAt": "2025-10-20T10:00:00Z",
  "updatedAt": "2025-10-22T15:30:00Z"
}
```

### List All Module Requests

```bash
# Get all open issues with "module: new" label
gh issue list \
  --label "module: new" \
  --state open \
  --json number,title,labels,reactions,comments,createdAt \
  --limit 100

# Filter by additional criteria
gh issue list \
  --label "module: new" \
  --state open \
  --search "NOT label:\"priority: high\""
```

### Search for Related Issues

```bash
# Find duplicate or related requests
gh issue list \
  --search "TOOL_NAME in:title label:\"module: new\"" \
  --state all \
  --json number,title,state
```

## Parsing Issue Body

Issue bodies usually follow the template in `.github/ISSUE_TEMPLATE/module-request.yml`. Extract fields:

```bash
# Full body text
BODY=$(gh issue view ISSUE_NUMBER --json body --jq '.body')

# Extract specific fields using pattern matching
TOOL_NAME=$(echo "$BODY" | sed -n '/### Name of the tool/,/###/p' | sed '1d;$d' | xargs)
TOOL_URL=$(echo "$BODY" | sed -n '/### Tool homepage or repository/,/###/p' | sed '1d;$d' | xargs)
DESCRIPTION=$(echo "$BODY" | sed -n '/### Tool description/,/###/p' | sed '1d;$d' | xargs)
```

**Better approach:** Use Claude's natural language understanding to extract structured data from the issue body directly.

## Checking for Example Files

Example files can be provided in multiple ways:

1. Uploaded directly to the issue (attachments)
2. Linked to or submitted in the MultiQC/test-data repository

```bash
# Check for uploaded files (attachments)
if echo "$BODY" | grep -q "!\[.*\](.*)" || echo "$BODY" | grep -q "\[.*\..*\](.*)" ; then
  echo "Has uploaded files"
fi

# Check for references to test-data repository
if echo "$BODY" | grep -q "github.com/MultiQC/test-data"; then
  echo "References test-data repository"
fi
```

**Scoring guideline:** Files uploaded to the issue or properly linked to the test-data repository should receive full points.

## Managing Labels

### Add Priority Labels

**IMPORTANT: Always remove existing priority labels before adding new ones to avoid label conflicts.**

```bash
# Step 1: Remove ALL existing priority labels first
gh issue edit ISSUE_NUMBER \
  --remove-label "waiting: example data" \
  --remove-label "module: prio-hold" \
  --remove-label "module: prio-low" \
  --remove-label "module: prio-medium" \
  --remove-label "module: prio-high" 2>/dev/null || true

# Step 2: Add the appropriate new priority label
gh issue edit ISSUE_NUMBER --add-label "module: prio-high"
```

**Note:** The `2>/dev/null || true` suppresses errors if a label doesn't exist on the issue.

### Label Reference

**Priority labels:**

- `module: prio-high` - Score ≥ 70
- `module: prio-medium` - Score 40-69
- `module: prio-low` - Score 20-39
- `module: prio-hold` - Score <= 20

**Status labels:**

- `waiting: example data` - Missing example files

### Add Status Labels

```bash
# Mark as needing examples
gh issue edit ISSUE_NUMBER --add-label "waiting: example data"
# Add priority label
gh issue edit ISSUE_NUMBER --add-label "module: prio-medium"
```

## Posting Comments

### Analysis Comment

```bash
# Post analysis results
gh issue comment ISSUE_NUMBER --body "$(cat <<'EOF'
## 📊 Module Request Analysis

[See analysis-templates.md for full template]
EOF
)"
```

### Using Heredoc for Multi-line Comments

```bash
ANALYSIS=$(cat <<'EOF'
Thanks for requesting a new MultiQC module! This is an automated triage review to help prioritise development work.

| Item               | Details                                      |
| ------------------ | -------------------------------------------- |
| **Tool**           | FastQC                                       |
| **Repository**     | https://github.com/s-andrews/FastQC (⭐ 450) |
| **Priority Score** | 75/100 🔴 **High Priority**                  |

<details>

### Score Breakdown

| Category                 | Score | Notes                             |
| ------------------------ | ----- | --------------------------------- |
| 🌟 Tool Popularity       | 20/25 | 450 stars, actively maintained    |
| 📦 Package Downloads     | 12/15 | 50K downloads/month on Bioconda   |
| 💬 Community Engagement  | 18/35 | 8 👍 reactions, 3 comments         |
| ✅ Request Quality       | 20/20 | Complete info, example files      |
| ⚙️ Technical Feasibility | 5/15  | Complex output format             |

**Total: 75/100**

### 🔴 High Priority: Ready for Development

This is an excellent module request. The tool is popular, actively maintained, and the request is complete with example files.

**What's great:**

- ✨ Comprehensive example files provided
- ✨ Clear description of expected metrics
- ✨ Widely used tool in the community

### Next Steps

No action needed. Follow this issue for progress updates.

---

<sup>

This analysis was performed automatically. Comment `@claude analyze-module` for re-evaluation.
Learn more: [Module Triage System Guide](https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md)

</sup>

</details>
EOF
)

gh issue comment ISSUE_NUMBER --body "$ANALYSIS"
```

## Project Board Operations

### Get Project Information

```bash
# List organization projects
gh project list --owner MultiQC

# Get project fields and options
gh project field-list PROJECT_NUMBER --owner MultiQC
```

### Add Issue to Project

```bash
# Add issue to project (returns item ID)
ITEM_ID=$(gh project item-add PROJECT_NUMBER \
  --owner MultiQC \
  --url "https://github.com/MultiQC/MultiQC/issues/ISSUE_NUMBER")
```

### Update Project Fields

```bash
# Update custom fields (e.g., priority score)
gh project item-edit \
  --project-id PROJECT_ID \
  --id $ITEM_ID \
  --field-id FIELD_ID \
  --text "75"

# Move to column (using status field)
gh project item-edit \
  --project-id PROJECT_ID \
  --id $ITEM_ID \
  --field-id STATUS_FIELD_ID \
  --option-id "High Priority"
```

**Note:** Project board updates are optional. Labels alone provide sufficient tracking.

## Batch Operations

### Triage All Open Requests

```bash
# Get all module requests
ISSUES=$(gh issue list \
  --label "module: new" \
  --state open \
  --json number \
  --jq '.[].number')

# Process each
for ISSUE_NUMBER in $ISSUES; do
  echo "Analyzing #$ISSUE_NUMBER..."
  # [Run analysis logic]
done
```

### Update Stale Requests

```bash
# Find requests with no activity in 180 days
gh issue list \
  --label "module: new" \
  --search "updated:<$(date -v-180d +%Y-%m-%d)" \
  --json number,title
```

## Rate Limiting

GitHub API has rate limits. Check remaining calls:

```bash
# Check rate limit status
gh api rate_limit

# Response shows remaining calls per endpoint
{
  "resources": {
    "core": {
      "limit": 5000,
      "remaining": 4823,
      "reset": 1640000000
    }
  }
}
```

**Best practices:**

- Cache results when processing multiple issues
- Use `--json` queries to minimize API calls
- Implement exponential backoff on errors
- Batch operations when possible

## Error Handling

### Issue Not Found

```bash
if ! gh issue view ISSUE_NUMBER &>/dev/null; then
  echo "Issue #$ISSUE_NUMBER not found"
  exit 1
fi
```

### Label Already Exists

```bash
# Safe to add label multiple times (idempotent)
gh issue edit ISSUE_NUMBER --add-label "priority: high"
```

### Permission Errors

```bash
# Check if authenticated
if ! gh auth status &>/dev/null; then
  echo "Error: Not authenticated with GitHub CLI"
  exit 1
fi
```

## Dry Run Mode

When operating in dry-run mode, output intended actions without executing:

```bash
DRY_RUN=true

if [ "$DRY_RUN" = true ]; then
  echo "[DRY RUN] Would add label: priority: high"
  echo "[DRY RUN] Would post comment with analysis"
else
  gh issue edit ISSUE_NUMBER --add-label "priority: high"
  gh issue comment ISSUE_NUMBER --body "$ANALYSIS"
fi
```

## Complete Example Workflow

```bash
#!/bin/bash
set -e

ISSUE_NUMBER=$1
DRY_RUN=${2:-false}

# 1. Fetch issue data
echo "Fetching issue #$ISSUE_NUMBER..."
ISSUE_DATA=$(gh issue view $ISSUE_NUMBER --json title,body,labels,reactions,comments)

# 2. Extract tool information
TOOL_URL=$(echo "$ISSUE_DATA" | jq -r '.body' | grep -o 'https://github.com/[^/]*/[^/ ]*')

# 3. Fetch tool metrics
if [ -n "$TOOL_URL" ]; then
  REPO=${TOOL_URL#https://github.com/}
  STARS=$(gh api repos/$REPO --jq '.stargazers_count')
  echo "Tool has $STARS stars"
fi

# 4. Calculate score
# [Run scoring logic - see scoring-criteria.md]

# 5. Add labels
if [ "$DRY_RUN" != "true" ]; then
  gh issue edit $ISSUE_NUMBER --add-label "priority: high"
  echo "Added priority label"
else
  echo "[DRY RUN] Would add label: priority: high"
fi

# 6. Post analysis
if [ "$DRY_RUN" != "true" ]; then
  gh issue comment $ISSUE_NUMBER --body "$ANALYSIS"
  echo "Posted analysis comment"
else
  echo "[DRY RUN] Would post analysis"
fi
```

## Testing

Test GitHub operations in dry-run mode:

```bash
# Test fetching issue
gh issue view 1234 --json title,body,labels

# Test label operations (on test issue)
gh issue edit TEST_ISSUE --add-label "priority: high"
gh issue edit TEST_ISSUE --remove-label "priority: high"

# Test comment formatting
echo "$ANALYSIS" | gh issue comment TEST_ISSUE --body-file -
```

## Workflow Integration

The GitHub Action workflow (`.github/workflows/module-requests.yml`) provides:

```yaml
env:
  GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}

steps:
  - uses: anthropics/claude-code-action@v1
    with:
      prompt: |
        Use the `triaging-module-requests` skill
        Mode: ${{ mode }}
        Issue: #${{ issue_number }}
```

This ensures Claude has access to `gh` CLI with appropriate permissions.
