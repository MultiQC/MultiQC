---
name: Triaging MultiQC Module Requests
description: Analyzes MultiQC module requests from GitHub issues, evaluates tool popularity via GitHub/PyPI/Conda APIs, assesses community engagement, and generates priority scores with recommendations. Use when analyzing issues labeled 'module: new' or when asked to prioritize module requests.
---

# Triaging MultiQC Module Requests

This skill helps evaluate and prioritize MultiQC module requests systematically.

## When to Use

- Analyzing GitHub issues with `module: new` label
- Responding to `@claude analyze-module` comments
- Bulk triage of pending module requests
- Prioritizing development roadmap

## Evaluation Process

### 1. Extract Request Information

From the issue, identify:
- **Tool name**: Primary identifier
- **Tool homepage**: GitHub/GitLab/Bitbucket repository or project page
- **Description**: What the tool does
- **Example files**: Uploaded log/output files (not copy-pasted text)
- **Use case**: Why user needs this module
- **Community demand**: Reactions, comments, duplicates

### 2. Evaluate Tool Popularity

**GitHub Repository Metrics** (if available):
```python
# Use GitHub API
api_url = f"https://api.github.com/repos/{owner}/{repo}"
metrics = {
    "stars": repo["stargazers_count"],
    "forks": repo["forks_count"],
    "recent_activity": repo["pushed_at"],
    "open_issues": repo["open_issues_count"]
}
```

Scoring:
- ⭐⭐⭐⭐⭐ 1000+ stars (25 points)
- ⭐⭐⭐⭐ 500-999 stars (20 points)
- ⭐⭐⭐ 100-499 stars (15 points)
- ⭐⭐ 50-99 stars (10 points)
- ⭐ <50 stars (5 points)

**Package Manager Metrics**:
- Check PyPI: `https://pypi.org/pypi/{package}/json`
- Check Conda: `https://api.anaconda.org/package/bioconda/{package}`
- Check Bioconductor: For R packages

Download metrics (if available):
- High usage (15 points): >10K downloads/month
- Medium usage (10 points): 1K-10K downloads/month
- Low usage (5 points): <1K downloads/month

### 3. Assess Community Engagement

**Issue Metrics**:
- 👍 reactions on issue (1 point each, max 10)
- Comments from unique users (2 points each, max 10)
- Linked duplicate issues (5 points each, max 15)
- Author provides example files (5 points)

### 4. Evaluate Request Quality

**Completeness** (max 20 points):
- [ ] Tool name provided (2 points)
- [ ] Tool homepage/repo URL (3 points)
- [ ] Clear description (2 points)
- [ ] Example output files uploaded (8 points)
- [ ] File format described (2 points)
- [ ] Expected visualizations suggested (3 points)

**Technical Feasibility** (max 15 points):
- Example files parseable (10 points)
- Clear metrics identified (5 points)

### 5. Calculate Priority Score

Total possible: 100 points

```
Total Score =
  Tool Popularity (0-25) +
  Package Downloads (0-15) +
  Community Engagement (0-35) +
  Request Quality (0-20) +
  Technical Feasibility (0-15)
```

**Priority Bands**:
- 🔴 **High Priority** (70-100): Urgent, high-value module
- 🟡 **Medium Priority** (40-69): Valuable, queue for development
- 🟢 **Low Priority** (20-39): Consider if time permits
- ⚪ **Hold** (<20): Needs more information or community interest

## Output Format

Post analysis as GitHub comment:

```markdown
## 📊 Module Request Analysis

### Tool Information
- **Name**: {tool_name}
- **Repository**: {repo_url}
- **Stars**: ⭐ {stars}
- **Downloads**: {downloads}/month

### Priority Score: {total}/100
- Tool Popularity: {popularity_score}/25
- Package Downloads: {downloads_score}/15
- Community Engagement: {engagement_score}/35
- Request Quality: {quality_score}/20
- Technical Feasibility: {feasibility_score}/15

### Priority Band: {🔴 High / 🟡 Medium / 🟢 Low / ⚪ Hold}

### Recommendation

{1-2 paragraph recommendation}

### Next Steps

{Specific actions:
- For high priority: assign to developer, estimate timeline
- For medium priority: add to roadmap
- For low priority: label and defer
- For hold: request more information}

---
*🤖 Automated analysis via Claude Code*
```

## Bulk Triage Mode

When processing multiple issues:

1. Fetch all open issues with `module: new` label
2. Filter out already-triaged (have priority label)
3. Analyze each issue
4. Sort by priority score
5. Generate summary report:

```markdown
## Weekly Module Request Triage

Analyzed {N} pending requests:

### High Priority ({count})
- #{issue_num} {tool_name} - Score: {score}/100
- ...

### Medium Priority ({count})
- ...

### Low Priority ({count})
- ...

### Needs Info ({count})
- ...
```

## Label Recommendations

Based on priority score, suggest labels:
- `priority: high` (score ≥ 70)
- `priority: medium` (score 40-69)
- `priority: low` (score 20-39)
- `needs: more info` (score < 20 or missing critical info)
- `needs: example files` (no uploads)

## Rate Limiting

When using APIs:
- GitHub: 60 requests/hour (unauthenticated), 5000/hour (authenticated)
- PyPI: No strict limits, be respectful
- Conda: No strict limits

Cache results to avoid re-fetching.

## Related Resources

- [Scoring Rubric Details](./scoring-rubric.md)
- [GitHub API Documentation](https://docs.github.com/rest)
- [MultiQC Module Guidelines](../../CLAUDE.md)
