# MultiQC Module Request Triage System

This guide explains how MultiQC manages and prioritizes module requests through an automated triage system.

## Overview

The triage system helps manage the growing number of module requests by automatically analyzing and prioritizing them based on objective criteria. This ensures the most valuable modules are developed first while giving contributors clear guidance on how to create successful requests.

## How It Works

### Automated Analysis

When you submit a module request or comment `@claude analyze-module` on an existing request, the system automatically:

1. **Evaluates tool popularity** via GitHub stars, package downloads, and community metrics
2. **Assesses request quality** based on completeness and example files
3. **Calculates a priority score** (0-100) using multiple weighted factors
4. **Assigns priority labels** and moves the issue through the project board
5. **Provides feedback** with specific recommendations for improvement

### Priority Scoring

Your request receives a score based on:

| Category                  | Max Points | What's Evaluated                              |
| ------------------------- | ---------- | --------------------------------------------- |
| **Tool Popularity**       | 25         | GitHub stars, forks, recent activity          |
| **Package Downloads**     | 15         | PyPI, Conda, or Bioconductor download metrics |
| **Community Engagement**  | 35         | 👍 reactions, comments, duplicate requests    |
| **Request Quality**       | 20         | Complete information, example files           |
| **Technical Feasibility** | 15         | Parseable output, clear metrics               |

**Priority Bands:**

- 🔴 **High Priority** (70-100): Urgent, high-value module
- 🟡 **Medium Priority** (40-69): Valuable, queued for development
- 🟢 **Low Priority** (20-39): Consider if time permits
- ⚪ **Hold** (<20): Needs more information or community interest

## Getting Your Request Prioritized

### Quick Wins

The fastest ways to increase your priority score:

1. **Upload example files** (+8 points) - Most critical factor!
   - Drag and drop actual files (don't copy/paste)
   - Include typical tool output
   - Use `.zip` for unsupported file types

2. **Choose popular tools** (+25 points)
   - Tools with >1000 GitHub stars get maximum points
   - Active maintenance and recent releases help

3. **Provide complete information** (+20 points)
   - Tool name and homepage
   - Clear description
   - Expected visualizations
   - Use case explanation

4. **Generate community interest** (+35 points)
   - Each 👍 reaction adds 1 point
   - Meaningful comments add 2 points each
   - Related requests add 5 points each

### Interactive Help

Get instant feedback on your request:

```
@claude analyze-module
```

This triggers a detailed analysis showing:

- Current priority score breakdown
- Specific recommendations for improvement
- Comparison to similar requests
- Next steps to increase priority

## Project Board Workflow

Module requests flow through a structured project board:

| Column                    | Description                  | How to Get Here                    |
| ------------------------- | ---------------------------- | ---------------------------------- |
| **Needs Analysis**        | New requests awaiting triage | Automatic when issue created       |
| **Needs Examples**        | Missing example files        | Any request without uploaded files |
| **Low Priority**          | Score < 40                   | Incomplete or niche requests       |
| **Medium Priority**       | Score 40-69                  | Good candidates needing refinement |
| **High Priority**         | Score ≥ 70                   | Ready for development              |
| **Ready for Development** | Maintainer-approved          | Manual promotion after review      |
| **In Development**        | Active work                  | When implementation begins         |
| **Complete**              | Module merged                | When PR is merged                  |

## Automated Actions

### New Request Analysis

When you open a module request:

- Instant analysis and priority scoring
- Automatic label assignment
- Project board placement
- Initial feedback comment

### Weekly Bulk Triage

Every Monday at 9 AM UTC:

- All open module requests are re-analyzed
- Priority scores are updated
- Stale requests are identified
- Board positions are synchronized

### On-Demand Analysis

Comment `@claude analyze-module` anytime to:

- Get current priority score
- Receive specific recommendations
- Understand blocking issues
- See comparative analysis

## Request Quality Guidelines

### What Makes a Great Request

✅ **Do:**

- Upload actual example files from the tool
- Provide the tool's GitHub/GitLab repository URL
- Explain why this module would be valuable
- Suggest which metrics should appear in general stats
- Describe expected visualizations

❌ **Don't:**

- Copy/paste file contents (formatting matters!)
- Request modules for unpopular or unmaintained tools
- Submit duplicate requests (search first!)
- Leave fields blank or incomplete

### Example Files

Example files are the most critical component:

- **Quality over quantity**: 2-3 representative examples beat 20 minimal ones
- **Real data**: Use actual tool output, not synthetic examples
- **Varied cases**: Include both typical and edge cases if possible
- **Proper format**: Upload files, don't paste text (whitespace matters!)
- **Reasonable size**: Truncate large files but keep format clear

## Understanding Analysis Results

When your request is analyzed, you'll receive a comment like:

```markdown
Thanks for requesting a new MultiQC module! This is an automated triage review to help prioritise development work.

| Item               | Details                                      |
| ------------------ | -------------------------------------------- |
| **Tool**           | FastQC                                       |
| **Repository**     | https://github.com/s-andrews/FastQC (⭐ 450) |
| **Priority Score** | 75/100 🔴 **High Priority**                  |

<details>

### Score Breakdown

| Category                 | Score | Notes                           |
| ------------------------ | ----- | ------------------------------- |
| 🌟 Tool Popularity       | 20/25 | 450 stars, actively maintained  |
| 📦 Package Downloads     | 12/15 | 50K downloads/month on Bioconda |
| 💬 Community Engagement  | 18/35 | 8 👍 reactions, 3 comments      |
| ✅ Request Quality       | 20/20 | Complete info, example files    |
| ⚙️ Technical Feasibility | 5/15  | Complex output format           |

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
```

## For Contributors

### Improving Existing Requests

If your request has low priority, you can improve it:

1. **Add example files** if missing (biggest impact!)
2. **Explain use case** - why is this module valuable?
3. **Engage community** - share with colleagues who use the tool
4. **Provide context** - link to publications, workflows, or pipelines

### Contributing Development

High-priority requests are great for contributors:

- Clear requirements and examples
- Community-vetted value
- Maintainer support available
- Good for portfolio building

See [CONTRIBUTING.md](../../.github/CONTRIBUTING.md) for development guidelines.

## For Maintainers

### Manual Overrides

Maintainers can override automated priority:

1. **Manually add priority labels** to adjust board placement
2. **Move issues** to "Ready for Development" when verified
3. **Close stale requests** with kind explanation
4. **Pin important requests** for visibility

### System Configuration

The triage system is implemented through:

- **Skills**: `.claude/skills/triaging-module-requests/` - Core triage logic
- **Workflows**: `.github/workflows/module-requests.yml` - Automation triggers
- **Commands**: `.claude/commands/new-module.md` - Module generation

For project board setup, see [module-triage-project-setup.md](./module-triage-project-setup.md).

## Benefits

This system provides:

- **Transparency**: Clear criteria for prioritization
- **Fairness**: Objective, consistent evaluation
- **Efficiency**: Automated processing of high volume
- **Guidance**: Actionable feedback for contributors
- **Tracking**: Visual workflow through project board

## Feedback and Improvements

The triage system is continuously evolving. If you have suggestions:

- Open a discussion in the [Seqera Community Forum](https://community.seqera.io/c/multiqc/6)
- Comment on the [triage system tracking issue](https://github.com/MultiQC/MultiQC/issues/3219)
- Propose changes via pull request

## Related Resources

- [Module Creation Guide](../../.claude/skills/creating-multiqc-modules/SKILL.md)
- [Contributing Guidelines](../../.github/CONTRIBUTING.md)
- [Project Board Setup](./module-triage-project-setup.md)
- [Example Module Requests](https://github.com/MultiQC/MultiQC/issues?q=label%3A%22module%3A+new%22+label%3A%22priority%3A+high%22)
