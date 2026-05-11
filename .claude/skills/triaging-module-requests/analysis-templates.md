# Analysis Comment Template

This document provides a unified template for analysis comments posted to module request issues.

## Main Template

Keep all language succinct, technical and to the point. Avoid being overly effusive.

```markdown
Thanks for requesting a new MultiQC module! This is an automated triage review to help prioritise development work.

| Item               | Details                                                                                    |
| ------------------ | ------------------------------------------------------------------------------------------ |
| **Tool**           | [Tool Name]                                                                                |
| **Repository**     | [GitHub URL] (⭐ [Stars if available])                                                     |
| **Priority Score** | [Score]/100 [Priority Band: 🔴 High ≥70 \| 🟡 Medium 40-69 \| 🟢 Low 20-39 \| ⚪ Hold <20] |

<details>

### Score Breakdown

| Category                 | Score | Notes        |
| ------------------------ | ----- | ------------ |
| 🌟 Tool Popularity       | XX/25 | [Brief note] |
| 📦 Package Downloads     | XX/15 | [Brief note] |
| 💬 Community Engagement  | XX/35 | [Brief note] |
| ✅ Request Quality       | XX/20 | [Brief note] |
| ⚙️ Technical Feasibility | XX/15 | [Brief note] |

**Total: [Score]/100**

### [Priority Band]: [Status Message]

[Customized feedback paragraph, mentioning any missing information from the issue]

**[Strengths/What's great:]**

- ✨ [Positive point 1]
- ✨ [Positive point 2]

**[Areas for improvement/Required actions:]**

- [Icon: 💡 for suggestions, ⚠️ for important, 🚫 for blocking] [Improvement 1] [(+X points) if applicable]
- [Icon] [Improvement 2]

### Next Steps

[Customized action items based on priority:

- High: No action needed, follow issue for updates
- Medium: Suggestions for improvement, re-evaluation process
- Low: How to add missing information and increase priority
- Hold: Required actions to move forward with examples]

[Closing sentiment: appreciation/encouragement]

---

<sup>
This analysis was performed automatically. Comment <code>@claude analyze-module</code> after making improvements for re-evaluation.
Learn more: <a href="https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md">Module Triage System Guide</a>
</sup>

</details>
```

## Customization Guidelines

When using this template, adapt these sections based on priority:

### High Priority (≥70)

- **Status:** "Ready for Development" or "Excellent Request"
- **Tone:** Enthusiastic, welcoming, appreciative
- **Focus:** Acknowledge strengths, thank requester
- **Next Steps:** Added to development queue, follow for updates

### Medium Priority (40-69)

- **Status:** "Good Candidate" or "Solid Request"
- **Tone:** Positive, helpful, collaborative
- **Focus:** Balance strengths with improvement suggestions
- **Next Steps:** How improvements can increase priority, re-evaluation process

### Low Priority (20-39)

- **Status:** "Needs Improvement" or "Has Potential"
- **Tone:** Encouraging, educational, supportive
- **Focus:** What's missing, how to improve
- **Next Steps:** Specific actions to provide missing info, re-evaluation instructions

### Hold Priority (<20)

- **Status:** "Critical Information Needed" or "On Hold"
- **Tone:** Patient, instructive, clear
- **Focus:** Blocking issues, required actions with examples
- **Next Steps:** Step-by-step instructions, offer help

## Re-Analysis Template

For subsequent analyses after improvements:

```markdown
## 🔄 Re-Analysis Results

**Previous Score:** [Old Score]/100 ([Old Priority])
**Current Score:** [New Score]/100 ([New Priority]) [Trend: ↗️ Improved | ↘️ Decreased | → Unchanged]

### What Changed

- [Improvement 1]: +X points
- [Improvement 2]: +X points
- [Change 3]: ±X points

### Progress

[Specific praise for improvements made]
[Remaining suggestions if not yet high priority]
```

## Special Case Templates

Add these sections when applicable:

### Missing Example Files

```markdown
### 📁 Missing Example Files

**Important:** Example files are the most critical component of a module request (+8 points).

**Please provide example files in one of these ways:**

- ✅ Drag and drop actual files from the tool into this issue
- ✅ Use `.zip` if GitHub doesn't support the file type
- ✅ Submit a PR to the [MultiQC/test-data](https://github.com/MultiQC/test-data) repository
- ✅ Link to files already in the test-data repository

**Requirements:**

- Include 2-3 representative examples (typical output)
- ❌ Do NOT copy/paste file contents (formatting and whitespace matter!)

Without example files, module development is very difficult. This is the #1 way to increase your request's priority.
```

### Tool Not Found

```markdown
### ⚠️ Tool Repository Not Found

The repository URL provided could not be accessed. This may be because:

- The repository is private
- The URL is incorrect
- The repository has been deleted

**Action needed:** Please provide a valid, public repository URL or tool homepage. This is required for assessing tool popularity and maintenance status.
```

### Related Requests

```markdown
### 🔄 Related Requests Found

This request is related to:

- #[ISSUE_NUMBER]: [Similar tool/request]

Consider coordinating with these requesters to combine efforts, increase community engagement, and share use cases. Related requests increase this request's priority (+5 points each).
```

### Already Implemented

```markdown
### ✅ Module May Already Exist

Our records show that [Tool Name] might already have a MultiQC module:

- Module: `multiqc.modules.[module_name]`
- Documentation: [Link]

Please check the [list of existing modules](https://docs.seqera.io/multiqc/modules/) before proceeding. If the existing module doesn't meet your needs, please explain what additional functionality you require.
```

## General Guidelines

When using this template:

1. **Be specific:** Replace all [placeholders] with actual data
2. **Be encouraging:** Always acknowledge what's good first
3. **Be actionable:** Provide concrete next steps, not vague suggestions
4. **Be concise:** Remove sections that aren't relevant
5. **Be accurate:** Double-check all score calculations
6. **Be helpful:** Anticipate questions and provide guidance

Always end on a constructive note and invite questions.
