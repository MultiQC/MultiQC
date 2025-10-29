# Analysis Comment Templates

This document provides templates for analysis comments posted to module request issues. Customize based on priority score and specific circumstances.

## General Template Structure

```markdown
## 📊 Module Request Analysis

**Tool:** [Tool Name]
**Repository:** [GitHub URL]
**Priority Score:** [Score]/100 [Priority Band Emoji and Label]

### Score Breakdown

| Category              | Score | Max          | Notes |
| --------------------- | ----- | ------------ | ----- |
| Tool Popularity       | XX/25 | [Brief note] |
| Package Downloads     | XX/15 | [Brief note] |
| Community Engagement  | XX/35 | [Brief note] |
| Request Quality       | XX/20 | [Brief note] |
| Technical Feasibility | XX/15 | [Brief note] |

### [Priority Band]: Recommendation

[Specific feedback paragraph]

### Next Steps

[Specific actionable items]

---

_This analysis was performed automatically. For questions or to request re-analysis, comment `@claude analyze-module`._
```

## High Priority Template (≥70)

```markdown
## 📊 Module Request Analysis

**Tool:** [Tool Name]
**Repository:** [GitHub URL] (⭐ [Stars])
**Priority Score:** [Score]/100 🔴 **High Priority**

### Score Breakdown

| Category                 | Score | Notes                             |
| ------------------------ | ----- | --------------------------------- |
| 🌟 Tool Popularity       | XX/25 | [Stars count, active maintenance] |
| 📦 Package Downloads     | XX/15 | [Download count/month]            |
| 💬 Community Engagement  | XX/35 | [Reactions, comments, duplicates] |
| ✅ Request Quality       | XX/20 | [Completeness, examples]          |
| ⚙️ Technical Feasibility | XX/15 | [Output format, metrics clarity]  |

**Total: [Score]/100**

### 🔴 High Priority: Ready for Development

This is an excellent module request! The tool is popular, actively maintained, and your request is complete with example files. This module would add significant value to MultiQC.

**What's great:**

- ✨ [Specific positive point 1]
- ✨ [Specific positive point 2]
- ✨ [Specific positive point 3]

### Next Steps

1. This request has been added to the **High Priority** development queue
2. [Any specific next action: maintainer review, contributor assignment, etc.]
3. Follow issue updates for progress

Thank you for this well-prepared request!

---

_This analysis was performed automatically. For questions, comment `@claude analyze-module`._
_Learn more: [Module Triage System Guide](https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md)_
```

## Medium Priority Template (40-69)

```markdown
## 📊 Module Request Analysis

**Tool:** [Tool Name]
**Repository:** [GitHub URL] (⭐ [Stars])
**Priority Score:** [Score]/100 🟡 **Medium Priority**

### Score Breakdown

| Category                 | Score | Notes        |
| ------------------------ | ----- | ------------ |
| 🌟 Tool Popularity       | XX/25 | [Brief note] |
| 📦 Package Downloads     | XX/15 | [Brief note] |
| 💬 Community Engagement  | XX/35 | [Brief note] |
| ✅ Request Quality       | XX/20 | [Brief note] |
| ⚙️ Technical Feasibility | XX/15 | [Brief note] |

**Total: [Score]/100**

### 🟡 Medium Priority: Good Candidate

This is a solid module request with [positive aspect]. With some improvements, this could move to high priority.

**Strengths:**

- ✨ [What's already good]
- ✨ [Another strong point]

**Quick ways to increase priority:**

- 💡 [Specific improvement 1] (+X points)
- 💡 [Specific improvement 2] (+X points)
- 💡 [Specific improvement 3] (+X points)

### Next Steps

1. Consider making the suggested improvements above
2. Share this request with colleagues who use [Tool] to gather more community support (each 👍 adds 1 point)
3. The request will be re-evaluated during weekly triage

---

_This analysis was performed automatically. Comment `@claude analyze-module` after making improvements for re-evaluation._
_Learn more: [Module Triage System Guide](https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md)_
```

## Low Priority Template (20-39)

```markdown
## 📊 Module Request Analysis

**Tool:** [Tool Name]
**Repository:** [GitHub URL if available]
**Priority Score:** [Score]/100 🟢 **Low Priority**

### Score Breakdown

| Category                 | Score | Notes        |
| ------------------------ | ----- | ------------ |
| 🌟 Tool Popularity       | XX/25 | [Brief note] |
| 📦 Package Downloads     | XX/15 | [Brief note] |
| 💬 Community Engagement  | XX/35 | [Brief note] |
| ✅ Request Quality       | XX/20 | [Brief note] |
| ⚙️ Technical Feasibility | XX/15 | [Brief note] |

**Total: [Score]/100**

### 🟢 Low Priority: Needs Improvement

This request has potential but needs more information or community interest to prioritize for development.

**Areas for improvement:**

- ⚠️ [Critical missing item 1]
- ⚠️ [Critical missing item 2]
- ℹ️ [Optional improvement]

**How to increase priority:**

The most impactful improvements:

1. 📁 **Upload example files** (not pasted text) → +8 points
2. 🔗 **Provide tool homepage/repository** → +3 points (+ enables popularity scoring)
3. 👍 **Gather community support** → +1 point per reaction

### Next Steps

1. Add the missing information described above
2. Edit your original request (don't add comments) to include example files
3. Comment `@claude analyze-module` when ready for re-evaluation

We're here to help! Feel free to ask questions if you need clarification.

---

_This analysis was performed automatically. Comment `@claude analyze-module` after making improvements._
_Learn more: [Module Triage System Guide](https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md)_
```

## Hold Priority Template (<20)

```markdown
## 📊 Module Request Analysis

**Tool:** [Tool Name if provided]
**Priority Score:** [Score]/100 ⚪ **On Hold**

### Score Breakdown

| Category                 | Score | Notes   |
| ------------------------ | ----- | ------- |
| 🌟 Tool Popularity       | XX/25 | [Issue] |
| 📦 Package Downloads     | XX/15 | [Issue] |
| 💬 Community Engagement  | XX/35 | [Issue] |
| ✅ Request Quality       | XX/20 | [Issue] |
| ⚙️ Technical Feasibility | XX/15 | [Issue] |

**Total: [Score]/100**

### ⚪ On Hold: Critical Information Needed

This request cannot be prioritized yet due to missing critical information.

**Required actions:**

- 🚫 [Blocking issue 1]
- 🚫 [Blocking issue 2]

### How to Move Forward

To get this request prioritized, you need to:

1. **[Most critical item]** - This is essential for any module development

   - How: [Specific instructions]
   - Why: [Brief explanation]

2. **[Second critical item]**

   - How: [Specific instructions]

3. **[Third item if applicable]**

### Example of a High-Quality Request

Check out [link to example request] to see what a complete, high-priority request looks like.

### Need Help?

If you're unsure how to provide this information:

- See the [Contributing Guide](https://github.com/MultiQC/MultiQC/blob/main/.github/CONTRIBUTING.md)
- Review the [Module Triage System Guide](https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md)
- Ask questions here and we'll guide you!

---

_This analysis was performed automatically. Comment `@claude analyze-module` after adding the required information._
```

## Special Case Templates

### Missing Example Files

Add this section to any priority level:

```markdown
### 📁 Missing Example Files

**Important:** Example files are the most critical component of a module request (+8 points).

**Please:**

- ✅ Drag and drop actual files from the tool into this issue
- ✅ Use `.zip` if GitHub doesn't support the file type
- ✅ Include 2-3 representative examples (typical output)
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

### Duplicate Request

```markdown
### 🔄 Related Requests Found

This request is related to:

- #[ISSUE_NUMBER]: [Similar tool/request]

Consider coordinating with these requesters to:

- Combine efforts and examples
- Increase community engagement
- Share use cases

Related requests increase this request's priority (+5 points each).
```

### Already Implemented

```markdown
### ✅ Module May Already Exist

Our records show that [Tool Name] might already have a MultiQC module:

- Module: `multiqc.modules.[module_name]`
- Documentation: [Link]

Please check the [list of existing modules](https://multiqc.info/modules/) before proceeding. If the existing module doesn't meet your needs, please explain what additional functionality you require.
```

## Re-Analysis Template

For subsequent analyses (e.g., after improvements):

```markdown
## 🔄 Re-Analysis Results

**Previous Score:** [Old Score]/100 ([Old Priority])
**Current Score:** [New Score]/100 ([New Priority]) [Trend: ↗️/↘️/→]

### What Changed

- [Improvement 1]: +X points
- [Improvement 2]: +X points
- [Change 3]: ±X points

[Rest of standard template for current priority level]

### Progress

Great work! [Specific praise for improvements made]
[Remaining suggestions if not yet high priority]
```

## Customization Guidelines

When using these templates:

1. **Be specific:** Replace all [placeholders] with actual data
2. **Be encouraging:** Always acknowledge what's good first
3. **Be actionable:** Provide concrete next steps, not vague suggestions
4. **Be concise:** Remove sections that aren't relevant
5. **Be accurate:** Double-check all score calculations
6. **Be helpful:** Anticipate questions and provide guidance

### Tone Guidelines

- **High Priority:** Enthusiastic, welcoming, appreciative
- **Medium Priority:** Positive, helpful, collaborative
- **Low Priority:** Encouraging, educational, supportive
- **Hold:** Patient, instructive, clear about requirements

Always end on a constructive note and invite questions.
