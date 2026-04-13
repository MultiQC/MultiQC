# Scoring Criteria for Module Request Triage

This document provides the detailed rubric for calculating priority scores (0-100) for MultiQC module requests.

## Score Categories

### 1. Tool Popularity (25 points max)

**GitHub Stars** (primary metric):

- ≥200 stars: 25 points
- 100-199 stars: 20 points
- 50-99 stars: 15 points
- 25-49 stars: 12 points
- 10-24 stars: 8 points
- 5-9 stars: 5 points
- <5 stars: 2 points

**Bonus factors** (add up to 3 points):

- Active maintenance: +2 (commit in last 3 months)
- High fork ratio: +1 (forks/stars > 0.1)

**Example calculation:**

- Tool with 75 stars, recent commit, 10 forks: 15 + 2 + 1 = 18/25

### 2. Package Downloads (15 points max)

Check in order of preference: PyPI → Conda → Bioconda → Docker pulls

**Monthly Downloads:**

- ≥1M downloads: 15 points
- 500K-999K: 13 points
- 100K-499K: 11 points
- 50K-99K: 9 points
- 10K-49K: 7 points
- 1K-9K: 5 points
- 100-999: 3 points
- <100: 1 point

**How to check:**

```bash
# PyPI stats (use scripts/fetch-tool-metrics.js)
node scripts/fetch-tool-metrics.js pypi PACKAGE_NAME

# Conda downloads (from anaconda.org)
node scripts/fetch-tool-metrics.js conda PACKAGE_NAME
```

**Example calculation:**

- PyPI package with 250K monthly downloads: 11/15

### 3. Community Engagement (35 points max)

**Reactions** (👍 on the issue):

- Each 👍 reaction: +1 point (max 15 points)
- Cap at 15 to prevent gaming

**Comments** (substantive, non-bot):

- Each meaningful comment: +2 points (max 10 points)
- Exclude bot comments and "+1" style comments

**Duplicate/Related Requests:**

- Each related open issue: +5 points (max 10 points)
- Search for similar tool names or use cases

**Example calculation:**

- 8 👍 reactions = 8 points
- 3 meaningful comments = 6 points
- 1 related request = 5 points
- Total: 19/35

### 4. Request Quality (20 points max)

**Required Fields Completed:**

- Tool name provided: 2 points
- Tool homepage/repository URL: 3 points
- Tool description: 2 points
- Data/plots suggestions: 2 points
- General stats suggestions: 2 points

**Example Files:**

- Example files uploaded (not pasted): 8 points
- Example files pasted as text: 4 points
- No example files: 0 points

**Additional Context:**

- Context field completed: 1 point

**Example calculation:**

- All required fields: 11 points
- Uploaded example files: 8 points
- Context provided: 1 point
- Total: 20/20 (perfect!)

### 5. Technical Feasibility (15 points max)

**Output Format** (assess from examples or description):

- Structured format (JSON/TSV/CSV): 8 points
- Semi-structured (key-value pairs): 5 points
- Unstructured/text: 2 points
- Unknown/no examples: 0 points

**Metrics Clarity:**

- Clear quantitative metrics: 4 points
- Qualitative or unclear metrics: 2 points
- No metrics evident: 0 points

**Parsing Complexity:**

- Simple parsing expected: 3 points
- Complex parsing needed: 1 point
- Very difficult/unclear: 0 points

**Example calculation:**

- TSV output format: 8 points
- Clear quantitative metrics: 4 points
- Simple parsing: 3 points
- Total: 15/15

## Score Calculation Process

1. **Fetch tool metadata:**

   ```bash
   # Use the metrics script
   node scripts/fetch-tool-metrics.js github OWNER/REPO
   node scripts/fetch-tool-metrics.js pypistats PACKAGE_NAME
   ```

2. **Parse issue body:**
   - Extract tool name, URL, description
   - Check for example files (attachments vs. pasted)
   - Count filled required fields

3. **Check community metrics:**

   ```bash
   gh issue view ISSUE_NUMBER --json reactions,comments
   gh issue list --label "module: new" --search "TOOL_NAME in:title"
   ```

4. **Calculate each category:**
   - Tool Popularity: Use GitHub API results
   - Package Downloads: Use PyPI/Conda stats
   - Community Engagement: Sum reactions + comments + duplicates
   - Request Quality: Count completed fields + example quality
   - Technical Feasibility: Assess from examples/description

5. **Sum total score (0-100)**

6. **Assign priority band:**
   - ≥70: High Priority 🔴
   - 40-69: Medium Priority 🟡
   - 20-39: Low Priority 🟢
   - <20: Hold ⚪

## Edge Cases

### Tool URL Not Provided

- Assign 0 for Tool Popularity
- Assign 0 for Package Downloads
- Note in feedback: "Please provide tool homepage/repository"

### Private or Deleted Repository

- Assign minimum scores for popularity (2 points)
- Note in analysis comment
- Request alternative URL or more information

### Multiple Package Names

- Use the package with highest downloads
- Note alternative packages in analysis

### Tool Without Package

- Score 0 for Package Downloads
- Don't penalize further; adjust feedback
- GitHub stars become more important

### Example Files in Comments

- Check original issue body and all comments
- Give partial credit (6 points vs. 8) if in comments
- Note: "Consider editing issue to add files to main body"

### Non-Bioinformatics Tools

- Apply same rubric objectively
- If clearly off-topic, note in analysis
- Don't artificially inflate or deflate score

## Scoring Examples

### Example 1: High Priority Request (Score: 88)

- **Tool:** STAR aligner (GitHub: 250 stars, active)
- **Tool Popularity:** 25 + 2 = 27/25 (capped at 25)
- **Downloads:** 150K/month PyPI = 11/15
- **Community:** 12 👍, 4 comments, 1 duplicate = 12 + 8 + 5 = 25/35
- **Quality:** All fields + uploaded files = 20/20
- **Feasibility:** TSV output, clear metrics = 15/15
- **Total:** 96/100 (capped at 25 for popularity) = 88/100 → 🔴 High Priority

### Example 2: Medium Priority Request (Score: 52)

- **Tool:** NewTool (GitHub: 30 stars, recent activity)
- **Tool Popularity:** 12 + 2 = 14/25
- **Downloads:** 5K/month Conda = 5/15
- **Community:** 3 👍, 1 comment, 0 duplicates = 3 + 2 = 5/35
- **Quality:** All fields + uploaded files = 20/20
- **Feasibility:** JSON output, clear metrics = 12/15
- **Total:** 56/100 → 🟡 Medium Priority

### Example 3: Low Priority Request (Score: 28)

- **Tool:** CustomScript (no GitHub, no package)
- **Tool Popularity:** 0/25
- **Downloads:** 0/15
- **Community:** 1 👍, 0 comments = 1/35
- **Quality:** Basic fields, no examples = 11/20
- **Feasibility:** Unknown format = 0/15
- **Total:** 12/100 → ⚪ Hold (needs more info)

## Improving Scores

**Quick wins for requesters:**

1. Upload example files (+8 points)
2. Complete all required fields (+11 points)
3. Share request for 👍 reactions (+1 per reaction)
4. Use popular, well-maintained tools (determined by tool choice)

**Feedback should always:**

- Show current score and breakdown
- Identify missing high-value items
- Provide specific, actionable next steps
- Be encouraging about what's already good

## Validation

After scoring, sanity check:

- Does the priority feel right intuitively?
- Are we comparing apples to apples?
- Is the feedback constructive and specific?
- Would this help a contributor improve their request?
