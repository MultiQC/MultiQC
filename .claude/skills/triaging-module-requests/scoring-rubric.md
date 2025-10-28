# Module Request Scoring Rubric

Detailed breakdown of the scoring system for MultiQC module request prioritization.

## Category 1: Tool Popularity (0-25 points)

### GitHub Stars

- **25 points**: 1000+ stars (widely adopted, battle-tested)
- **20 points**: 500-999 stars (popular in community)
- **15 points**: 100-499 stars (established tool)
- **10 points**: 50-99 stars (emerging tool)
- **5 points**: <50 stars (niche or new tool)
- **0 points**: Not on GitHub or no public repo

### Additional Popularity Indicators

- Part of major pipeline/framework (+3): nf-core, Snakemake, CWL
- Cited in recent papers (+2): Published within 2 years
- Actively maintained (+2): Commits within 3 months

## Category 2: Package Downloads (0-15 points)

### Monthly Download Metrics

**Bioconda** (primary source for bioinformatics):

- **15 points**: >10,000 downloads/month
- **12 points**: 5,000-10,000 downloads/month
- **10 points**: 1,000-5,000 downloads/month
- **7 points**: 500-1,000 downloads/month
- **5 points**: <500 downloads/month

**PyPI** (Python packages):

- Similar thresholds as Bioconda
- Use if not available on Bioconda

**NPM** (JavaScript/Node tools):

- Adjust thresholds ×10 (higher volume)

**Other Indicators**:

- Docker pulls (+2): >1000 pulls
- Singularity usage (+1): Available in registries

## Category 3: Community Engagement (0-35 points)

### Issue Reactions (0-10 points)

- 1 point per 👍 reaction (max 10)
- 0.5 points per other positive reaction

### Community Discussion (0-10 points)

- 2 points per unique commenter (max 10)
- Comments showing:
  - Similar use case
  - Request for same tool
  - Offer to help test

### Linked Issues (0-15 points)

- 5 points per duplicate request (max 15)
- Shows repeated independent demand

## Category 4: Request Quality (0-20 points)

### Information Completeness (0-12 points)

**Required Fields**:

- Tool name (2 points)
- Tool homepage/repository (3 points)
- Clear tool description (2 points)
- Use case explanation (2 points)
- File format description (3 points)

### Example Files (0-8 points)

**Uploaded Files** (not copy-paste):

- 8 points: 3+ diverse example files
- 6 points: 2 example files
- 4 points: 1 example file
- 0 points: No uploads or only pasted snippets

**File Quality**:

- Real tool output (not synthetic)
- Representative of typical use
- Different scenarios/parameters

## Category 5: Technical Feasibility (0-15 points)

### Data Parseability (0-10 points)

**File Format**:

- 10 points: Structured format (JSON, TSV, YAML)
- 7 points: Semi-structured (key-value pairs)
- 4 points: Unstructured but parseable (regex patterns)
- 0 points: Binary or proprietary format

**Parsing Complexity**:

- Simple (+2): Single-pass parsing
- Moderate (+1): Multi-section parsing
- Complex (+0): State machine or complex logic

### Metrics Clarity (0-5 points)

- 5 points: Clear numeric metrics identified
- 3 points: Metrics extractable with analysis
- 1 point: Unclear what to extract
- 0 points: No useful metrics apparent

## Priority Band Calculation

```python
def calculate_priority_band(total_score):
    if total_score >= 70:
        return "🔴 High Priority"
    elif total_score >= 40:
        return "🟡 Medium Priority"
    elif total_score >= 20:
        return "🟢 Low Priority"
    else:
        return "⚪ Hold"
```

### High Priority (70-100)

**Characteristics**:

- Popular, widely-used tool
- Strong community demand
- Complete, high-quality request
- Clear implementation path

**Action**: Prioritize for next development cycle

### Medium Priority (40-69)

**Characteristics**:

- Moderately popular tool OR
- High quality request for niche tool OR
- Strong community demand for emerging tool

**Action**: Add to roadmap, implement when capacity allows

### Low Priority (20-39)

**Characteristics**:

- Niche tool with limited usage OR
- Incomplete request for popular tool OR
- Feasible but low demand

**Action**: Keep open, revisit if demand increases

### Hold (<20)

**Characteristics**:

- Insufficient information
- No example files
- Unclear use case
- Very low interest

**Action**: Request more information before proceeding

## Special Cases

### Override Criteria

**Immediate High Priority** (regardless of score):

- Requested by major MultiQC contributor
- Part of funded collaboration
- Critical for major pipeline integration

**Automatic Defer**:

- Tool deprecated or abandoned
- Functionality already covered by existing module
- License incompatibility

### Tie-Breaking

When scores are equal:

1. Request with example files wins
2. More recent tool version wins
3. Better-maintained tool wins
4. First requested wins

## Example Scoring

### Example 1: High Priority Tool

**samtools ampliconstats**

- GitHub stars: 1500+ (25 pts)
- Conda downloads: 50K/month (15 pts)
- Issue reactions: 8 👍 (8 pts)
- Comments: 4 users (8 pts)
- Complete info + examples (20 pts)
- Clear TSV format (10 pts)

**Total: 86/100** → 🔴 High Priority

### Example 2: Medium Priority Tool

**obscure-aligner**

- GitHub stars: 120 (15 pts)
- Conda downloads: 2K/month (10 pts)
- Issue reactions: 2 👍 (2 pts)
- Comments: 1 user (2 pts)
- Complete info + examples (20 pts)
- Complex log format (6 pts)

**Total: 55/100** → 🟡 Medium Priority

### Example 3: Needs Info

**new-tool-request**

- GitHub stars: 200 (15 pts)
- No package metrics (0 pts)
- Issue reactions: 0 (0 pts)
- No comments (0 pts)
- Missing examples (8 pts)
- Unknown feasibility (0 pts)

**Total: 23/100** → 🟢 Low Priority (but mark "needs: example files")

## Calibration Notes

Scoring should be:

- **Objective**: Based on measurable metrics
- **Reproducible**: Same inputs → same score
- **Transparent**: Clear rationale for each point
- **Adjustable**: Weights can be tuned based on outcomes

Review scoring effectiveness quarterly and adjust weights if needed.
