# MultiQC Claude Code Skills

This directory contains specialized skills that help Claude Code work effectively with MultiQC development, particularly for module request management and module creation.

## Overview

Skills are modular capabilities that Claude automatically invokes based on context. They replace the previous GitHub Actions-based workflows with more flexible, locally-executable capabilities.

## Available Skills

### 1. Triaging Module Requests
**Location**: `triaging-module-requests/`

**Purpose**: Systematically evaluate and prioritize MultiQC module requests from GitHub issues.

**Automatically invoked when**:
- Analyzing issues with `module: new` label
- Responding to `@claude analyze-module` comments
- Asked to prioritize module requests

**Key capabilities**:
- Fetches tool popularity metrics from GitHub, PyPI, Conda
- Evaluates community engagement (reactions, comments, duplicates)
- Assesses request quality and completeness
- Calculates priority scores (0-100 scale)
- Generates structured analysis reports
- Recommends appropriate priority labels

**Supporting files**:
- `scoring-rubric.md` - Detailed breakdown of scoring system

**Example usage**:
```
You: "Can you analyze issue #1234 which is a module request?"
Claude: [Automatically uses triaging skill]
- Fetches GitHub repo stats
- Checks package download metrics
- Evaluates issue quality
- Generates priority score
- Posts analysis comment
```

---

### 2. Creating MultiQC Modules
**Location**: `creating-multiqc-modules/`

**Purpose**: Generate complete, standards-compliant MultiQC module implementations from issue specifications.

**Automatically invoked when**:
- Creating new MultiQC modules
- Converting example log files into code
- Asked to implement module requests
- Need to follow MultiQC coding standards

**Key capabilities**:
- Downloads and analyzes example log files
- Identifies data structures and parsing strategies
- Generates complete module implementation
- Creates comprehensive tests
- Updates configuration files (search_patterns.yaml, pyproject.toml)
- Follows all MultiQC standards and best practices
- Creates pull requests with proper formatting

**Supporting files**:
- `templates/module_template.py` - Production-ready module template
- `examples/simple-module.md` - Example of simple TSV-based module
- `examples/complex-module.md` - Example of complex JSON-based module with multiple plots

**Example usage**:
```
You: "Create a module for the tool in issue #5678"
Claude: [Automatically uses module creation skill]
- Downloads example files from issue
- Analyzes file format and structure
- Generates parsing code
- Creates visualizations
- Writes tests
- Submits PR
```

---

### 3. Analyzing Bioinformatics Log Files
**Location**: `analyzing-log-files/`

**Purpose**: Examine unfamiliar log file formats to understand structure and design parsing strategies.

**Automatically invoked when**:
- Encountering new/unfamiliar log formats
- Designing parsing logic for modules
- Analyzing example files from issues
- Troubleshooting parsing problems

**Key capabilities**:
- Identifies file format (JSON, TSV, key-value, mixed, etc.)
- Detects section structure and markers
- Locates sample name extraction points
- Identifies key metrics and their types
- Suggests appropriate parsing strategies
- Recommends visualization types
- Plans general statistics metrics

**Supporting files**:
- `common-patterns.md` - Reference for frequently encountered formats

**Example usage**:
```
You: "I'm not sure how to parse this log file format"
Claude: [Automatically uses log analysis skill]
- Examines file structure
- Identifies format patterns
- Suggests regex or parsing approach
- Recommends plot types
- Identifies metrics for general stats
```

## How Skills Work

### Automatic Invocation

Skills are **automatically discovered and used** by Claude based on:
1. **Context**: Keywords in conversation ("module request", "parse log file")
2. **Descriptions**: Each skill has a detailed description Claude matches against tasks
3. **File patterns**: Working with specific files or directories

You don't need to explicitly call skills - Claude invokes them when appropriate.

### Skill Structure

Each skill follows this pattern:
```
skill-name/
├── SKILL.md              # Main skill definition with YAML frontmatter
├── supporting-file.md    # Additional reference documentation
├── templates/            # Code templates
├── examples/             # Usage examples
└── scripts/              # Helper scripts (if needed)
```

### SKILL.md Format

```markdown
---
name: Skill Name in Gerund Form
description: What the skill does and when to use it
allowed-tools: [Read, Write, Edit, Bash]  # Optional
---

# Instructions

[Detailed instructions for Claude]

## Examples
[Usage examples]

## Related Resources
[Links to supporting files]
```

## Benefits Over GitHub Actions

### ✅ Local Execution
- Works in CLI and CI/CD
- No GitHub API rate limits
- Faster iteration

### ✅ More Flexible
- Automatic invocation based on context
- No manual workflow triggers
- Adapts to conversation flow

### ✅ Better Maintained
- Single source of truth
- Version controlled with code
- Easy to update and test

### ✅ Composable
- Skills work together
- Progressive disclosure of information
- Modular capabilities

## Migrating from GitHub Actions

### Old Workflow (GitHub Actions)
```yaml
# .github/workflows/module-requests.yml
on:
  issue_comment:
    types: [created]
# ... complex YAML configuration
```

### New Approach (Skills)
```
You: "@claude analyze this module request"
Claude: [Uses triaging skill automatically]
```

The skills replace these workflows:
- ❌ `.github/workflows/module-requests.yml` → ✅ `triaging-module-requests/`
- ❌ `.github/workflows/module-generation.yml` → ✅ `creating-multiqc-modules/`
- ❌ `.claude/commands/new-module.md` → ✅ `creating-multiqc-modules/` (automatic)

## Usage Examples

### Example 1: Triaging a Module Request

```
You: "Can you evaluate issue #3456 for priority?"

Claude: [Triaging skill activates]
1. Fetches issue details
2. Analyzes tool popularity (GitHub stars, downloads)
3. Evaluates community engagement
4. Assesses request quality
5. Generates score: 78/100 (High Priority)
6. Posts analysis comment on issue

Result: Comprehensive priority analysis posted
```

### Example 2: Creating a Module

```
You: "Implement the module for issue #7890"

Claude: [Module creation skill activates]
1. Downloads example files from issue
2. [Log analysis skill activates automatically]
   - Identifies TSV format
   - Finds sample name in filename
   - Detects key metrics
3. Generates module code
4. Creates tests
5. Updates configuration
6. Commits and creates PR

Result: Ready-to-review pull request
```

### Example 3: Understanding a Log File

```
You: "Help me understand this weird log format"

Claude: [Log analysis skill activates]
1. Examines file structure
2. Identifies multi-section format
3. Suggests parsing strategy
4. Recommends visualization types
5. Proposes general stats metrics

Result: Clear parsing plan and implementation guidance
```

## Testing Skills

### Verify Skills are Discoverable

```bash
# List all skills
ls -la .claude/skills/

# Check skill YAML frontmatter
head -10 .claude/skills/*/SKILL.md
```

### Test Automatic Invocation

Try these prompts:

1. **Triaging**: "Analyze module request issue #1234"
2. **Module Creation**: "Create a module for samtools"
3. **Log Analysis**: "Help me understand this log file format"

Claude should automatically use the appropriate skill.

## Customization

### Adding Project-Specific Context

Edit `CLAUDE.md` in project root to add:
- Project-specific patterns
- Common module types
- Team preferences
- Style guidelines

### Extending Skills

To add new capabilities:

1. Create new skill directory:
   ```bash
   mkdir -p .claude/skills/new-skill-name
   ```

2. Create SKILL.md with frontmatter:
   ```markdown
   ---
   name: New Skill Name
   description: What it does and when to use it
   ---

   # Instructions
   [Your instructions]
   ```

3. Add supporting files as needed

4. Test with Claude

## Troubleshooting

### Skill Not Being Invoked

**Problem**: Claude doesn't use a skill when expected

**Solutions**:
1. Make description more specific in YAML frontmatter
2. Use keywords from description in your prompt
3. Explicitly mention the task type
4. Check that SKILL.md has valid YAML frontmatter

### Skill Conflicts

**Problem**: Wrong skill activates

**Solutions**:
1. Make descriptions mutually exclusive
2. Add domain-specific keywords
3. Be more specific in prompt

### Skill Too Complex

**Problem**: Claude struggles with skill instructions

**Solutions**:
1. Break into multiple simpler skills
2. Move details to supporting files
3. Add more concrete examples
4. Simplify workflow steps

## Documentation

### Reference Documentation

- [Skills Overview](../docs/skills.md)
- [Best Practices](../docs/best-practices.md)
- [MultiQC Module Guidelines](../../CLAUDE.md)

### Skill-Specific Docs

- [Triaging Scoring Rubric](triaging-module-requests/scoring-rubric.md)
- [Module Template](creating-multiqc-modules/templates/module_template.py)
- [Log Format Patterns](analyzing-log-files/common-patterns.md)

## Contributing

When improving skills:

1. ✅ Keep main SKILL.md under 500 lines
2. ✅ Use progressive disclosure (link to supporting files)
3. ✅ Add concrete examples
4. ✅ Test with multiple scenarios
5. ✅ Document changes
6. ✅ Update version/date in SKILL.md

## Future Enhancements

Potential additions:
- **Module Testing** skill - Automated test generation and validation
- **Documentation Generation** skill - Auto-generate module documentation
- **Performance Analysis** skill - Profile and optimize module code
- **Integration Testing** skill - Test module against real data

## Support

If you encounter issues:
1. Check skill YAML frontmatter is valid
2. Verify file paths in supporting files
3. Test with explicit skill mention
4. Review skill description clarity
5. Check Claude Code logs

For questions or improvements, contact the MultiQC team or open an issue.

---

**Last Updated**: 2025-10-16
**Version**: 1.0.0
**Maintained by**: MultiQC Development Team
