# Claude Skills Quick Start Guide

This guide helps you get started with the new MultiQC Claude skills system.

## What Changed?

### Before (GitHub Actions)

- Module triage workflow: `.github/workflows/module-requests.yml`
- Module generation workflow: `.github/workflows/module-generation.yml`
- Slash command: `.claude/commands/new-module.md`
- Manual triggers and workflow dispatch

### After (Claude Skills)

- Automatic skill invocation based on context
- Works locally in Claude Code CLI
- Simpler, more flexible workflow
- All context in `.claude/skills/`

## Quick Examples

### 1. Triage a Module Request

**Old way**:

```
Comment on issue: "@claude analyze-module"
Wait for GitHub Actions to run
```

**New way**:

```
You: "Can you analyze issue #1234 for priority?"
Claude: [Automatically uses triaging skill and responds immediately]
```

### 2. Create a Module

**Old way**:

```
Manually trigger workflow
Provide issue number
Wait for Actions to complete
```

**New way**:

```
You: "Create a module for the tool in issue #5678"
Claude: [Automatically downloads files, analyzes, and creates module]
```

### 3. Understand a Log File

**New capability** (wasn't possible before):

```
You: "Help me understand how to parse this log file"
Claude: [Analyzes format and suggests parsing strategy]
```

## Available Skills

### 🎯 Triaging Module Requests

**Trigger words**: "triage", "prioritize", "analyze module request", "evaluate issue"

**What it does**:

- Fetches tool popularity metrics
- Evaluates community engagement
- Calculates priority score (0-100)
- Generates analysis report

**Example**:

```
"What's the priority of issue #3456?"
```

### 🔨 Creating MultiQC Modules

**Trigger words**: "create module", "implement module", "new module", "generate module"

**What it does**:

- Downloads example files
- Analyzes log format
- Generates complete module code
- Creates tests
- Submits PR

**Example**:

```
"Implement the samtools stats module from issue #7890"
```

### 🔍 Analyzing Log Files

**Trigger words**: "parse", "understand format", "log file", "analyze output"

**What it does**:

- Identifies file format
- Suggests parsing strategy
- Recommends visualizations
- Plans general stats

**Example**:

```
"How should I parse this tool's output?"
```

## Testing Your Setup

### Verify Skills are Installed

```bash
# List all skills
ls .claude/skills/

# Should show:
# analyzing-log-files/
# creating-multiqc-modules/
# triaging-module-requests/
# README.md
```

### Test Skill Discovery

Ask Claude:

```
"What skills do you have available for MultiQC?"
```

Claude should mention the three skills.

### Test Auto-Invocation

Try each skill with a natural prompt:

1. **Triaging**: "Analyze the module request in issue #1234"
2. **Module Creation**: "Create a module for fastqc"
3. **Log Analysis**: "How do I parse a TSV file?"

## Common Workflows

### Workflow 1: New Module Request Comes In

**Steps**:

1. New issue created with `module: new` label
2. Ask Claude: "Can you triage issue #XXXX?"
3. Claude automatically:
   - Fetches tool popularity
   - Evaluates request quality
   - Generates priority score
   - Posts analysis

### Workflow 2: Implement High-Priority Module

**Steps**:

1. High-priority module identified
2. Ask Claude: "Create the module for issue #XXXX"
3. Claude automatically:
   - Downloads example files
   - Analyzes format
   - Generates code
   - Creates tests
   - Submits PR

### Workflow 3: Debug Parsing Issue

**Steps**:

1. Module not parsing correctly
2. Ask Claude: "Help me fix the parsing in this module"
3. Claude automatically:
   - Examines log file
   - Identifies format issues
   - Suggests corrections
   - Updates code

## Tips for Best Results

### ✅ Do

1. **Be specific about what you want**

   - Good: "Analyze issue #1234 for module priority"
   - Bad: "Look at this issue"

2. **Provide context**

   - Good: "Create a module for samtools stats from issue #5678"
   - Bad: "Make a module"

3. **Reference files directly**

   - Good: "Parse the example file at path/to/file.txt"
   - Bad: "Parse that file"

4. **Use natural language**
   - Skills respond to conversational prompts
   - No need for special syntax

### ❌ Don't

1. **Don't try to manually invoke skills**

   - Let Claude choose automatically
   - Trust the automatic discovery

2. **Don't combine too many tasks**

   - Break complex work into steps
   - Let Claude use multiple skills

3. **Don't override skill decisions**
   - Skills follow best practices
   - Trust the automated workflow

## Troubleshooting

### Skill Not Being Used

**Symptom**: Claude doesn't seem to use a skill

**Solutions**:

1. Use keywords from skill descriptions
2. Be more specific about the task
3. Mention "module request" or "log file" explicitly

**Example**:

```
❌ "What about this?"
✅ "Analyze this module request for priority"
```

### Wrong Skill Activates

**Symptom**: Claude uses unexpected skill

**Solutions**:

1. Be more specific about your intent
2. Mention the specific task type
3. Provide more context

**Example**:

```
❌ "Check this file"
✅ "Parse this samtools log file and extract metrics"
```

### Need to See What's Happening

**Check**:

```bash
# View skill instructions
cat .claude/skills/triaging-module-requests/SKILL.md

# Check supporting files
cat .claude/skills/triaging-module-requests/scoring-rubric.md
```

## Migration from GitHub Actions

### Deprecated Workflows

These files are **no longer needed**:

- `.github/workflows/module-requests.yml` (in feature branches)
- `.github/workflows/module-generation.yml` (in feature branches)
- `.claude/commands/new-module.md` (in feature branches)

### What to Keep

These files remain useful:

- `.github/workflows/claude.yml` - General @claude integration
- `CLAUDE.md` - Project documentation for Claude
- `.claude/settings.local.json` - Your personal settings

### Clean Up

If you have old branches with Actions:

```bash
# These workflows are replaced by skills
# You can clean them up or leave them as reference
```

## Advanced Usage

### Combining Skills

Skills work together automatically:

```
You: "Create a module for the tool in issue #7890"

Claude:
1. [Triaging skill] - Checks if high priority
2. [Log analysis skill] - Examines example files
3. [Module creation skill] - Generates implementation
```

### Customizing Skills

Edit SKILL.md files to:

- Add project-specific patterns
- Include team preferences
- Add new examples
- Refine descriptions

### Creating New Skills

Follow the template:

```bash
mkdir .claude/skills/new-skill-name
cat > .claude/skills/new-skill-name/SKILL.md << 'EOF'
---
name: New Skill Name
description: What it does and when to use it
---

# Instructions
[Your instructions here]
EOF
```

## Getting Help

### Documentation

- Full guide: `.claude/skills/README.md`
- Skills overview: `.claude/docs/skills.md`
- Best practices: `.claude/docs/best-practices.md`

### Examples

- Simple module: `.claude/skills/creating-multiqc-modules/examples/simple-module.md`
- Complex module: `.claude/skills/creating-multiqc-modules/examples/complex-module.md`
- Log patterns: `.claude/skills/analyzing-log-files/common-patterns.md`

### Support

- Open an issue on GitHub
- Ask Claude directly: "How do I use the triaging skill?"
- Review skill documentation: `cat .claude/skills/*/SKILL.md`

## What's Next?

Now that skills are set up:

1. **Try them out**: Test with real module requests
2. **Provide feedback**: What works? What needs improvement?
3. **Extend them**: Add project-specific patterns
4. **Share knowledge**: Help team members understand the new system

## Success Metrics

You'll know skills are working when:

✅ Module triage happens in seconds, not minutes
✅ Module creation is consistent and complete
✅ Log file parsing is systematic and reliable
✅ Less time spent on routine tasks
✅ More time for feature development

---

**Questions?** Ask Claude: "How do MultiQC skills work?"
