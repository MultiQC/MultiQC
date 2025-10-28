# Claude Code Agent Skills

Agent Skills are modular capabilities that extend Claude's functionality. They are automatically invoked by Claude based on context, not user commands.

## Key Characteristics

- **Modular**: Each Skill focuses on a specific capability
- **Discoverable**: Claude automatically finds and uses Skills based on their descriptions
- **Flexible**: Can be personal (individual use), project-based (team use), or plugin-based
- **Composable**: Consist of instructions, optional scripts, and supporting resources

## Skill Structure

Skills are stored in organized folders with a `SKILL.md` file:

```
.claude/skills/{skill-name}/
├── SKILL.md              # Main skill definition (required)
├── supporting-file.md    # Additional reference files
├── scripts/              # Helper scripts
└── templates/            # Code templates
```

## SKILL.md Format

Every skill requires a `SKILL.md` file with YAML frontmatter:

```markdown
---
name: Processing PDFs
description: Extracts text and metadata from PDF files using pdfplumber. Use when the user needs to analyze or extract information from PDF documents.
---

# Instructions

[Your detailed instructions here]

## Examples

[Usage examples]

## Related Resources

- [reference-file.md](./reference-file.md)
```

### YAML Frontmatter Fields

- **name** (required): Skill name in gerund form (e.g., "Processing PDFs")
- **description** (required): When and how to use this skill (third person)
- **allowed-tools** (optional): Restrict which tools Claude can use

### Description Best Practices

Write clear, specific descriptions that help Claude understand:

- **What** the skill does
- **When** to use it
- **What domain** it applies to

Examples:

- ✅ "Analyzes Python codebases to identify performance bottlenecks using profiling tools"
- ❌ "Helps with code" (too vague)

## Skill Creation Process

1. **Create Directory**: `mkdir -p .claude/skills/your-skill-name`
2. **Write SKILL.md**: Include frontmatter and instructions
3. **Add Supporting Files**: Scripts, templates, reference docs
4. **Test**: Verify Claude discovers and uses the skill appropriately

## Skill Types

### Personal Skills

- Location: `~/.claude/skills/`
- Scope: Available to you across all projects
- Use for: Personal workflows, preferences, common patterns

### Project Skills

- Location: `.claude/skills/` (project root)
- Scope: Available to all project contributors
- Use for: Project-specific patterns, team workflows
- Version control: Check into git

### Plugin Skills

- Distribution: NPM packages
- Scope: Shareable across organizations
- Use for: Standardized capabilities, company-wide patterns

## Best Practices

### 1. Keep Skills Focused

Each skill should handle one specific capability. Break complex tasks into multiple skills.

### 2. Use Progressive Disclosure

- Keep main SKILL.md under 500 lines
- Move detailed reference material to separate files
- Reference additional files using relative links

### 3. Write Concise Instructions

- Only include information Claude doesn't already know
- Challenge each sentence: "Does Claude really need this?"
- Remove redundant explanations

### 4. Include Examples

Show Claude concrete examples of:

- Input/output formats
- Common usage patterns
- Edge cases to handle

### 5. Create Clear Workflows

For multi-step processes:

- Number the steps
- Use checklists for tracking
- Create feedback loops
- Include validation steps

### 6. Handle Errors Explicitly

- Provide error handling scripts
- Define fallback behaviors
- Create verifiable outputs

## Tool Restriction

Control which tools Claude can use in a skill:

```yaml
---
name: Safe File Analysis
description: Analyzes files without making modifications
allowed-tools: [Read, Grep, Glob]
---
```

Available tools: Read, Write, Edit, Bash, WebFetch, WebSearch, Task, etc.

## Advanced Features

### Supporting Scripts

Create utility scripts to handle complex operations:

```
.claude/skills/data-analysis/
├── SKILL.md
└── scripts/
    ├── process_csv.py
    └── validate_output.sh
```

Reference in SKILL.md:

````markdown
Run the processing script:

```bash
python .claude/skills/data-analysis/scripts/process_csv.py input.csv
```
````

```

### Code Templates

Provide templates for code generation:

```

.claude/skills/api-client/
├── SKILL.md
└── templates/
├── client.py.template
└── test.py.template

````

### Reference Documentation

Link to additional context:

```markdown
For detailed API specifications, see [api-reference.md](./api-reference.md)
````

## Troubleshooting

### Skill Not Being Used

**Problem**: Claude doesn't invoke your skill when expected

**Solutions**:

1. Make description more specific
2. Add domain keywords to description
3. Include example trigger scenarios
4. Test with explicit mentions

### Multiple Skills Conflicting

**Problem**: Two skills trigger for the same task

**Solutions**:

1. Make descriptions mutually exclusive
2. Define clear boundaries in descriptions
3. Merge related skills

### Skill Too Complex

**Problem**: Claude struggles to follow skill instructions

**Solutions**:

1. Break into multiple simpler skills
2. Use progressive disclosure for details
3. Add more examples
4. Simplify workflow steps

## Examples

### Simple Skill (Inline)

```markdown
---
name: Formatting Commit Messages
description: Formats git commit messages following conventional commits specification. Use when creating commits.
---

# Conventional Commits Format

Format: `<type>(<scope>): <subject>`

Types:

- feat: New feature
- fix: Bug fix
- docs: Documentation
- refactor: Code refactoring
- test: Tests
- chore: Maintenance

Example: `feat(auth): add OAuth2 support`
```

### Complex Skill (With Supporting Files)

```
.claude/skills/api-documentation/
├── SKILL.md
├── openapi-template.yaml
├── examples/
│   ├── rest-api.md
│   └── graphql-api.md
└── scripts/
    └── validate-openapi.py
```

## Sharing Skills

### Via Git

Commit project skills to repository:

```bash
git add .claude/skills/
git commit -m "Add project skills"
```

### Via NPM Plugin

Package skills for distribution:

```json
{
  "name": "@company/claude-skills",
  "version": "1.0.0",
  "claudePlugin": {
    "skills": ["skills/"]
  }
}
```

## Related Documentation

- [Best Practices Guide](./best-practices.md)
- [Tool Reference](https://docs.claude.com/en/docs/claude-code/tools)
- [Plugin Development](https://docs.claude.com/en/docs/claude-code/plugins)
