# Agent Skills Best Practices

This guide provides comprehensive best practices for authoring effective Claude Code Agent Skills.

## Core Principles

### 1. Conciseness is Key

**Only add context Claude doesn't already know**

❌ Bad:

```markdown
Python is a programming language. To create a virtual environment,
which is an isolated Python environment that helps manage dependencies...
```

✅ Good:

````markdown
Create virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
```
````

````

**Challenge every sentence**: "Does Claude really need this?"

### 2. Skill Structure Best Practices

#### Naming Convention
- Use gerund form: "Processing PDFs", "Analyzing Logs", "Creating Modules"
- Be specific, not generic: "Parsing FASTQ Files" not "File Processing"

#### Description Guidelines
- Write in third person
- Explain what AND when to use
- Include domain keywords
- Be concrete, not abstract

Examples:

✅ Good:
```yaml
description: Analyzes MultiQC module requests from GitHub issues, evaluates tool popularity via GitHub API, and generates priority scores. Use when triaging module request issues.
````

❌ Bad:

```yaml
description: Helps with issues
```

### 3. Progressive Disclosure

**Organize content across multiple files**

```
.claude/skills/complex-skill/
├── SKILL.md              # Main entry point (<500 lines)
├── reference/
│   ├── api-details.md    # Detailed API reference
│   └── examples.md       # Extensive examples
└── scripts/
    └── helper.py         # Utility scripts
```

**In SKILL.md, reference additional files**:

```markdown
For detailed API specifications, see [api-details.md](./reference/api-details.md)
```

**Keep main SKILL.md under 500 lines**

### 4. Workflow Design

#### Clear Sequential Steps

✅ Good:

```markdown
## Implementation Process

1. **Analyze Issue**

   - Extract tool name and homepage
   - Download example files
   - Identify data patterns

2. **Generate Module**

   - Create directory structure
   - Write parsing code
   - Generate visualizations

3. **Create Tests**

   - Write unit tests
   - Verify output format
   - Run test suite

4. **Submit PR**
   - Commit changes
   - Push to branch
   - Create pull request
```

#### Use Checklists for Tracking

```markdown
## Pre-submission Checklist

- [ ] Module follows naming convention
- [ ] Tests pass with 100% coverage
- [ ] Documentation is complete
- [ ] Entry points registered
- [ ] Example files included
```

#### Implement Feedback Loops

```markdown
After generating code:

1. Run linters and fix issues
2. Run tests
3. If tests fail, analyze failures and regenerate
4. Repeat until tests pass
```

### 5. Code and Script Best Practices

#### Handle Errors Explicitly

❌ Bad:

```python
def parse_file(path):
    with open(path) as f:
        return json.load(f)
```

✅ Good:

```python
def parse_file(path):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        log.error(f"File not found: {path}")
        return None
    except json.JSONDecodeError as e:
        log.error(f"Invalid JSON in {path}: {e}")
        return None
```

#### Provide Utility Scripts

````markdown
Use this validation script to verify output:

```bash
python .claude/skills/module-creation/scripts/validate_module.py multiqc/modules/newtool/
```
````

````

#### Create Verifiable Intermediate Outputs

```markdown
After parsing, write intermediate JSON:
```python
with open('parsed_data.json', 'w') as f:
    json.dump(data, f, indent=2)
````

Verify the structure matches expected schema before proceeding.

````

#### Package Dependencies

```markdown
Required dependencies:
```bash
pip install requests beautifulsoup4 pydantic
````

````

### 6. Examples and Documentation

#### Include Concrete Examples

```markdown
## Example: Parsing FastQC Output

Input file (`fastqc_data.txt`):
````

##FastQC 0.11.9

> > Basic Statistics pass
> > Filename sample1.fastq
> > Total Sequences 1000000

````

Expected output:
```python
{
    "sample_name": "sample1",
    "total_sequences": 1000000,
    "version": "0.11.9"
}
````

````

#### Show Edge Cases

```markdown
## Edge Cases to Handle

1. **Missing version**: Use "unknown"
2. **Duplicate samples**: Append suffix `_1`, `_2`
3. **Empty files**: Skip with warning
4. **Malformed data**: Log error and continue
````

## Anti-Patterns to Avoid

### 1. ❌ Offering Too Many Options

Bad:

```markdown
You can either:

- Option A: Do this...
- Option B: Or do that...
- Option C: Or maybe this...
```

Good:

```markdown
Do this:
[Single clear instruction]
```

### 2. ❌ Using Windows-Style Paths

Bad:

```markdown
Read `C:\Users\...\file.txt`
```

Good:

```markdown
Read `/path/to/file.txt`
Use forward slashes on all platforms
```

### 3. ❌ Including Time-Sensitive Information

Bad:

```markdown
As of January 2024, use Python 3.12...
```

Good:

```markdown
Use Python 3.9+
```

### 4. ❌ Punting Error Handling to Claude

Bad:

```markdown
If something goes wrong, figure it out
```

Good:

```markdown
If parsing fails:

1. Check file format matches expected pattern
2. Try alternate parsing strategy (see fallback.md)
3. If still failing, skip file and log warning
```

### 5. ❌ Over-Explaining Basic Concepts

Bad:

```markdown
GitHub is a platform for hosting code repositories. A repository,
or "repo" for short, is a collection of files...
```

Good:

````markdown
Clone the repository:

```bash
git clone <url>
```
````

````

## Workflow for Creating Skills

### Recommended Process

1. **Complete Task Without Skill**
   - Work through task manually
   - Document what works
   - Note repeated patterns

2. **Identify Reusable Patterns**
   - What steps are always the same?
   - What context is consistently needed?
   - What errors commonly occur?

3. **Draft Skill**
   - Write initial SKILL.md
   - Include examples from your work
   - Add error handling for issues you hit

4. **Review for Conciseness**
   - Remove redundant explanations
   - Challenge every paragraph
   - Move details to supporting files

5. **Improve Information Architecture**
   - Organize into logical sections
   - Create supporting files for depth
   - Link related resources

6. **Test on Similar Tasks**
   - Try skill on different examples
   - Verify Claude discovers it appropriately
   - Check if skill is sufficient

7. **Iterate Based on Observation**
   - Watch how Claude uses skill
   - Add missing context
   - Remove unused information

## Testing Skills

### Test Discoverability

Ask Claude to perform related tasks WITHOUT mentioning the skill:
- Does Claude automatically use it?
- Does it trigger at appropriate times?
- Are there false positives?

### Test Multiple Models

If available, test with different Claude models to ensure compatibility.

### Create Evaluation Scenarios

```markdown
## Test Cases

### Scenario 1: Simple Module Request
- Issue: #123 (samtools stats)
- Expected: Generate complete module
- Verify: Tests pass, follows standards

### Scenario 2: Complex Log Format
- Issue: #456 (multi-section output)
- Expected: Parse all sections
- Verify: All metrics captured

### Scenario 3: Missing Information
- Issue: #789 (no example files)
- Expected: Request more info
- Verify: Helpful guidance provided
````

### Observe Navigation Patterns

Watch how Claude:

- Discovers the skill
- Navigates supporting files
- Uses provided scripts
- Handles errors

## Skill Maintenance

### Version History

Track changes in SKILL.md:

```markdown
## Changelog

### 2024-01-15

- Added support for JSON log formats
- Improved error messages
- Updated examples

### 2024-01-01

- Initial version
```

### Keep Dependencies Updated

```markdown
## Dependencies

Last verified: 2024-01-15

- Python 3.9+
- requests>=2.28.0
- pydantic>=2.0.0
```

### Regular Review

Schedule periodic reviews:

- Are instructions still current?
- Have patterns changed?
- Is skill still being used?
- Can it be simplified?

## Examples of Well-Structured Skills

### Example 1: Simple, Focused Skill

````markdown
---
name: Formatting Docstrings
description: Formats Python docstrings using Google style guide. Use when writing or updating Python documentation.
---

# Google Style Docstrings

Format:

```python
def function(arg1: str, arg2: int) -> bool:
    """Short description.

    Longer description if needed.

    Args:
        arg1: Description of arg1
        arg2: Description of arg2

    Returns:
        Description of return value

    Raises:
        ValueError: When validation fails
    """
```
````

Follow PEP 257 conventions.

```

### Example 2: Complex Skill with Structure

```

.claude/skills/api-testing/
├── SKILL.md # Main workflow (<300 lines)
├── reference/
│ ├── http-methods.md # HTTP verb reference
│ ├── status-codes.md # Status code meanings
│ └── auth-patterns.md # Authentication patterns
├── templates/
│ ├── test-case.py # Test template
│ └── mock-response.json # Mock data template
└── scripts/
├── run-tests.sh # Test runner
└── validate-spec.py # OpenAPI validator

```

## Related Resources

- [Skills Documentation](./skills.md)
- [Claude Code Tools Reference](https://docs.claude.com/en/docs/claude-code/tools)
- [Example Skills Repository](https://github.com/anthropics/claude-code-examples)
```
