# Module Triage Project Setup Guide

This guide explains how to set up a GitHub project to track and manage module requests as part of the automated triage system.

## Project Structure

### Recommended Project Setup

1. **Create a new GitHub Project** in the MultiQC organization
   - Name: "Module Request Triage"
   - Description: "Automated tracking and prioritization of MultiQC module requests"
   - Template: "Board" view

2. **Configure Project Columns**

   | Column Name               | Purpose                                    | Automation                                    |
   | ------------------------- | ------------------------------------------ | --------------------------------------------- |
   | **Needs Analysis**        | New requests awaiting initial triage       | Auto-add new issues with "module: new" label  |
   | **Low Priority**          | Score < 40                                 | Auto-move when `priority: low` label added    |
   | **Medium Priority**       | Score 40-69                                | Auto-move when `priority: medium` label added |
   | **High Priority**         | Score ≥ 70                                 | Auto-move when `priority: high` label added   |
   | **Needs Examples**        | Missing example files                      | Auto-move when `needs-examples` label added   |
   | **Ready for Development** | Complete requests ready for implementation | Manual move after maintainer review           |
   | **In Development**        | Active module development                  | Manual move when work begins                  |
   | **Complete**              | Finished modules                           | Auto-move when issue closed                   |

3. **Project Automation Rules**

   **Auto-add items:**

   ```
   When: Issue is opened
   Filters: Label contains "module: new"
   Then: Add to project in "Needs Analysis"
   ```

   **Auto-move by priority:**

   ```
   When: Label is added to item
   Filters: Label is "priority: low"
   Then: Move to "Low Priority"
   ```

   ```
   When: Label is added to item
   Filters: Label is "priority: medium"
   Then: Move to "Medium Priority"
   ```

   ```
   When: Label is added to item
   Filters: Label is "priority: high"
   Then: Move to "High Priority"
   ```

   **Auto-move for missing examples:**

   ```
   When: Label is added to item
   Filters: Label is "needs-examples"
   Then: Move to "Needs Examples"
   ```

   **Auto-archive completed:**

   ```
   When: Issue is closed
   Filters: Label contains "module: new"
   Then: Move to "Complete"
   ```

## Custom Fields (Optional)

Add these custom fields to track additional metadata:

- **Priority Score** (Number): Automated score from triage system
- **Tool Popularity** (Number): GitHub stars count
- **Has Examples** (Boolean): Whether example files are provided
- **Last Analysis** (Date): When last triage analysis was performed

## Labels for Triage System

Ensure these labels exist in the repository:

### Priority Labels

- `priority: high` (Red) - Score ≥ 70
- `priority: medium` (Yellow) - Score 40-69
- `priority: low` (Gray) - Score < 40

### Status Labels

- `needs-triage` (Purple) - Awaiting automated analysis
- `needs-examples` (Orange) - Missing example files
- `stale` (Gray) - No activity for 6+ months
- `ready-for-dev` (Green) - Complete and ready for implementation

### Quality Labels

- `good-first-module` (Light green) - Good for new contributors
- `complex-module` (Dark blue) - Requires advanced implementation
- `popular-tool` (Gold) - Tool has high GitHub stars

## Workflow Integration

The unified `module-requests.yml` workflow handles:

1. **Individual Analysis**: Immediate analysis of new requests and `@claude analyze-module` commands
2. **Weekly Bulk Triage**: Comprehensive review of all open requests every Monday
3. **Project Synchronization**: Automatic issue tracking and board management
4. **Flexible Execution**: Manual triggers with different modes (analyze-single, triage-all, dry-run)

## Manual Project Management

### Maintainer Actions

1. **Review High Priority** column weekly
2. **Move items to "Ready for Development"** after verification
3. **Assign items in "In Development"** when work begins
4. **Update custom fields** with additional context

### Community Engagement

- **Pin important issues** that need community input
- **Add milestone labels** for release planning
- **Use project discussions** for coordination

## Analytics and Reporting

The project provides visibility into:

- **Request volume trends** (new requests per month)
- **Priority distribution** (high/medium/low breakdown)
- **Processing time** (time from request to development)
- **Success rate** (requests that become modules)

## Integration with Automated Triage

### Workflow Permissions

The unified workflow includes all necessary permissions:

```yaml
permissions:
  issues: write # Add labels, post comments
  projects: write # Project board management
  contents: read # Repository access
```

### Environment Variables (if using organization-level project)

```yaml
env:
  PROJECT_ID: "PVT_kwDOAxxxxxxx" # Get from project URL
  ORGANIZATION: "MultiQC"
```

## Troubleshooting

### Common Issues

1. **Items not auto-adding to project**
   - Check automation rules are enabled
   - Verify label filters match exactly
   - Ensure workflow has `projects: write` permission

2. **Cards not moving between columns**
   - Verify automation rules for label-based moves
   - Check that labels are being added correctly by workflows
   - Manual moves may be needed initially

3. **Missing priority labels**
   - Run the unified workflow manually: `workflow_dispatch` with `triage-all` mode
   - Check that `module-requests.yml` is running on new issues

### Manual Sync

If automation fails, manually add issues to the project:

```bash
# Use GitHub CLI to bulk-add issues
gh project item-add PROJECT_NUMBER --url ISSUE_URL
```

## Benefits

This project structure provides:

- **Visual workflow** for module request lifecycle
- **Automated prioritization** based on objective criteria
- **Clear backlogs** for maintainers and contributors
- **Progress tracking** from request to implementation
- **Community transparency** into development priorities

The project becomes the single source of truth for module request status and facilitates better coordination between maintainers and contributors.
