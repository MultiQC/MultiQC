# MultiQC Contributing Guidelines

Hi there! Many thanks for taking an interest in improving MultiQC.

I try to manage the required tasks for MultiQC using GitHub issues, you probably came to this page when creating one. Most issues come in two flavours - either reporting a problem or requesting a feature. Please use the template prefilled into new issues when this is the case as it saves time. The most common reason for long-running issues is module requests without any example log files for example.

## Module Request Triage System

We've implemented an automated triage system for module requests to help manage the growing number of requests and ensure the most valuable modules are prioritized:

### How Module Requests Are Prioritized

Module requests are automatically scored and prioritized based on several criteria:

1. **Tool Popularity** (major factor): GitHub stars, forks, recent development activity
2. **Community Engagement**: 👍 reactions on the issue, meaningful comments and discussion
3. **Request Quality**: Complete information, quality example files, clear use case
4. **Technical Feasibility**: Well-structured tool output, parseable format
5. **Maintenance Status**: Active development, recent releases, documentation quality

### Automatic Actions

- **Weekly Triage**: Every Monday, an automated system analyzes all open module requests
- **Priority Labels**: Issues are tagged with `priority: high/medium/low` based on scoring
- **Stale Request Closure**: Low-priority requests with no activity for 6+ months may be automatically closed
- **Feedback & Guidance**: All analyzed requests receive constructive feedback for improvement

### Getting Your Request Prioritized

To increase your module request's priority score:

- **Provide complete information**: Fill all required fields thoroughly  
- **Upload quality example files**: Real tool outputs (not copy-pasted text)
- **Include tool repository**: Help us assess popularity and maintenance status
- **Explain the value**: Why would this module benefit the MultiQC community?
- **Engage constructively**: Participate in discussions, provide additional context
- **Use `@claude analyze-module`**: Get instant automated feedback on your request

### Interactive Analysis

You can get immediate analysis of any module request by commenting `@claude analyze-module` on the issue. This provides:

- Technical feasibility assessment
- Tool popularity metrics
- Priority scoring explanation  
- Specific improvement suggestions
- Comparison with existing modules

Closed requests can always be reopened if there's renewed interest or significant improvements to the tool/request.

However, don't be put off by this template - other more general issues and suggestions are welcome! Contributions to the code are even more welcome ;)

> _If you need help using MultiQC then the best place to go is the Seqera community forum where you can questions in the MultiQC category: [https://community.seqera.io](https://community.seqera.io/c/multiqc/6)_

## Contribution workflow

If you'd like to write some code for MultiQC, the standard workflow
is as follows:

1. Check that there isn't already an issue about your idea in the
   [MultiQC issues](https://github.com/MultiQC/MultiQC/issues) to avoid
   duplicating work.
   - Feel free to add a new issue here for the same reason.
2. Fork the MultiQC repository to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request and wait for the code to be reviewed and merged.

If you're not used to this workflow with git, you can start with some [basic docs from GitHub](https://help.github.com/articles/fork-a-repo/).

When it comes to MultiQC, please consult the [MultiQC documentation](https://docs.seqera.io/multiqc) and don't hesitate to get in touch on the [community forums](https://community.seqera.io) for help and feedback.

A few pointers to bear in mind:

- New modules should be _fast_
  - MultiQC modules parse log files, they _don't_ calculating new metrics (typically).
- New modules must scale well
  - Try to imagine what will happen if someone runs your module with 5000 samples

### Review workflow

Once you've submitted a new pull request, here's what you can expect from me:

- I usually don't look at your code at all until the automated tests pass
  - The tests use example data in the [test-data](https://github.com/MultiQC/test-data) repository, so you'll need some files there before the PR will go any further.
  - The tests use [GitHub Actions](https://github.com/features/actions) so should also run automatically on your fork.
- First pass - I go through and give feedback just by reading the code
- Second pass - I download and run your code, usually more feedback
- Merge! Once we're both happy, I merge into the main codebase.
