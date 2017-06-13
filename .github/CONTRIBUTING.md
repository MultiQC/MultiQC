# MultiQC Contributing Guidelines

Hi there! Many thanks for taking an interest in improving MultiQC.

I try to manage the required tasks for MultiQC using GitHub issues, you probably came to this page when creating one. Most issues come in two flavours - either reporting a problem or requesting a bug. Please use the template prefilled into new issues when this is the case as it saves time. The most common reason for long-running issues is module requests without any example log files for example.

However, don't be put off by this template - other more general issues and suggestions are welcome! Contributions to the code are even more welcome ;)

> _If you need help using MultiQC then the best place to go is the Gitter chatroom where you can ask me questions directly: https://gitter.im/ewels/MultiQC_

## Contribution workflow
If you'd like to write some code for MultiQC, the standard workflow
is as follows:

1. Check that there isn't already an issue about your idea in the
   [MultiQC issues](https://github.com/ewels/MultiQC/issues) to avoid
   duplicating work.
    * Feel free to add a new issue here for the same reason.
2. Fork the MultiQC repository to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request and wait for the code to be reviewed and merged.

If you're not used to this workflow with git, you can start with some [basic docs from GitHub](https://help.github.com/articles/fork-a-repo/) or even their [excellent interactive tutorial](https://try.github.io/).

When it comes to MultiQC, please consult the [MultiQC documentation](http://multiqc.info/docs/) and don't hesitate to get in touch on [gitter](https://gitter.im/ewels/MultiQC) for help and feedback.

A few pointers to bear in mind:

* New modules should be _fast_
    * MultiQC modules parse log files, they _don't_ calculating new metrics (typically).
* New modules must scale well
    * Try to imagine what will happen if someone runs your module with 5000 samples
* Code must run on both Python 2 and 3

### Review workflow
Once you've submitted a new pull request, here's what you can expect from me:

* I usually don't look at your code at all until the automated tests pass
   * The tests use example data in the [MultiQC_TestData](https://github.com/ewels/MultiQC_TestData) repository, so you'll need some files there before the PR will go any further.
   * You can set up [Travis](https://travis-ci.org) to run the same tests on your fork really easily - just enable the repo.
* First pass - I go through and give feedback just by reading the code
* Second pass - I download and run your code, usually more feedback
* Merge! Once we're both happy, I merge into the main codebase.
