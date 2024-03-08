name: Bug report
description: File a bug report if MultiQC is breaking / not behaving in the way you expect
body:
  - type: textarea
    id: description
    attributes:
      label: Description of bug
      description: |
        A clear and concise description of what the bug is.
        **Please do note paste in the contents of example files** - upload them as a file.
    validations:
      required: true

  - type: textarea
    id: error-file
    attributes:
      label: File that triggers the error
      description: |
        Please drag and drop (and upload to the GitHub issue) an input file that can be used to replicate the error.
        ***Please do not copy and paste log contents***, as important whitespace can change.
        If the file type is not allowed, please compress into a `.zip` file.
      placeholder: |
        [ Drag and drop an example file here to upload ]
        Please do not copy and paste log contents! Use a .zip file if needed.

  - type: textarea
    id: log
    attributes:
      label: MultiQC Error log
      description: |
        Please paste your **full command** and MultiQC log.
        Please do not truncate it, as the MultiQC version and command may be needed to help.
      render: console
      placeholder: |
        $ multiqc .

          /// MultiQC 🔍 | v1.10.1

        |           multiqc | Search path : /path/to/my/data
        |         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 2352/2352
        |            fastqc | Found 498 reports
        |           multiqc | Compressing plot data
        |           multiqc | Report      : multiqc_report.html
        |           multiqc | Data        : multiqc_data
        |           multiqc | MultiQC complete

  - type: checkboxes
    id: checklist
    attributes:
      label: Before submitting
      description: >-
        Please ensure your bug report fulfills all of the following requirements.
      options:
        - label: >-
            I have read the [troubleshooting documentation](https://multiqc.info/docs/usage/troubleshooting/).
          required: true
        - label: >-
            I am using the latest release of MultiQC.
          required: true
        - label: >-
            I have included a full MultiQC log, not truncated.
          required: true
        - label: >-
            I have attached an input file (**.zip** if necessary) that triggers the error.
          required: true
