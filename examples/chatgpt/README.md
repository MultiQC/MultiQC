# MultiQC reports in ChatGPT

This example provides one MCP tool, `run_multiqc`. It accepts files attached to a
ChatGPT conversation, runs MultiQC, and returns links to the interactive HTML
report and full parsed JSON, plus a statistics preview for ChatGPT to summarize.
The report is downloaded and opened in a browser.

## Run locally

From the MultiQC repository, with [uv](https://docs.astral.sh/uv/) installed:

```bash
uv run --project examples/chatgpt multiqc-chatgpt --public-url https://YOUR-TUNNEL-HOST
```

The example uses Python 3.10 or newer and the MultiQC checkout. It binds to
`127.0.0.1:8000`; use `--port` to change the port. Connect an HTTPS tunnel or
authenticated ingress to that local port, and use its origin for `--public-url`.
Starting the command does not create a tunnel or connect a ChatGPT account.

In ChatGPT developer mode, add the server at `https://YOUR-TUNNEL-HOST/mcp`.
See the [official connection instructions](https://developers.openai.com/plugins/build/app-quickstart#connect-your-mcp-server-in-chatgpt).
Attach your QC logs, select the connection, and ask:

> Create a MultiQC report from these attachments. Give me the report and JSON
> download links, then summarize the observed metrics and any missing evidence.

Original file names help MultiQC recognize outputs. Individual tool logs,
FastQC ZIP files, and MultiQC custom content files work as usual. A generic ZIP
of an analysis directory is not unpacked. If no supported analysis is found, the
tool returns an error.

## Output and limits

The tool returns `report_url`, `data_url`, `sample_count`, `modules`, and
`general_stats`. The preview contains up to 20 samples and 40 metrics per sample.
`general_stats_sample_count` counts samples with general statistics;
`general_stats_truncated` identifies an incomplete preview. General statistics
use raw values; the HTML report can display scaled units. Use the JSON download
for additional analysis. The preview does not establish sample QC pass/fail or
clinical conclusions.

Each request accepts 1 to 50 files, up to 100 MiB in total. File downloads must
use HTTPS URLs under `oaiusercontent.com`, supplied through ChatGPT's
[file input contract](https://developers.openai.com/plugins/reference#define-file-inputs).
Redirects are rejected and download credentials are not forwarded. Downloads
have a two minute deadline; MultiQC runs have a five minute deadline. Two runs
can execute concurrently in separate processes. Uploaded files are removed after
each run. Generated artifacts remain in a temporary directory until shutdown.

This is a single-user development example with no authentication. Report URLs
contain unguessable identifiers and grant access to anyone who has the link.
Use your own data in a private development connection. A shared deployment needs
authentication, per-user artifact access, and artifact expiry at its ingress or
service boundary before it is exposed. MultiQC's AI summaries and MegaQC uploads
are disabled for these runs.

## Check the example

```bash
uv run --project examples/chatgpt --extra dev pytest examples/chatgpt/test_server.py
```

The check uses synthetic QC data and a simulated attachment download, then runs
the real MultiQC parser and renderer through an MCP client. It also downloads the
generated artifacts through the HTTP routes. A live ChatGPT connection requires
the tunnel and account setup above.
