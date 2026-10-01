"""A standalone MCP example for downloadable MultiQC reports in ChatGPT."""

import asyncio
import json
import logging
import shutil
import sys
from pathlib import Path, PureWindowsPath
from tempfile import TemporaryDirectory
from typing import Annotated, Any, Dict, List
from urllib.parse import urlsplit
from uuid import uuid4

import click
import httpx
from mcp.server.fastmcp import FastMCP
from mcp.server.transport_security import TransportSecuritySettings
from mcp.types import CallToolResult, ResourceLink, TextContent, ToolAnnotations
from pydantic import BaseModel, ConfigDict, Field
from starlette.requests import Request
from starlette.responses import FileResponse, Response

MAX_INPUT_BYTES = 100 * 1024 * 1024
RUN_TIMEOUT = 300
PREVIEW_SAMPLES = 20
PREVIEW_METRICS = 40

# Use MultiQC's public API in a fresh process because its report state is global.
WORKER = """
import json
import sys
from html import escape
from pathlib import Path
import multiqc
from multiqc.plots.table_object import Cell
from multiqc.utils.util_functions import dump_json

inputs, output, title = sys.argv[1:]
result = multiqc.run(inputs, cfg=multiqc.ClConfig(
    output_dir=output, filename="multiqc_report.html", title=escape(title),
    template="default", make_report=True, make_data_dir=True, data_format="json",
    plots_force_interactive=True, strict=True, no_ai=True, no_megaqc_upload=True, no_version_check=True,
    cl_config=[json.dumps({"output_fn_name": "multiqc_report.html", "data_dir_name": "multiqc_data",
                          "data_dump_file": True, "data_dump_file_write_raw": True,
                          "show_analysis_paths": False, "make_pdf": False, "development": False,
                          "zip_data_dir": False, "export_plots": False})],
))
if result.sys_exit_code or result.message:
    sys.exit(result.message or "MultiQC failed; see the server log")
stats = {sample: {key: value.raw if isinstance(value, Cell) else value
                  for key, value in row.items()}
         for sample, row in multiqc.get_general_stats_data().items()}
with (Path(output) / "summary.json").open("w", encoding="utf-8") as handle:
    dump_json({"sample_count": len(multiqc.list_samples()),
               "modules": multiqc.list_modules(), "general_stats": stats}, handle)
"""


class ChatGPTFile(BaseModel):
    """The file schema required by ChatGPT's openai/fileParams metadata."""

    model_config = ConfigDict(extra="forbid")
    download_url: str
    file_id: str
    mime_type: str = ""
    file_name: str = ""


class ReportResult(BaseModel):
    report_url: str
    data_url: str
    sample_count: int
    modules: List[str]
    general_stats: Dict[str, Dict[str, Any]]
    general_stats_sample_count: int
    general_stats_truncated: bool


async def download_files(files: List[ChatGPTFile], inputs: Path) -> None:
    """Stage bounded uploads without following redirects or forwarding credentials."""
    size = 0
    async with httpx.AsyncClient(trust_env=False, follow_redirects=False, timeout=30) as client:
        for index, file in enumerate(files):
            url = urlsplit(file.download_url)
            host = url.hostname or ""
            if (
                url.scheme != "https"
                or not host.endswith(".oaiusercontent.com")
                or url.username is not None
                or url.password is not None
                or url.port not in (None, 443)
            ):
                raise ValueError("Files must use HTTPS download URLs from oaiusercontent.com")
            name = file.file_name or file.file_id
            if not name or name in (".", "..") or Path(name).name != name or PureWindowsPath(name).name != name:
                raise ValueError("File names must be base names without directories")
            if "\x00" in name:
                raise ValueError("File names must not contain null bytes")
            directory = inputs / str(index)
            directory.mkdir()
            try:
                async with client.stream("GET", file.download_url) as response:
                    response.raise_for_status()
                    with (directory / name).open("wb") as handle:
                        async for chunk in response.aiter_bytes(64 * 1024):
                            size += len(chunk)
                            if size > MAX_INPUT_BYTES:
                                raise ValueError("Uploads exceed the 100 MiB total limit")
                            handle.write(chunk)
            except httpx.HTTPError:
                # Signed download URLs must not appear in error messages or logs.
                raise ValueError(f"Could not download file {index + 1}; attach it again and retry") from None


def create_server(reports: Path, public_url: str, port: int = 8000) -> FastMCP:
    """Create a single-user development server; reports expire when it stops."""
    origin = urlsplit(public_url)
    if (
        origin.scheme not in ("http", "https")
        or (origin.scheme == "http" and origin.hostname not in ("localhost", "127.0.0.1"))
        or not origin.hostname
        or origin.username is not None
        or origin.password is not None
        or origin.path not in ("", "/")
        or origin.query
        or origin.fragment
    ):
        raise ValueError("public-url must be an HTTP or HTTPS origin without a path or credentials")
    public_url = public_url.rstrip("/")
    # HTTPX's informational request logs include signed attachment URLs.
    logging.getLogger("httpx").setLevel(logging.WARNING)
    server = FastMCP(
        "MultiQC",
        instructions=(
            "Run MultiQC on attached bioinformatics tool outputs. Return the report and data download links, "
            "then summarize the observed metrics. General statistics are raw values; the report may scale them "
            "for display. The preview can omit samples or metrics; use the full JSON for further analysis. "
            "Missing metrics are unknown. Do not invent QC thresholds or clinical conclusions."
        ),
        port=port,
        stateless_http=True,
        json_response=True,
        transport_security=TransportSecuritySettings(
            allowed_hosts=["127.0.0.1:*", "localhost:*", origin.netloc],
            allowed_origins=["http://127.0.0.1:*", "http://localhost:*", public_url],
        ),
    )
    jobs = asyncio.Semaphore(2)

    @server.tool(
        title="Create a MultiQC report",
        annotations=ToolAnnotations(readOnlyHint=True, destructiveHint=False, openWorldHint=True),
        meta={"openai/fileParams": ["files"]},
    )
    async def run_multiqc(
        files: Annotated[List[ChatGPTFile], Field(min_length=1, max_length=50)],
        title: Annotated[str, Field(max_length=200)] = "MultiQC report",
    ) -> Annotated[CallToolResult, ReportResult]:
        """Analyze attached QC logs or FastQC ZIP files and return report/data downloads plus raw metrics.

        Preserve original file names. Upload individual tool outputs, not an archive of an analysis directory.
        Reports remain available until the server stops. The statistics preview is limited to 20 samples and
        40 metrics per sample; the full parsed JSON contains the remaining data.
        """
        async with jobs:
            report_id = uuid4().hex
            output = reports / report_id
            output.mkdir()
            try:
                with TemporaryDirectory(prefix="multiqc-inputs-") as temporary:
                    inputs = Path(temporary) / "inputs"
                    inputs.mkdir()
                    try:
                        await asyncio.wait_for(download_files(files, inputs), timeout=120)
                    except asyncio.TimeoutError:
                        raise ValueError("File downloads exceeded the two minute deadline") from None
                    process = await asyncio.create_subprocess_exec(
                        sys.executable,
                        "-c",
                        WORKER,
                        str(inputs),
                        str(output),
                        title,
                        cwd=temporary,
                        stdout=asyncio.subprocess.DEVNULL,
                        stderr=asyncio.subprocess.PIPE,
                    )
                    try:
                        _, stderr = await asyncio.wait_for(process.communicate(), timeout=RUN_TIMEOUT)
                    except BaseException as error:
                        if process.returncode is None:
                            process.kill()
                        await process.wait()
                        if isinstance(error, asyncio.TimeoutError):
                            raise ValueError("MultiQC exceeded the five minute deadline") from None
                        raise
                    if process.returncode:
                        raise ValueError(
                            f"MultiQC could not create a report: {stderr.decode(errors='replace')[-2000:]}"
                        )
                if (
                    not (output / "multiqc_report.html").is_file()
                    or not (output / "multiqc_data" / "multiqc_data.json").is_file()
                ):
                    raise ValueError("MultiQC did not produce both the report and parsed JSON")
                summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
                stats = summary["general_stats"]
                preview = {
                    sample: dict(list(row.items())[:PREVIEW_METRICS])
                    for sample, row in list(stats.items())[:PREVIEW_SAMPLES]
                }
                result = ReportResult(
                    report_url=f"{public_url}/reports/{report_id}/multiqc_report.html",
                    data_url=f"{public_url}/reports/{report_id}/multiqc_data.json",
                    sample_count=summary["sample_count"],
                    modules=summary["modules"],
                    general_stats=preview,
                    general_stats_sample_count=len(stats),
                    general_stats_truncated=len(stats) > PREVIEW_SAMPLES
                    or any(len(row) > PREVIEW_METRICS for row in stats.values()),
                )
                return CallToolResult(
                    structuredContent=result.model_dump(),
                    content=[
                        TextContent(type="text", text=result.model_dump_json()),
                        ResourceLink(
                            type="resource_link",
                            name="multiqc_report.html",
                            uri=result.report_url,
                            mimeType="text/html",
                        ),
                        ResourceLink(
                            type="resource_link",
                            name="multiqc_data.json",
                            uri=result.data_url,
                            mimeType="application/json",
                        ),
                    ],
                )
            except BaseException:
                shutil.rmtree(output)
                raise

    @server.custom_route("/reports/{report_id}/{filename}", methods=["GET"])
    async def download_report(request: Request) -> Response:
        report_id = request.path_params["report_id"]
        filename = request.path_params["filename"]
        if len(report_id) != 32 or any(c not in "0123456789abcdef" for c in report_id):
            return Response(status_code=404)
        if filename not in ("multiqc_report.html", "multiqc_data.json"):
            return Response(status_code=404)
        path = reports / report_id / filename
        if filename == "multiqc_data.json":
            path = reports / report_id / "multiqc_data" / filename
        if not path.is_file():
            return Response(status_code=404)
        return FileResponse(
            path,
            filename=filename,
            headers={"Cache-Control": "no-store", "X-Content-Type-Options": "nosniff"},
        )

    return server


@click.command()
@click.option("--public-url", required=True, help="Public HTTPS origin of the tunnel or authenticated ingress.")
@click.option("--port", default=8000, type=click.IntRange(1, 65535), show_default=True)
def main(public_url: str, port: int) -> None:
    """Serve the MultiQC tool on localhost for a private ChatGPT development connection."""
    # ponytail: reports live until shutdown; add expiry for a long-running service.
    with TemporaryDirectory(prefix="multiqc-reports-") as temporary:
        create_server(Path(temporary), public_url, port).run(transport="streamable-http")


if __name__ == "__main__":
    main()
