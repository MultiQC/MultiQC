"""Exercise the MCP contract, real report rendering, and upload boundaries."""

import asyncio
import logging
from pathlib import Path
from typing import List
from urllib.parse import urlsplit

import httpx
import multiqc_chatgpt as integration
import pytest
from jsonschema import validate
from mcp import ClientSession
from mcp.client.streamable_http import streamable_http_client
from mcp.shared.memory import create_connected_server_and_client_session

NANOQ = b"""Nanoq Read Summary
Number of reads: 100
Number of bases: 200000
N50 read length: 3000
Longest read: 10000
Shortest read: 500
Mean read length: 2000
Median read length: 1800
Mean read quality: 15
Median read quality: 14
Read length thresholds
> 1000 90 90.0%
> 5000 20 20.0%
Read quality thresholds
> 10 80 80.0%
> 20 20 20.0%
"""


def attachments(prefix: str, count: int = 1) -> List[dict]:
    return [
        {
            "download_url": f"https://files.oaiusercontent.com/{prefix}{index}?sig=private",
            "file_id": f"file-{prefix}{index}",
            "file_name": f"{prefix}{index}.txt",
        }
        for index in range(count)
    ]


def mock_downloads(monkeypatch, content=NANOQ, status=200):
    original = httpx.AsyncClient
    transport = httpx.MockTransport(lambda request: httpx.Response(status, content=content))
    monkeypatch.setattr(integration.httpx, "AsyncClient", lambda **kwargs: original(transport=transport, **kwargs))


def test_http_report_round_trip_and_isolated_runs(tmp_path: Path, monkeypatch, caplog):
    caplog.set_level(logging.INFO)
    server = integration.create_server(tmp_path, "http://localhost:8000")
    app = server.streamable_http_app()
    # Keep the actual ASGI client separate from the simulated OpenAI downloads.
    client = httpx.AsyncClient(transport=httpx.ASGITransport(app=app), base_url="http://localhost:8000")
    mock_downloads(monkeypatch)

    async def check():
        async with server.session_manager.run(), client:
            async with streamable_http_client("http://localhost:8000/mcp", http_client=client) as (read, write, _):
                async with ClientSession(read, write) as session:
                    await session.initialize()
                    tool = (await session.list_tools()).tools[0]
                    assert tool.name == "run_multiqc"
                    assert tool.meta is not None
                    assert tool.outputSchema is not None
                    assert tool.meta["openai/fileParams"] == ["files"]
                    file_schema = tool.inputSchema["$defs"]["ChatGPTFile"]
                    assert set(file_schema["required"]) == {"download_url", "file_id"}
                    assert set(file_schema["properties"]) == {"download_url", "file_id", "mime_type", "file_name"}
                    large, small = await asyncio.gather(
                        session.call_tool("run_multiqc", {"files": attachments("large", 21)}),
                        session.call_tool("run_multiqc", {"files": attachments("small"), "title": "<QC & audit>"}),
                    )
                    assert not large.isError, large
                    assert not small.isError, small
                    assert large.structuredContent is not None
                    assert small.structuredContent is not None
                    for result in (large, small):
                        validate(result.structuredContent, tool.outputSchema)
                        assert [part.type for part in result.content] == ["text", "resource_link", "resource_link"]
                    assert large.structuredContent["sample_count"] == 21
                    assert large.structuredContent["general_stats_sample_count"] == 21
                    assert large.structuredContent["general_stats_truncated"]
                    assert len(large.structuredContent["general_stats"]) == 20
                    output = small.structuredContent
                    assert output["sample_count"] == 1
                    assert output["modules"] == ["nanoq"]
                    assert not output["general_stats_truncated"]
                    assert set(output["general_stats"]) == {"small0"}
                    assert output["general_stats"]["small0"]["nanoq.Number of reads"] == 100
                    assert output["report_url"] != large.structuredContent["report_url"]
                    html = await client.get(output["report_url"])
                    assert html.status_code == 200
                    assert "attachment" in html.headers["content-disposition"]
                    assert "&lt;QC &amp; audit&gt;" in html.text
                    assert "small0" in html.text
                    data = await client.get(output["data_url"])
                    assert data.status_code == 200
                    parsed = data.json()
                    assert parsed["report_general_stats_data"]["nanoq"]["small0"]["Number of reads"] == 100
                    assert (await client.get("/reports/invalid/multiqc_report.html")).status_code == 404
                    assert (
                        await client.get(
                            urlsplit(output["report_url"]).path.replace("multiqc_report.html", "summary.json")
                        )
                    ).status_code == 404

    asyncio.run(asyncio.wait_for(check(), timeout=60))
    assert len(list(tmp_path.iterdir())) == 2
    assert "sig=private" not in caplog.text


@pytest.mark.parametrize(
    "changes",
    [
        {"download_url": "http://files.oaiusercontent.com/file"},
        {"download_url": "https://127.0.0.1/file"},
        {"download_url": "https://files.oaiusercontent.com.evil.example/file"},
        {"download_url": "https://user:secret@files.oaiusercontent.com/file"},
        {"download_url": "https://files.oaiusercontent.com:8080/file"},
        {"file_name": "../escape.txt"},
        {"file_name": "C:\\escape.txt"},
        {"file_name": "bad\x00.txt"},
    ],
)
def test_reject_unsafe_uploads(tmp_path: Path, monkeypatch, changes):
    mock_downloads(monkeypatch)
    files = attachments("unsafe")
    files[0].update(changes)

    async def check():
        server = integration.create_server(tmp_path, "http://localhost:8000")
        async with create_connected_server_and_client_session(server) as session:
            result = await session.call_tool("run_multiqc", {"files": files})
            assert result.isError
            assert result.structuredContent is None

    asyncio.run(check())
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("failure", ["size", "redirect", "http_error", "no_analysis", "timeout", "empty", "too_many"])
def test_failed_runs_leave_no_reports(tmp_path: Path, monkeypatch, failure):
    files = attachments("failure")
    content, status = NANOQ, 200
    if failure == "size":
        monkeypatch.setattr(integration, "MAX_INPUT_BYTES", 10)
    elif failure == "redirect":
        status = 302
    elif failure == "http_error":
        status = 404
    elif failure == "no_analysis":
        content = b"Unrecognized analysis output\n"
    elif failure == "timeout":
        monkeypatch.setattr(integration, "RUN_TIMEOUT", 0.001)
    elif failure == "empty":
        files = []
    elif failure == "too_many":
        files = attachments("failure", 51)
    mock_downloads(monkeypatch, content, status)

    async def check():
        server = integration.create_server(tmp_path, "http://localhost:8000")
        async with create_connected_server_and_client_session(server) as session:
            result = await session.call_tool("run_multiqc", {"files": files})
            assert result.isError
            assert result.structuredContent is None
            if failure in ("redirect", "http_error"):
                assert "sig=private" not in str(result.content)

    asyncio.run(check())
    assert not list(tmp_path.iterdir())
