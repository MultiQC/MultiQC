import io
import json
from typing import Any, Dict, Optional
from unittest.mock import MagicMock

import pytest

from multiqc import config
from multiqc.core import ai

pytest.importorskip("boto3")


def _bedrock_client(monkeypatch, model: str) -> ai.AWSBedrockClient:
    monkeypatch.setattr(config, "ai_model", model)
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-west-2")
    client = ai.AWSBedrockClient()
    client.client = MagicMock()
    response_body = {"stop_reason": "end_turn", "content": [{"type": "text", "text": "summary"}]}
    client.client.invoke_model.return_value = {"body": io.BytesIO(json.dumps(response_body).encode())}
    return client


@pytest.mark.parametrize(
    "model, expected_thinking",
    [
        ("us.anthropic.claude-sonnet-4-6", {"type": "disabled"}),
        ("us.anthropic.claude-sonnet-5-5", {"type": "between_tools"}),
        ("anthropic.claude-opus-5-5", None),
        ("global.anthropic.claude-fable-5-1", None),
    ],
)
def test_bedrock_thinking_disabled_per_model(
    monkeypatch, model: str, expected_thinking: Optional[Dict[str, Any]]
) -> None:
    monkeypatch.setattr(config, "ai_extended_thinking", False)
    client = _bedrock_client(monkeypatch, model)

    response = client.interpret_report_short("report")
    assert response.interpretation is not None
    assert response.interpretation.summary == "summary"

    body = json.loads(client.client.invoke_model.call_args.kwargs["body"])
    assert body.get("thinking") == expected_thinking
