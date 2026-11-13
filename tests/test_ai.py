from unittest.mock import MagicMock

import pytest

from multiqc import config
from multiqc.core import ai

botocore_exceptions = pytest.importorskip("botocore.exceptions")
pytest.importorskip("boto3")


def _bedrock_client(monkeypatch, invoke_model_error: Exception) -> ai.AWSBedrockClient:
    monkeypatch.setattr(config, "ai_model", "us.anthropic.claude-sonnet-4-6")
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-west-2")
    client = ai.AWSBedrockClient()
    client.client = MagicMock()
    client.client.invoke_model.side_effect = invoke_model_error
    return client


def test_bedrock_model_end_of_life_returns_stub(monkeypatch):
    error = botocore_exceptions.ClientError(
        {
            "Error": {
                "Code": "ResourceNotFoundException",
                "Message": "This model version has reached the end of its life. Please refer to the AWS documentation for more details.",
            }
        },
        "InvokeModel",
    )
    client = _bedrock_client(monkeypatch, error)

    response = client.interpret_report_short("report")
    assert response.interpretation is not None
    assert response.interpretation.summary == ai.MODEL_END_OF_LIFE_STUB_CONTENT

    response = client.interpret_report_full("report")
    assert response.interpretation is not None
    assert response.interpretation.summary == ai.MODEL_END_OF_LIFE_STUB_CONTENT


def test_bedrock_other_client_errors_raise(monkeypatch):
    error = botocore_exceptions.ClientError(
        {"Error": {"Code": "ResourceNotFoundException", "Message": "Model not found"}},
        "InvokeModel",
    )
    client = _bedrock_client(monkeypatch, error)

    with pytest.raises(botocore_exceptions.ClientError):
        client.interpret_report_short("report")
