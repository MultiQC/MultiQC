from unittest.mock import MagicMock, patch

from multiqc import config
from multiqc.core.ai import AnthropicClient, is_reasoning_model


def test_is_reasoning_model_claude_5():
    assert is_reasoning_model("claude-sonnet-5-5")
    assert is_reasoning_model("claude-opus-5-5")
    assert is_reasoning_model("claude-haiku-5-5")
    assert is_reasoning_model("claude-sonnet-4-5")
    assert not is_reasoning_model("claude-3-haiku-20240307")
    assert not is_reasoning_model("gpt-4o")


def test_anthropic_client_default_model(monkeypatch):
    monkeypatch.setattr(config, "ai_model", None)
    client = AnthropicClient(api_key="test-key")
    assert client.model == "claude-sonnet-5-5"


def test_anthropic_client_custom_model(monkeypatch):
    monkeypatch.setattr(config, "ai_model", "claude-3-haiku-20240307")
    client = AnthropicClient(api_key="test-key")
    assert client.model == "claude-3-haiku-20240307"


def test_anthropic_client_query_omits_temperature_for_claude_5(monkeypatch):
    monkeypatch.setattr(config, "ai_model", "claude-sonnet-5-5")
    client = AnthropicClient(api_key="test-key")

    mock_response = {
        "content": [{"text": "Sample summary"}],
        "model": "claude-sonnet-5-5",
    }
    with patch.object(client, "_request_with_error_handling_and_retries", return_value=mock_response) as mock_request:
        result = client._query("Summarize this")
        assert result.content == "Sample summary"
        assert result.model == "claude-sonnet-5-5"

        mock_request.assert_called_once()
        body = mock_request.call_args[1]["body"]
        assert body["model"] == "claude-sonnet-5-5"
        assert "temperature" not in body


def test_anthropic_client_query_includes_temperature_for_legacy_claude(monkeypatch):
    monkeypatch.setattr(config, "ai_model", "claude-3-haiku-20240307")
    client = AnthropicClient(api_key="test-key")

    mock_response = {
        "content": [{"text": "Sample summary"}],
        "model": "claude-3-haiku-20240307",
    }
    with patch.object(client, "_request_with_error_handling_and_retries", return_value=mock_response) as mock_request:
        result = client._query("Summarize this")
        assert result.content == "Sample summary"

        mock_request.assert_called_once()
        body = mock_request.call_args[1]["body"]
        assert body["model"] == "claude-3-haiku-20240307"
        assert body.get("temperature") == 0.0
