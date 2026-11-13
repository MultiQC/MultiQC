import io
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

from multiqc import config, report
from multiqc.core import ai, ai_models
from multiqc.core.ai_models import parse_model

CASES = json.loads((Path(__file__).parent / "ai_model_cases.json").read_text())
JS_MODELS = Path(__file__).parent.parent / "multiqc/templates/default/src/js/ai-models.js"


class FakeResponse:
    def __init__(self, status: int, data: Any):
        self.status_code = status
        self.ok = 200 <= status < 300
        self._data = data
        self.text = json.dumps(data)

    def json(self):
        return self._data


class FakeApi:
    """Records requests and answers them from a queue of responses."""

    def __init__(self, *responses: FakeResponse, get_data: Any = None):
        self.responses = list(responses)
        self.requests: List[Dict[str, Any]] = []
        self.get_urls: List[str] = []
        self.get_data = get_data

    def post(self, url, headers=None, json=None, timeout=None):
        self.requests.append({"url": url, "headers": headers, "body": json})
        return self.responses.pop(0)

    def get(self, url, headers=None, timeout=None):
        self.get_urls.append(url)
        if self.get_data is None:
            return FakeResponse(404, {"error": {"message": "not found"}})
        return FakeResponse(200, self.get_data)


@pytest.fixture
def fake_api(monkeypatch):
    def install(*responses: FakeResponse, get_data: Any = None) -> FakeApi:
        api = FakeApi(*responses, get_data=get_data)
        monkeypatch.setattr(ai.requests, "post", api.post)
        monkeypatch.setattr(ai.requests, "get", api.get)
        monkeypatch.setattr(ai.time, "sleep", lambda _: None)
        return api

    return install


@pytest.fixture(autouse=True)
def quiet(monkeypatch):
    monkeypatch.setattr(config, "verbose", True)  # skip the progress spinner


def claude_response(*blocks: Dict[str, Any], stop_reason="end_turn", model="claude-opus-5-5") -> FakeResponse:
    return FakeResponse(200, {"model": model, "stop_reason": stop_reason, "content": list(blocks)})


def openai_response(content: Optional[str], finish_reason="stop", model="gpt-6.1-sol") -> FakeResponse:
    return FakeResponse(
        200, {"model": model, "choices": [{"message": {"content": content}, "finish_reason": finish_reason}]}
    )


TEXT = {"type": "text", "text": "Summary"}


# Model rules: Python and JS against the shared fixture


@pytest.mark.parametrize("case", CASES["models"], ids=lambda c: c["model"] or "<empty>")
def test_model_rules(case):
    info = parse_model(case["model"])
    assert info.vendor == case["vendor"]
    assert info.tier == case["tier"]
    assert list(info.version) == case["version"]
    assert ai_models.rule_context_window(info) == case["context_window"]
    assert ai_models.claude_thinking_style(info) == case["claude_thinking"]
    assert ai_models.claude_supports_server_fallbacks(info) == case["server_fallbacks"]
    assert ai_models.bedrock_api_version(info) == case["bedrock_api_version"]
    assert ai_models.claude_default_max_tokens(info) == case["claude_max_tokens"]


@pytest.mark.parametrize("case", CASES["context_overflow"], ids=lambda c: c["message"][:40])
def test_parse_context_overflow(case):
    overflow = ai_models.parse_context_overflow(case["message"])
    if case.get("overflow") is False:
        assert overflow is None
    else:
        assert overflow == ai_models.ContextOverflow(limit=case["limit"], actual=case["actual"])


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_js_model_rules_match_python():
    script = f"""
import * as m from {json.dumps(JS_MODELS.as_uri())};
const cases = {json.dumps(CASES)};
const models = cases.models.map((c) => {{
  const info = m.parseModel(c.model);
  return {{
    model: c.model, vendor: info.vendor, tier: info.tier, version: info.version,
    context_window: m.ruleContextWindow(info), claude_thinking: m.claudeThinkingStyle(info),
    server_fallbacks: m.claudeSupportsServerFallbacks(info), bedrock_api_version: m.bedrockApiVersion(info),
    claude_max_tokens: m.claudeDefaultMaxTokens(info),
  }};
}});
const overflow = cases.context_overflow.map((c) => m.parseContextOverflow(c.message));
console.log(JSON.stringify({{ models, overflow }}));
"""
    result = subprocess.run(["node", "--input-type=module", "-e", script], capture_output=True, text=True, check=True)
    output = json.loads(result.stdout)
    assert output["models"] == CASES["models"]
    for case, overflow in zip(CASES["context_overflow"], output["overflow"]):
        if case.get("overflow") is False:
            assert overflow is None
        else:
            assert overflow == {"limit": case["limit"], "actual": case["actual"]}


def test_shrink_budget():
    # Provider counted 1.5x MultiQC's estimate: scale the budget to fit, with a margin
    overflow = ai_models.ContextOverflow(limit=1_000_000, actual=1_200_000)
    assert ai_models.shrink_budget(900_000, 800_000, overflow) == int(1_000_000 / 1.5 * 0.95)
    # Limit only
    assert ai_models.shrink_budget(900_000, 800_000, ai_models.ContextOverflow(500_000, None)) == 450_000
    # Nothing known: halve
    assert ai_models.shrink_budget(900_000, 800_000, ai_models.ContextOverflow(None, None)) == 450_000
    # Always shrinks
    assert ai_models.shrink_budget(100, 10, ai_models.ContextOverflow(1_000_000, 10)) == 99


def test_models_endpoint_helpers():
    assert ai_models.models_url_for_endpoint("http://host:8000/v1/chat/completions") == "http://host:8000/v1/models"
    assert ai_models.models_url_for_endpoint("https://gw/v1/messages") == "https://gw/v1/models"
    assert ai_models.models_url_for_endpoint("https://gw/llm/invoke") == "https://gw/llm/models"
    vllm = {"data": [{"id": "other", "max_model_len": 1}, {"id": "qwen3", "max_model_len": 131072}]}
    assert ai_models.find_context_window(ai_models.find_model_entry(vllm, "qwen3")) == 131072
    openrouter = {"data": [{"id": "x/y", "context_length": 400000, "top_provider": {"context_length": 200000}}]}
    assert ai_models.find_context_window(ai_models.find_model_entry(openrouter, "x/y")) == 400000
    litellm = {"data": [{"model_name": "z", "id": "z", "model_info": {"max_input_tokens": 64000}}]}
    assert ai_models.find_context_window(ai_models.find_model_entry(litellm, "z")) == 64000
    assert ai_models.find_context_window(ai_models.find_model_entry({"data": [{"id": "a"}]}, "a")) is None
    assert ai_models.find_model_entry(None, "a") is None


# Anthropic


def test_anthropic_default_request(fake_api):
    api = fake_api(claude_response(TEXT), get_data={"id": "claude-opus-5-5", "max_input_tokens": 1_000_000})
    client = ai.AnthropicClient("key")
    response = client.query("prompt")

    assert response.content == "Summary"
    request = api.requests[0]
    assert request["url"] == "https://api.anthropic.com/v1/messages"
    assert request["headers"]["anthropic-version"] == "2026-10-21"
    assert request["headers"]["anthropic-beta"] == ai_models.ANTHROPIC_FALLBACK_BETA
    assert request["headers"]["x-api-key"] == "key"
    assert request["body"] == {
        "model": "claude-opus-5-5",
        "max_tokens": 16000,
        "messages": [{"role": "user", "content": "prompt"}],
        "fallbacks": "default",
    }
    assert api.get_urls == ["https://api.anthropic.com/v1/models/claude-opus-5-5"]
    assert client.context_window() == 1_000_000


def test_anthropic_model_from_config(fake_api):
    config.ai_model = "claude-haiku-5-5"
    api = fake_api(claude_response(TEXT))
    ai.AnthropicClient("key").query("prompt")
    body = api.requests[0]["body"]
    assert body["model"] == "claude-haiku-5-5"
    assert "fallbacks" not in body  # Haiku has no server-side fallback
    assert "anthropic-beta" not in api.requests[0]["headers"]


def test_anthropic_ignores_non_claude_model():
    config.ai_model = "gpt-6.1-sol"
    assert ai.AnthropicClient("key").model == ai_models.DEFAULT_MODELS["anthropic"]


@pytest.mark.parametrize(
    "model, thinking, max_tokens",
    [
        ("claude-opus-5-5", {"type": "adaptive"}, 16000),
        ("claude-sonnet-4-6", {"type": "adaptive"}, 16000),
        ("claude-sonnet-4-5", {"type": "enabled", "budget_tokens": 20000}, 24096),
        ("claude-3-5-sonnet-20241022", None, 4096),
    ],
)
def test_anthropic_extended_thinking(fake_api, model, thinking, max_tokens):
    config.ai_model = model
    config.ai_extended_thinking = True
    config.ai_thinking_budget_tokens = 20000
    api = fake_api(claude_response(TEXT))
    ai.AnthropicClient("key").query("prompt")
    body = api.requests[0]["body"]
    assert body.get("thinking") == thinking
    assert body["max_tokens"] == max_tokens


def test_anthropic_effort_and_output_tokens(fake_api):
    config.ai_reasoning_effort = "xhigh"
    config.ai_max_completion_tokens = 32000
    api = fake_api(claude_response(TEXT), get_data={"id": "claude-opus-5-5", "max_tokens": 20000})
    ai.AnthropicClient("key").query("prompt")
    body = api.requests[0]["body"]
    assert body["output_config"] == {"effort": "xhigh"}
    assert body["max_tokens"] == 20000  # capped at the model's maximum from the models API


def test_anthropic_response_with_thinking_and_fallback(fake_api):
    fake_api(
        claude_response(
            {"type": "thinking", "thinking": ""},
            {"type": "text", "text": "Partial answer from the declining model"},
            {"type": "fallback", "from": {"model": "claude-opus-5-5"}, "to": {"model": "claude-opus-4-8"}},
            {"type": "thinking", "thinking": ""},
            {"type": "text", "text": "Final "},
            {"type": "text", "text": "answer"},
            model="claude-opus-4-8",
        )
    )
    response = ai.AnthropicClient("key").query("prompt")
    assert response.content == "Final answer"
    assert response.model == "claude-opus-4-8"


def test_anthropic_refusal(fake_api):
    fake_api(claude_response(stop_reason="refusal"))
    assert ai.AnthropicClient("key").query("prompt").content == ai.REFUSAL_STUB_CONTENT


def test_anthropic_out_of_output_tokens(fake_api):
    fake_api(claude_response({"type": "thinking", "thinking": ""}, stop_reason="max_tokens"))
    with pytest.raises(ValueError, match="ai_max_completion_tokens"):
        ai.AnthropicClient("key").query("prompt")


def test_rejected_fields_are_dropped(fake_api):
    """Fields a model rejects are removed one at a time, and the beta header goes with `fallbacks`."""
    config.ai_extra_query_options = {"temperature": 0.0, "model": "claude-opus-5-5"}
    api = fake_api(
        FakeResponse(
            400,
            {
                "type": "error",
                "error": {"type": "invalid_request_error", "message": "`temperature` is deprecated for this model."},
            },
        ),
        FakeResponse(
            400,
            {
                "type": "error",
                "error": {"type": "invalid_request_error", "message": "fallbacks: Extra inputs are not permitted"},
            },
        ),
        claude_response(TEXT),
    )
    assert ai.AnthropicClient("key").query("prompt").content == "Summary"
    assert "temperature" in api.requests[0]["body"]
    assert "temperature" not in api.requests[1]["body"]
    assert "fallbacks" in api.requests[1]["body"]
    assert "fallbacks" not in api.requests[2]["body"]
    assert "anthropic-beta" not in api.requests[2]["headers"]


def test_max_tokens_reduced_to_model_limit(fake_api):
    api = fake_api(
        FakeResponse(
            400,
            {"error": {"message": "max_tokens: 16000 > 8192, which is the maximum allowed number of output tokens"}},
        ),
        claude_response(TEXT),
    )
    ai.AnthropicClient("key").query("prompt")
    assert api.requests[1]["body"]["max_tokens"] == 8192


def test_unrepairable_400_is_not_retried(fake_api):
    api = fake_api(FakeResponse(400, {"error": {"message": "messages.0.content: Input should be a valid list"}}))
    with pytest.raises(ai.ProviderError):
        ai.AnthropicClient("key").query("prompt")
    assert len(api.requests) == 1


def test_transient_errors_are_retried(fake_api):
    api = fake_api(
        FakeResponse(429, {"error": {"message": "rate limited"}}),
        FakeResponse(529, {"error": {"message": "overloaded"}}),
        claude_response(TEXT),
    )
    assert ai.AnthropicClient("key").query("prompt").content == "Summary"
    assert len(api.requests) == 3


def test_missing_model(fake_api):
    fake_api(FakeResponse(404, {"type": "error", "error": {"type": "not_found_error", "message": "model: claude-x"}}))
    with pytest.raises(ai.ModelUnavailableError):
        ai.AnthropicClient("key").query("prompt")


def test_context_overflow_is_raised(fake_api):
    fake_api(FakeResponse(400, {"error": {"message": "prompt is too long: 1050000 tokens > 1000000 maximum"}}))
    with pytest.raises(ai.ContextOverflowError) as e:
        ai.AnthropicClient("key").query("prompt")
    assert e.value.overflow == ai_models.ContextOverflow(limit=1_000_000, actual=1_050_000)


# OpenAI and custom endpoints


def test_openai_default_request(fake_api):
    api = fake_api(openai_response("Summary"))
    client = ai.OpenAiClient("key")
    assert client.query("prompt").content == "Summary"
    request = api.requests[0]
    assert request["headers"]["Authorization"] == "Bearer key"
    assert request["body"] == {"model": "gpt-6.1-sol", "messages": [{"role": "user", "content": "prompt"}]}
    assert client.context_window() == 922_000
    assert api.get_urls == []  # OpenAI's models endpoint doesn't report context windows


def test_openai_configured_reasoning(fake_api):
    config.ai_model = "o3"
    config.ai_reasoning_effort = "high"
    config.ai_max_completion_tokens = 50000
    api = fake_api(openai_response("Summary"))
    ai.OpenAiClient("key").query("prompt")
    body = api.requests[0]["body"]
    assert body["reasoning_effort"] == "high"
    assert body["max_completion_tokens"] == 50000
    assert "temperature" not in body


def test_openai_rejected_param(fake_api):
    config.ai_model = "gpt-4o"
    config.ai_reasoning_effort = "high"
    api = fake_api(
        FakeResponse(400, {"error": {"message": "Unsupported parameter", "param": "reasoning_effort"}}),
        openai_response("Summary"),
    )
    ai.OpenAiClient("key").query("prompt")
    assert "reasoning_effort" not in api.requests[1]["body"]


def test_openai_empty_answer_after_running_out_of_tokens(fake_api):
    fake_api(openai_response(None, finish_reason="length"))
    with pytest.raises(ValueError, match="ai_max_completion_tokens"):
        ai.OpenAiClient("key").query("prompt")


def test_custom_endpoint_probes_context_window(fake_api, monkeypatch):
    config.ai_provider = "custom"
    config.ai_summary = True
    config.ai_model = "qwen3"
    config.ai_custom_endpoint = "http://localhost:8000/v1/chat/completions"
    monkeypatch.setenv("OPENAI_API_KEY", "key")
    api = fake_api(openai_response("Summary"), get_data={"data": [{"id": "qwen3", "max_model_len": 131072}]})
    client = ai.get_llm_client()
    assert isinstance(client, ai.OpenAiClient)
    assert client.context_window() == 131072
    assert api.get_urls == ["http://localhost:8000/v1/models"]
    client.query("prompt")
    assert api.requests[0]["body"] == {"model": "qwen3", "messages": [{"role": "user", "content": "prompt"}]}


def test_custom_endpoint_unknown_context_window(fake_api, monkeypatch):
    config.ai_model = "llama3.1:8b"
    fake_api(get_data={"data": [{"id": "llama3.1:8b"}]})
    client = ai.OpenAiClient("key", endpoint="http://localhost:11434/v1/chat/completions")
    assert client.context_window() == ai_models.FALLBACK_CONTEXT_WINDOW
    config.ai_custom_context_window = 32000
    assert ai.OpenAiClient("key", endpoint="http://x/v1/chat/completions").context_window() == 32000


def test_custom_gateway_headers_and_response_path(fake_api, monkeypatch):
    """An internal gateway with its own auth header and response shape (MultiQC#3161)."""
    config.ai_provider = "custom"
    config.ai_summary = True
    config.ai_model = "gateway-model"
    config.ai_custom_endpoint = "https://gateway.example.com/llm/invoke"
    config.ai_custom_headers = {"x-api-key": "${GATEWAY_TOKEN}", "Authorization": None}
    config.ai_custom_response_path = "completion.content"
    config.ai_extra_query_options = {"model": None}
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setenv("GATEWAY_TOKEN", "secret")
    api = fake_api(FakeResponse(200, {"completion": {"content": "Summary"}}))

    client = ai.get_llm_client()
    assert client is not None
    assert client.query("prompt").content == "Summary"
    request = api.requests[0]
    assert request["headers"] == {"Content-Type": "application/json", "x-api-key": "secret"}
    assert request["body"] == {"messages": [{"role": "user", "content": "prompt"}]}


def test_custom_headers_missing_env_var(monkeypatch):
    config.ai_provider = "custom"
    config.ai_summary = True
    config.ai_model = "m"
    config.ai_custom_endpoint = "https://gateway.example.com/v1/chat/completions"
    config.ai_custom_headers = {"x-api-key": "${NOT_SET_ANYWHERE}"}
    monkeypatch.delenv("NOT_SET_ANYWHERE", raising=False)
    with pytest.raises(RuntimeError, match="NOT_SET_ANYWHERE"):
        ai.get_llm_client()


def test_custom_anthropic_style_endpoint(fake_api, monkeypatch):
    config.ai_provider = "custom"
    config.ai_summary = True
    config.ai_model = "claude-sonnet-5-5"
    config.ai_custom_endpoint = "https://gateway.example.com/v1/messages"
    config.ai_custom_api_style = "anthropic"
    monkeypatch.setenv("ANTHROPIC_API_KEY", "key")
    api = fake_api(claude_response(TEXT))
    client = ai.get_llm_client()
    assert isinstance(client, ai.AnthropicClient)
    client.query("prompt")
    request = api.requests[0]
    assert request["url"] == "https://gateway.example.com/v1/messages"
    assert request["headers"]["x-api-key"] == "key"
    assert "fallbacks" not in request["body"]  # only sent to the Claude API itself


class EchoClient(ai.Client):
    name = "echo"
    title = "Echo"
    model = "echo-1"

    def query(self, prompt: str) -> ai.LlmResponse:
        return ai.LlmResponse(content=f"{len(prompt)} characters", model=self.model)


def test_custom_client_class():
    config.ai_provider = "custom"
    config.ai_summary = True
    config.ai_custom_client = "tests.test_ai:EchoClient"
    client = ai.get_llm_client()
    assert isinstance(client, EchoClient)
    assert client.interpret_report_short("data").interpretation.summary.endswith("characters")
    assert client.context_window() == ai_models.FALLBACK_CONTEXT_WINDOW


def test_custom_client_class_must_subclass_client():
    config.ai_provider = "custom"
    config.ai_summary = True
    config.ai_custom_client = "tests.test_ai:FakeApi"
    with pytest.raises(RuntimeError, match="subclass"):
        ai.get_llm_client()


# AWS Bedrock


def bedrock_client(monkeypatch, model: str, *responses: Any):
    boto3 = pytest.importorskip("boto3")
    from unittest.mock import MagicMock

    monkeypatch.setattr(config, "ai_model", model)
    monkeypatch.setattr(boto3, "client", lambda **kwargs: MagicMock())
    client = ai.AWSBedrockClient()
    client.client.invoke_model.side_effect = list(responses)
    return client


def bedrock_body(*blocks, stop_reason="end_turn"):
    return {"body": io.BytesIO(json.dumps({"stop_reason": stop_reason, "content": list(blocks)}).encode())}


@pytest.mark.parametrize(
    "model, api_version",
    [
        ("us.anthropic.claude-sonnet-4-6", "bedrock-2023-05-31"),
        ("anthropic.claude-opus-4-7", "2026-10-21"),
        ("global.anthropic.claude-sonnet-5-5", "2026-10-21"),
    ],
)
def test_bedrock_request(monkeypatch, model, api_version):
    client = bedrock_client(monkeypatch, model, bedrock_body({"type": "thinking", "thinking": ""}, TEXT))
    assert client.query("prompt").content == "Summary"
    body = json.loads(client.client.invoke_model.call_args.kwargs["body"])
    assert body["anthropic_version"] == api_version
    assert "thinking" not in body  # each model runs with its own default
    assert "fallbacks" not in body
    assert "model" not in body


def test_bedrock_end_of_life(monkeypatch):
    botocore = pytest.importorskip("botocore.exceptions")
    error = botocore.ClientError(
        {
            "Error": {
                "Code": "ResourceNotFoundException",
                "Message": "This model version has reached the end of its life.",
            },
            "ResponseMetadata": {"HTTPStatusCode": 404},
        },
        "InvokeModel",
    )
    client = bedrock_client(monkeypatch, "anthropic.claude-3-haiku-20240307-v1:0", error)
    with pytest.raises(ai.ModelUnavailableError):
        client.query("prompt")


# Prompt building and the summary flow


class CountingClient(EchoClient):
    """Counts one token per character, with a small context window."""

    window = 10_000
    calls: List[str]

    def __init__(self):
        super().__init__()
        self.calls = []

    def n_tokens(self, text: str) -> int:
        return len(text)

    def context_window(self) -> int:
        return self.window

    def output_tokens(self) -> int:
        return 0


def metadata_with_sections(sizes: Dict[str, int]) -> ai.AiReportMetadata:
    return ai.AiReportMetadata(
        tools={"tool": ai.AiToolMetadata(name="Tool")},
        sections={
            name: ai.AiSectionMetadata(name=name, module_anchor="tool", content="x" * size)
            for name, size in sizes.items()
        },
    )


def test_build_prompt_fits_sections_greedily():
    """A section that doesn't fit is skipped, but later sections that do fit are still included."""
    client = CountingClient()
    metadata = metadata_with_sections({"small": 1000, "huge": 50_000, "medium": 3000})
    build = ai.build_prompt(metadata, "system", client, budget=6000)
    assert not build.too_long
    assert "x" * 1000 in build.prompt
    assert "x" * 3000 in build.prompt
    assert "x" * 50_000 not in build.prompt
    assert "left out to fit the model's context window: Tool: huge" in build.prompt
    # The estimate is taken before blank lines are collapsed, so it never undercounts
    assert len("system") + len(build.prompt) <= build.n_tokens <= 6000


def test_build_prompt_without_client_includes_everything():
    metadata = metadata_with_sections({"a": 100, "b": 200})
    build = ai.build_prompt(metadata, "system")
    assert "x" * 100 in build.prompt and "x" * 200 in build.prompt
    assert "left out" not in build.prompt


def test_summary_retries_with_less_data_on_overflow(monkeypatch):
    config.ai_summary = True
    client = CountingClient()
    budgets: List[int] = []

    def interpret(prompt: str):
        budgets.append(len(prompt))
        if len(budgets) == 1:
            raise ai.ContextOverflowError(ai_models.ContextOverflow(limit=100, actual=200), "prompt is too long")
        return ai.InterpretationResponse(interpretation=ai.InterpretationOutput(summary="Fits now"), model="m")

    monkeypatch.setattr(ai, "get_llm_client", lambda: client)
    monkeypatch.setattr(client, "interpret_report_short", interpret)
    ai.add_ai_summary_to_report()
    assert len(budgets) == 2
    assert report.ai_global_summary == "Fits now"


def test_summary_stub_for_retired_model(monkeypatch):
    config.ai_summary = True
    client = CountingClient()

    def interpret(prompt: str):
        raise ai.ModelUnavailableError("Model 'echo-1' is not available")

    monkeypatch.setattr(ai, "get_llm_client", lambda: client)
    monkeypatch.setattr(client, "interpret_report_short", interpret)
    ai.add_ai_summary_to_report()
    assert "echo-1" in report.ai_global_summary
    assert "ai_model" in report.ai_global_summary


def test_js_config(monkeypatch):
    config.ai_summary = True
    client = CountingClient()
    monkeypatch.setattr(ai, "get_llm_client", lambda: client)
    ai.add_ai_summary_to_report()
    assert report.ai_js_config["defaultModels"] == ai_models.DEFAULT_MODELS
    assert report.ai_js_config["contextWindows"] == {"echo-1": 10_000}
    assert report.ai_js_config["anthropicApiVersion"] == ai_models.ANTHROPIC_API_VERSION


def test_rejected_field_is_not_sent_again(fake_api):
    config.ai_model = "qwen3"
    config.ai_reasoning_effort = "high"
    api = fake_api(
        FakeResponse(400, {"error": {"message": "reasoning_effort is not supported", "param": "reasoning_effort"}}),
        openai_response("First"),
        openai_response("Second"),
    )
    client = ai.OpenAiClient("key", endpoint="http://localhost:8000/v1/chat/completions")
    client.query("prompt")
    client.query("prompt")
    assert [("reasoning_effort" in r["body"]) for r in api.requests] == [True, False, False]


def test_payload_too_large_is_an_overflow(fake_api):
    fake_api(FakeResponse(413, {"error": {"message": "Request Entity Too Large"}}))
    with pytest.raises(ai.ContextOverflowError) as e:
        ai.OpenAiClient("key").query("prompt")
    assert e.value.overflow == ai_models.ContextOverflow(limit=None, actual=None)
