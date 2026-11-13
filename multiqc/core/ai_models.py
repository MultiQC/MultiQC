"""
Model rules for AI summaries.

MultiQC never keeps a list of individual model names. Instead it reads the vendor, tier and
version out of a model ID and applies "at least version X" rules, so that new releases in a
known family work without code changes. Anything that is not recognised gets the most
conservative treatment: a minimal request and a modest context window.

The same rules are implemented for the browser in
`multiqc/templates/default/src/js/ai-models.js`. Both are tested against
`tests/data/ai_model_cases.json`, so keep them in sync.
"""

import logging
import re
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Tuple

logger = logging.getLogger(__name__)

# Defaults used when `ai_model` is not set. Bump these when a provider releases a new generation.
DEFAULT_MODELS: Dict[str, str] = {
    "anthropic": "claude-opus-5-5",
    "openai": "gpt-6.1-sol",
}

# Context window assumed for models that MultiQC can't recognise and whose endpoint doesn't report one.
FALLBACK_CONTEXT_WINDOW = 128_000

ANTHROPIC_API_VERSION = "2026-10-21"
ANTHROPIC_FALLBACK_BETA = "server-side-fallback-2026-07-01"
# Bedrock accepts the new API version only for Claude 4.7 and later.
BEDROCK_LEGACY_API_VERSION = "bedrock-2023-05-31"


@dataclass(frozen=True)
class ModelInfo:
    """What MultiQC could work out from a model ID."""

    vendor: Optional[str]  # "anthropic", "openai", or None if not recognised
    tier: str  # Claude: "opus", "sonnet", "haiku", "fable", ...; OpenAI: "gpt" or "o"
    version: Tuple[int, int]  # (major, minor); (0, 0) if unknown
    variant: str = ""  # anything after the version, e.g. "-sol", "o" (for gpt-4o), "-mini"

    def at_least(self, major: int, minor: int = 0) -> bool:
        return self.version >= (major, minor)


UNKNOWN_MODEL = ModelInfo(vendor=None, tier="", version=(0, 0))

# claude-opus-4-5, claude-sonnet-4-5-20250929, claude-opus-5, us.anthropic.claude-opus-4-8-v1:0,
# claude-opus-4-5@20251101 (Vertex), anthropic/claude-sonnet-4.5 (routers). The minor version is one
# or two digits, so a date suffix such as -20250514 is not mistaken for it.
_CLAUDE_RE = re.compile(r"(?<![a-z])claude-(?P<tier>[a-z]+)-(?P<major>\d{1,2})(?:[-.](?P<minor>\d{1,2}))?(?!\d)")
# Older style: claude-3-7-sonnet-latest, claude-3-haiku-20240307, anthropic/claude-3.5-sonnet
_CLAUDE_OLD_RE = re.compile(r"(?<![a-z])claude-(?P<major>\d)(?:[-.](?P<minor>\d))?-(?P<tier>[a-z]+)")
# o1, o3-mini, o4-mini-deep-research, openai/o3
_OPENAI_O_RE = re.compile(r"(?:^|/)o(?P<major>\d+)(?P<variant>(?:[-:@].*)?)$")
# gpt-4o, gpt-4.1-mini, gpt-5-chat-latest, gpt-6.1-sol, openai/gpt-5.5
_OPENAI_GPT_RE = re.compile(r"(?:^|/)gpt-(?P<major>\d+)(?:\.(?P<minor>\d+))?(?P<variant>.*)$")


def parse_model(model: Optional[str]) -> ModelInfo:
    """Work out the vendor, tier and version of a model from its ID."""
    if not model:
        return UNKNOWN_MODEL
    name = model.strip().lower()

    if m := _CLAUDE_RE.search(name):
        return ModelInfo("anthropic", m["tier"], (int(m["major"]), int(m["minor"] or 0)), name[m.end() :])
    if m := _CLAUDE_OLD_RE.search(name):
        return ModelInfo("anthropic", m["tier"], (int(m["major"]), int(m["minor"] or 0)), name[m.end() :])
    if "claude" in name:
        return ModelInfo("anthropic", "", (0, 0))
    if m := _OPENAI_O_RE.search(name):
        return ModelInfo("openai", "o", (int(m["major"]), 0), m["variant"])
    if m := _OPENAI_GPT_RE.search(name):
        return ModelInfo("openai", "gpt", (int(m["major"]), int(m["minor"] or 0)), m["variant"])
    return UNKNOWN_MODEL


def rule_context_window(info: ModelInfo) -> Optional[int]:
    """
    Input context window for a recognised model family, or None if MultiQC doesn't know.
    Used when the provider can't be asked. Too high is recovered from (the provider's error
    gives the real limit), too low just sends less data.
    """
    if info.vendor == "anthropic":
        if info.at_least(4, 6):
            return 1_000_000
        if info.at_least(3):
            return 200_000
        return None
    if info.vendor == "openai":
        if info.tier == "o":
            return 200_000
        if info.at_least(5, 5):
            return 922_000
        if info.at_least(5):
            return 272_000
        if info.version == (4, 1):
            return 1_000_000
        if info.at_least(4):
            return 128_000
        return None
    return None


def claude_thinking_style(info: ModelInfo) -> Optional[str]:
    """How a Claude model takes extended thinking: "adaptive" (4.6+), "budget" (3.7 to 4.5), or None."""
    if info.vendor != "anthropic":
        return None
    if info.at_least(4, 6):
        return "adaptive"
    if info.at_least(3, 7):
        return "budget"
    return None


def claude_default_max_tokens(info: ModelInfo) -> int:
    """
    Output token limit to request from a Claude model when `ai_max_completion_tokens` is not set.
    Thinking tokens count against it, so it needs headroom on models that think by default.
    """
    if info.vendor == "anthropic" and info.version != (0, 0) and not info.at_least(3, 7):
        return 4096
    return 16_000


def claude_supports_server_fallbacks(info: ModelInfo) -> bool:
    """
    Whether to ask the Claude API to re-run a refused request on another model. Offered for
    Claude 5 and later, except Haiku. If a model turns out not to accept it, the request is
    retried without it.
    """
    return info.vendor == "anthropic" and info.at_least(5) and info.tier != "haiku"


def bedrock_api_version(info: ModelInfo) -> str:
    """The `anthropic_version` body field for a Bedrock request."""
    if info.vendor == "anthropic" and info.at_least(4, 7):
        return ANTHROPIC_API_VERSION
    return BEDROCK_LEGACY_API_VERSION


# Fields that OpenAI-compatible servers use to report a model's context window, in order of preference.
_CONTEXT_WINDOW_KEYS = (
    "max_input_tokens",  # Anthropic models API, LiteLLM
    "max_model_len",  # vLLM
    "context_length",  # OpenRouter
    "max_context_length",  # LM Studio
    "context_window",
    "input_token_limit",
    "inputTokenLimit",  # Google
)


def find_context_window(entry: Any) -> Optional[int]:
    """Look for a context window field in a model description returned by a provider's models endpoint."""
    if not isinstance(entry, dict):
        return None
    candidates: Iterable[Dict[str, Any]] = [entry] + [v for v in entry.values() if isinstance(v, dict)]
    for key in _CONTEXT_WINDOW_KEYS:
        for candidate in candidates:
            value = candidate.get(key)
            if isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0:
                return int(value)
    return None


def find_model_entry(listing: Any, model: str) -> Optional[Dict[str, Any]]:
    """Find a model in a models endpoint response: an OpenAI-style `{"data": [...]}` listing, or a single model."""
    if not isinstance(listing, dict):
        return None
    if listing.get("id") == model:
        return listing
    for key in ("data", "models"):
        entries = listing.get(key)
        if isinstance(entries, list):
            for entry in entries:
                if isinstance(entry, dict) and model in (entry.get("id"), entry.get("name"), entry.get("model")):
                    return entry
    return None


def models_url_for_endpoint(endpoint: str) -> str:
    """Guess the models listing URL that sits next to a chat endpoint, e.g. .../v1/chat/completions -> .../v1/models"""
    base = endpoint.split("?")[0].rstrip("/")
    for suffix in ("/chat/completions", "/completions", "/messages", "/responses"):
        if base.endswith(suffix):
            return base[: -len(suffix)] + "/models"
    return base.rsplit("/", 1)[0] + "/models"


@dataclass(frozen=True)
class ContextOverflow:
    """What a provider said when a prompt was too long for the model."""

    limit: Optional[int]  # the model's input limit, if the error gave it
    actual: Optional[int]  # the prompt's size in the provider's tokens, if the error gave it


_OVERFLOW_PATTERNS = (
    # Anthropic, Bedrock: "prompt is too long: 1050000 tokens > 1000000 maximum"
    re.compile(r"(?P<actual>\d+)\s*tokens?\s*>\s*(?P<limit>\d+)"),
    # Anthropic (older): "input length and `max_tokens` exceed context limit: 197000 + 8192 > 200000"
    re.compile(r"(?P<actual>\d+)\s*\+\s*(?P<reserve>\d+)\s*>\s*(?P<limit>\d+)"),
    # OpenAI, vLLM: "maximum context length is 128000 tokens. However, your messages resulted in 130000 tokens"
    re.compile(
        r"maximum context length is (?P<limit>\d+)\D+?(?:resulted in|requested|has|contains)\D*?(?P<actual>\d+)",
        re.IGNORECASE,
    ),
    # Google: "The input token count (1200000) exceeds the maximum number of tokens allowed (1048576)"
    re.compile(r"\((?P<actual>\d+)\)\s*exceeds the maximum number of tokens allowed\s*\((?P<limit>\d+)\)"),
    # Limit only
    re.compile(r"maximum context length is (?P<limit>\d+)", re.IGNORECASE),
)

_OVERFLOW_HINTS = re.compile(
    r"context[_ ]length|context window|context limit|prompt is too long|input is too long|too many tokens|"
    r"maximum context|exceeds the maximum number of tokens|reduce the length",
    re.IGNORECASE,
)


def parse_context_overflow(message: str) -> Optional[ContextOverflow]:
    """If an error message says the prompt was too long, return what it says about the limit."""
    if not _OVERFLOW_HINTS.search(message):
        return None
    text = message.replace(",", "")
    for pattern in _OVERFLOW_PATTERNS:
        if m := pattern.search(text):
            groups = m.groupdict()
            limit = int(groups["limit"]) - int(groups.get("reserve") or 0)
            actual = int(groups["actual"]) if groups.get("actual") else None
            return ContextOverflow(limit=limit, actual=actual)
    return ContextOverflow(limit=None, actual=None)


def shrink_budget(budget: int, estimated_tokens: int, overflow: ContextOverflow) -> int:
    """
    A smaller prompt budget after a context overflow. When the provider reported both its limit and
    the prompt's size, use the ratio of its token count to MultiQC's estimate, which also corrects
    for tokenizer differences. Otherwise halve the budget.
    """
    if overflow.limit and overflow.actual and estimated_tokens:
        tokens_per_estimated_token = overflow.actual / estimated_tokens
        new_budget = int(overflow.limit / tokens_per_estimated_token * 0.95)
    elif overflow.limit:
        new_budget = int(overflow.limit * 0.9)
    else:
        new_budget = budget // 2
    return max(1, min(new_budget, budget - 1))
