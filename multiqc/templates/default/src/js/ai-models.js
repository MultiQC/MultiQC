////////////////////////////////////////////////
// AI model rules
//
// A port of multiqc/core/ai_models.py, for requests made from the report. MultiQC doesn't keep a
// list of model names: it reads the vendor, tier and version from a model ID and applies
// "at least version X" rules. Both implementations are tested against tests/data/ai_model_cases.json.
////////////////////////////////////////////////

export const FALLBACK_CONTEXT_WINDOW = 128000;
export const ANTHROPIC_API_VERSION = "2026-10-21";
export const ANTHROPIC_FALLBACK_BETA = "server-side-fallback-2026-07-01";
export const BEDROCK_LEGACY_API_VERSION = "bedrock-2023-05-31";

const UNKNOWN_MODEL = { vendor: null, tier: "", version: [0, 0], variant: "" };

const CLAUDE_RE = /(?<![a-z])claude-([a-z]+)-(\d{1,2})(?:[-.](\d{1,2}))?(?!\d)/;
const CLAUDE_OLD_RE = /(?<![a-z])claude-(\d)(?:[-.](\d))?-([a-z]+)/;
const OPENAI_O_RE = /(?:^|\/)o(\d+)((?:[-:@].*)?)$/;
const OPENAI_GPT_RE = /(?:^|\/)gpt-(\d+)(?:\.(\d+))?(.*)$/;

export function parseModel(model) {
  if (!model) return UNKNOWN_MODEL;
  const name = model.trim().toLowerCase();
  let m;
  if ((m = CLAUDE_RE.exec(name))) {
    return {
      vendor: "anthropic",
      tier: m[1],
      version: [+m[2], +(m[3] || 0)],
      variant: name.slice(m.index + m[0].length),
    };
  }
  if ((m = CLAUDE_OLD_RE.exec(name))) {
    return {
      vendor: "anthropic",
      tier: m[3],
      version: [+m[1], +(m[2] || 0)],
      variant: name.slice(m.index + m[0].length),
    };
  }
  if (name.includes("claude")) return { vendor: "anthropic", tier: "", version: [0, 0], variant: "" };
  if ((m = OPENAI_O_RE.exec(name))) return { vendor: "openai", tier: "o", version: [+m[1], 0], variant: m[2] };
  if ((m = OPENAI_GPT_RE.exec(name))) {
    return { vendor: "openai", tier: "gpt", version: [+m[1], +(m[2] || 0)], variant: m[3] };
  }
  return UNKNOWN_MODEL;
}

function atLeast(info, major, minor = 0) {
  return info.version[0] > major || (info.version[0] === major && info.version[1] >= minor);
}

export function ruleContextWindow(info) {
  if (info.vendor === "anthropic") {
    if (atLeast(info, 4, 6)) return 1000000;
    if (atLeast(info, 3)) return 200000;
    return null;
  }
  if (info.vendor === "openai") {
    if (info.tier === "o") return 200000;
    if (atLeast(info, 5, 5)) return 922000;
    if (atLeast(info, 5)) return 272000;
    if (info.version[0] === 4 && info.version[1] === 1) return 1000000;
    if (atLeast(info, 4)) return 128000;
    return null;
  }
  return null;
}

export function claudeThinkingStyle(info) {
  if (info.vendor !== "anthropic") return null;
  if (atLeast(info, 4, 6)) return "adaptive";
  if (atLeast(info, 3, 7)) return "budget";
  return null;
}

export function claudeDefaultMaxTokens(info) {
  if (info.vendor === "anthropic" && !(info.version[0] === 0 && info.version[1] === 0) && !atLeast(info, 3, 7)) {
    return 4096;
  }
  return 16000;
}

export function claudeSupportsServerFallbacks(info) {
  return info.vendor === "anthropic" && atLeast(info, 5) && info.tier !== "haiku";
}

export function bedrockApiVersion(info) {
  return info.vendor === "anthropic" && atLeast(info, 4, 7) ? ANTHROPIC_API_VERSION : BEDROCK_LEGACY_API_VERSION;
}

const OVERFLOW_PATTERNS = [
  /(?<actual>\d+)\s*tokens?\s*>\s*(?<limit>\d+)/,
  /(?<actual>\d+)\s*\+\s*(?<reserve>\d+)\s*>\s*(?<limit>\d+)/,
  /maximum context length is (?<limit>\d+)\D+?(?:resulted in|requested|has|contains)\D*?(?<actual>\d+)/i,
  /\((?<actual>\d+)\)\s*exceeds the maximum number of tokens allowed\s*\((?<limit>\d+)\)/,
  /maximum context length is (?<limit>\d+)/i,
];

const OVERFLOW_HINTS =
  /context[_ ]length|context window|context limit|prompt is too long|input is too long|too many tokens|maximum context|exceeds the maximum number of tokens|reduce the length/i;

// If an error message says the prompt was too long, return {limit, actual} (either can be null)
export function parseContextOverflow(message) {
  if (!message || !OVERFLOW_HINTS.test(message)) return null;
  const text = message.replace(/,/g, "");
  for (const pattern of OVERFLOW_PATTERNS) {
    const m = pattern.exec(text);
    if (m) {
      const g = m.groups;
      return {
        limit: parseInt(g.limit) - parseInt(g.reserve || "0"),
        actual: g.actual ? parseInt(g.actual) : null,
      };
    }
  }
  return { limit: null, actual: null };
}

// A smaller prompt budget after a context overflow. See shrink_budget() in ai_models.py
export function shrinkBudget(budget, estimatedTokens, overflow) {
  let newBudget;
  if (overflow.limit && overflow.actual && estimatedTokens) {
    newBudget = Math.floor((overflow.limit / (overflow.actual / estimatedTokens)) * 0.95);
  } else if (overflow.limit) {
    newBudget = Math.floor(overflow.limit * 0.9);
  } else {
    newBudget = Math.floor(budget / 2);
  }
  return Math.max(1, Math.min(newBudget, budget - 1));
}

// The first optional field of a request body that an error message blames, or null
export function findRejectedField(body, message, param, requiredFields) {
  const keys = Object.keys(body).filter((key) => !requiredFields.includes(key));
  const named = param ? param.split(".")[0].split("[")[0] : null;
  if (named && keys.includes(named)) return named;
  let best = null;
  let bestPosition = Infinity;
  for (const key of keys) {
    const m = new RegExp(`\\b${key.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}\\b`).exec(message || "");
    if (m && m.index < bestPosition) {
      best = key;
      bestPosition = m.index;
    }
  }
  return best;
}
