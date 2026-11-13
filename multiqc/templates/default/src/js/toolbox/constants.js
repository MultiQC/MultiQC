////////////////////////////////////////////////
// MultiQC Report Toolbox Constants
////////////////////////////////////////////////

// Make constants available globally
window.mqc_colours = [
  "#e41a1c",
  "#377eb8",
  "#4daf4a",
  "#984ea3",
  "#ff7f00",
  "#a9a904",
  "#a65628",
  "#f781bf",
  "#999999",
];
window.zip_threshold = 8;

// Add these constants at the top of the file
window.AI_PROVIDERS = {
  seqera: {
    name: "Seqera AI",
    apiKeysUrl: "https://cloud.seqera.io/tokens",
  },
  anthropic: {
    name: "Anthropic",
    // Defaults come from the MultiQC version that built the report, see DEFAULT_MODELS in ai_models.py
    get defaultModel() {
      return window.aiModelConfig?.defaultModels?.anthropic || "claude-opus-5-5";
    },
    suggestedModels: ["claude-opus-5-5", "claude-sonnet-5-5", "claude-haiku-5-5"],
    apiKeysUrl: "https://console.anthropic.com/settings/keys",
    modelsUrl: "https://platform.claude.com/docs/en/about-claude/models/overview",
  },
  openai: {
    name: "OpenAI",
    get defaultModel() {
      return window.aiModelConfig?.defaultModels?.openai || "gpt-6.1-sol";
    },
    suggestedModels: ["gpt-6.1-sol", "gpt-6.1-terra", "gpt-6-luna"],
    apiKeysUrl: "https://platform.openai.com/api-keys",
    modelsUrl: "https://developers.openai.com/api/docs/models",
  },
  aws_bedrock: {
    name: "AWS Bedrock",
    modelsUrl: "https://docs.aws.amazon.com/bedrock/latest/userguide/models-supported.html",
  },
  custom: {
    name: "Custom",
    defaultModel: "",
  },
  clipboard: {
    name: "Copy prompts",
  },
  none: {
    name: "Remove AI buttons",
  },
};

window.AI_PROVIDER_GROUPS = {
  "In-report summaries": ["seqera", "anthropic", "openai", "custom"],
  Alternatives: ["clipboard", "none"],
};

window.AUTO_SAVE_PREFIX = "autosave_";
