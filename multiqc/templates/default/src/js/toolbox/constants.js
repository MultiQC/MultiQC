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
    defaultModel: "claude-3-5-sonnet-latest",
    suggestedModels: ["claude-3-5-sonnet-latest", "claude-3-5-haiku-latest"],
    apiKeysUrl: "https://console.anthropic.com/settings/keys",
    modelsUrl: "https://docs.anthropic.com/en/docs/intro-to-claude#model-options",
  },
  openai: {
    name: "OpenAI",
    defaultModel: "gpt-4o",
    suggestedModels: ["gpt-4o", "gpt-4o-mini"],
    apiKeysUrl: "https://platform.openai.com/api-keys",
    modelsUrl: "https://platform.openai.com/docs/models",
  },
  aws_bedrock: {
    name: "AWS Bedrock",
    modelsUrl: "https://docs.anthropic.com/en/docs/intro-to-claude#model-options",
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
