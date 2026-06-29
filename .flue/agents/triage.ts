import type { FlueContext } from "@flue/runtime";
import { local } from "@flue/runtime/node";
import * as v from "valibot";

export const triggers = {};

interface TriagePayload {
  number?: number;
  title?: string;
  body?: string;
  author?: string;
}

export default async function ({ init, payload }: FlueContext) {
  const issue = (payload ?? {}) as TriagePayload;

  const harness = await init({
    sandbox: local(),
    model: "anthropic/claude-sonnet-4-6",
  });
  const session = await harness.session();

  const { data } = await session.prompt(
    [
      "You are triaging a new GitHub issue for the MultiQC repository.",
      "MultiQC is a Python tool that aggregates bioinformatics QC reports.",
      "",
      `Issue #${issue.number ?? "?"} by ${issue.author ?? "unknown"}`,
      `Title: ${issue.title ?? "(no title)"}`,
      "",
      "Body:",
      issue.body ?? "(no body)",
      "",
      "Classify the issue and propose next steps.",
    ].join("\n"),
    {
      result: v.object({
        summary: v.string(),
        category: v.picklist(["bug", "feature-request", "module-request", "question", "docs", "other"]),
        suggested_labels: v.array(v.string()),
        priority: v.picklist(["low", "medium", "high"]),
        next_steps: v.string(),
      }),
    },
  );

  return data;
}
