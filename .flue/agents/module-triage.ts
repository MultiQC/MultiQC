import { readFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";
import type { FlueContext } from "@flue/runtime";
import { local } from "@flue/runtime/node";
import * as v from "valibot";

export const triggers = {};

interface ModuleTriagePayload {
  number?: number;
  title?: string;
  body?: string;
  author?: string;
  labels?: string[];
}

const here = dirname(fileURLToPath(import.meta.url));
const projectRoot = resolve(here, "..");

async function loadMarkdown(relativePath: string): Promise<string> {
  return readFile(resolve(projectRoot, relativePath), "utf8");
}

export default async function ({ init, payload }: FlueContext) {
  const issue = (payload ?? {}) as ModuleTriagePayload;
  const [repoContext, skillSpec] = await Promise.all([
    loadMarkdown("context.md"),
    loadMarkdown("skills/module-triage.md"),
  ]);

  const harness = await init({
    sandbox: local(),
    model: "anthropic/claude-sonnet-4-6",
  });
  const session = await harness.session();

  const { data } = await session.prompt(
    [
      "# Repository context",
      repoContext,
      "",
      "# Task",
      skillSpec,
      "",
      "# Module request",
      `#${issue.number ?? "?"} by ${issue.author ?? "unknown"}`,
      `Title: ${issue.title ?? "(no title)"}`,
      `Labels: ${(issue.labels ?? []).join(", ") || "(none)"}`,
      "",
      "Body:",
      issue.body ?? "(no body)",
    ].join("\n"),
    {
      result: v.object({
        tool: v.object({
          name: v.string(),
          purpose: v.string(),
          output_formats: v.array(v.string()),
        }),
        existing_coverage: v.string(),
        adoption_signal: v.picklist(["weak", "moderate", "strong"]),
        implementation_complexity: v.picklist(["low", "medium", "high"]),
        request_completeness: v.picklist(["incomplete", "partial", "complete"]),
        suggested_priority_label: v.picklist(["module: prio-low", "module: prio-medium", "module: prio-high"]),
        summary: v.string(),
        next_step: v.string(),
      }),
    },
  );

  return data;
}
