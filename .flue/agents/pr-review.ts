import { readFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";
import type { FlueContext } from "@flue/runtime";
import { local } from "@flue/runtime/node";
import * as v from "valibot";

export const triggers = {};

interface ReviewPayload {
  number?: number;
  title?: string;
  body?: string;
  author?: string;
  diff?: string;
  changed_files?: string[];
}

const here = dirname(fileURLToPath(import.meta.url));
const projectRoot = resolve(here, "..");

async function loadMarkdown(relativePath: string): Promise<string> {
  return readFile(resolve(projectRoot, relativePath), "utf8");
}

export default async function ({ init, payload }: FlueContext) {
  const pr = (payload ?? {}) as ReviewPayload;
  const [repoContext, skillSpec] = await Promise.all([loadMarkdown("context.md"), loadMarkdown("skills/pr-review.md")]);

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
      "# Pull request",
      `#${pr.number ?? "?"} by ${pr.author ?? "unknown"}`,
      `Title: ${pr.title ?? "(no title)"}`,
      "",
      "Body:",
      pr.body ?? "(no body)",
      "",
      "Changed files:",
      (pr.changed_files ?? []).map((f) => `- ${f}`).join("\n") || "(none)",
      "",
      "Diff:",
      "```diff",
      pr.diff ?? "(no diff supplied)",
      "```",
    ].join("\n"),
    {
      result: v.object({
        summary: v.string(),
        findings: v.array(
          v.object({
            category: v.picklist(["bug", "style", "perf", "security", "test-coverage", "docs"]),
            severity: v.picklist(["low", "medium", "high"]),
            file: v.optional(v.string()),
            note: v.string(),
          }),
        ),
        suggestions: v.array(v.string()),
        verdict: v.picklist(["approve", "comment", "request-changes"]),
      }),
    },
  );

  return data;
}
