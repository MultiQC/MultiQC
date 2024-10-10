import os
from typing import List, Optional
import openai  # type: ignore
from openai.types.chat.chat_completion_message import ChatCompletionMessage  # type: ignore

from multiqc import config, report
from multiqc.plots.plotly.plot import Plot
from multiqc.plots.table_object import DataTable
from multiqc.types import Anchor, Section

from dotenv import load_dotenv  # type: ignore

load_dotenv()


class LLMClient:
    def __init__(self, model: str):
        self.model = model

    def chat(self, prompt: str) -> Optional[str]:
        raise NotImplementedError("Not implemented")


SYSTEM_PROMPT = """\
You are an expert in bioinformatics, sequencing technologies, and genomics data analysis. You are given a MultiQC summary of the quality control data from a bioinformatics workflow.

Your task is to eyeball the data and give a very short and helpful summary of the results to the customer. The summary should only warn about what a human would otherwise miss or leave unnoticed.

Be concise: only mention samples that stand out as problematic.

Point to the title of the section to give the reader context.

Use HTML to format lists, paragraphs, and style to text, but not for anything else.

Do no add any headers.

The data is presented below.
"""


class OpenAIClient(LLMClient):
    def __init__(self, model: str, token: str):
        super().__init__(model)
        self.client = openai.OpenAI(api_key=token)
        self.history: List = [{"role": "system", "content": SYSTEM_PROMPT}]

    def chat(self, prompt: str) -> Optional[str]:
        self.history.append({"role": "user", "content": prompt})
        response = self.client.chat.completions.create(
            model=self.model,
            messages=self.history,
            temperature=0.0,
        )
        ai_message: ChatCompletionMessage = response.choices[0].message
        self.history.append(ai_message)
        return ai_message.content


class AnthropicClient(LLMClient):
    def __init__(self, model: str, token: str):
        try:
            import anthropic  # type: ignore
        except ImportError:
            raise ImportError(
                "anthropic package is not installed, make sure to install MultiQC with `pip install multiqc[anthropic]`"
            )

        super().__init__(model)
        self.client = anthropic.Anthropic(api_key=token)
        self.history: List = []

    def chat(self, prompt: str) -> Optional[str]:
        self.history += {"role": "user", "content": prompt}

        response = self.client.messages.create(
            model=self.model,
            messages=self.history,
            temperature=0.0,
        )
        ai_message = response.content[0]
        self.history.append(ai_message)
        return ai_message.text


def get_llm_client() -> Optional[LLMClient]:
    if not config.ai_summary:
        return None

    openai_token = os.environ.get("OPENAI_API_KEY")
    if openai_token:
        return OpenAIClient(
            model=os.environ.get("OPENAI_MODEL", "gpt-4o"),
            token=openai_token,
        )
    anthropic_token = os.environ.get("ANTHROPIC_API_KEY")
    if anthropic_token:
        return AnthropicClient(
            model=os.environ.get("ANTHROPIC_MODEL", "claude-3-5-sonnet-20240620"),
            token=anthropic_token,
        )

    return None


def generate_ai_summary() -> Optional[str]:
    if not (llm := get_llm_client()):
        return None

    prompt = ""
    if report.general_stats_plot:
        prompt += f"""
**Title** MultiQC General Statistics
**Description**: Overview of key QC metrics for each sample.
**Data** {report.general_stats_plot.data_for_ai_prompt()}
"""

    for section in report.get_all_sections():
        if section.plot_anchor and section.plot_anchor in report.plot_by_id:
            plot = report.plot_by_id[section.plot_anchor]
            if plot_prompt := plot.data_for_ai_prompt():
                prompt += f"""
----------------------

**Tool**: {section.module}
**Title** {plot.pconfig.title}
**Description**: {section.description}
{f"**Extra plot description**: {section.helptext}" if section.helptext else ""}
**Data**
{plot_prompt}
                """

    if not prompt:
        return None

    summary = llm.chat(prompt)
    if summary:
        # insert emoji inside the first paragraph of the summary
        summary = summary.replace("<p>", "<p class='first-line'>✨ ", 1)
    return summary
