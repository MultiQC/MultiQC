import os
import sys
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union
import openai  # type: ignore
from openai.types.chat.chat_completion_message import ChatCompletionMessage
from pydantic.main import BaseModel  # type: ignore

from multiqc import config, report
from multiqc.plots.plotly.plot import Plot
from multiqc.types import Anchor, Section

from dotenv import load_dotenv  # type: ignore

if TYPE_CHECKING:
    import anthropic  # type: ignore

load_dotenv()


class SectionSummary(BaseModel):
    name: str
    summary: str


class SummaryResponse(BaseModel):
    section_summaries: List[SectionSummary]


PER_SECTION_SUMMARY_PROMPT = """\
You are an expert in bioinformatics, sequencing technologies, and genomics data analysis. You are given MultiQC report sections with the quality control data from a bioinformatics workflow.

Your task is to analyse the data and give a short and concise summary of each section. Be concise. Do not waste words. Only point to key issues.

Do no add any headers.

Use HTML to format lists, paragraphs, and add color and style to text, but not for anything else. 

Report sections are presented below.
"""
OVERALL_SUMMARY_PROMPT = """\
Now give an overall report summary. Make it very concise.

Do no add any headers.

Use HTML to format lists, paragraphs, and add color and style to text, but not for anything else.
"""


class LLMClient:
    def __init__(self, model: str):
        self.model = model
        self.history: List = [{"role": "system", "content": PER_SECTION_SUMMARY_PROMPT}]

    def generate_summaries(self, content: str) -> Tuple[Optional[SummaryResponse], Optional[str]]:
        raise NotImplementedError("Not implemented")


class OpenAIClient(LLMClient):
    def __init__(self, model: str, token: str):
        super().__init__(model)
        self.client = openai.OpenAI(api_key=token)

    def generate_summaries(self, content: str) -> Tuple[Optional[SummaryResponse], Optional[str]]:
        # First ask for a summary for each section
        self.history.append({"role": "user", "content": content})
        per_section_completion = self.client.beta.chat.completions.parse(
            model=self.model,
            messages=self.history,
            temperature=0.0,
            response_format=SummaryResponse,
        )
        per_section_summaries = per_section_completion.choices[0].message.parsed

        # Now asking for an overall summary
        self.history.append({"role": "assistant", "content": per_section_completion.choices[0].message.content})
        self.history.append({"role": "user", "content": OVERALL_SUMMARY_PROMPT})
        response = self.client.chat.completions.create(
            model=self.model,
            messages=self.history,
            temperature=0.0,
        )
        overall_summary_message: ChatCompletionMessage = response.choices[0].message
        self.history.append(overall_summary_message)
        return per_section_summaries, overall_summary_message.content


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

    def generate_summaries(self, content: str) -> Tuple[Optional[SummaryResponse], Optional[str]]:
        self.history += {"role": "user", "content": content}
        response = self.client.messages.create(
            model=self.model,
            messages=self.history,
            temperature=0.0,
        )
        ai_message = response.content[0]
        self.history.append(ai_message)
        self.history.append({"role": "user", "content": OVERALL_SUMMARY_PROMPT})
        response = self.client.messages.create(
            model=self.model,
            messages=self.history,
            temperature=0.0,
        )
        return None, ai_message.text


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


def add_ai_summary_to_report():
    if not (llm := get_llm_client()):
        return None

    content_for_llm = ""
    if report.general_stats_plot:
        content_for_llm += f"""
**id**: general_stats
**Title** MultiQC General Statistics
**Description**: Overview of key QC metrics for each sample.

{report.general_stats_plot.data_for_ai_prompt()}\
"""

    for section in report.get_all_sections():
        if section.plot_anchor and section.plot_anchor in report.plot_by_id:
            plot = report.plot_by_id[section.plot_anchor]
            if plot_prompt := plot.data_for_ai_prompt():
                content_for_llm += f"""
----------------------

**id**: {section.anchor}
**Tool**: {section.module} ({section.module_info})
**Title** {section.name}{f" {section.description}" if section.description else ""}
{f"**Extra plot description**: {section.helptext}\n" if section.helptext else ""}
{plot_prompt}\
"""

    if content_for_llm == "":
        return None

    # print(prompt)
    # sys.exit()

    tooltip = f"This block is AI-generated. Take with a grain of salt. Model: {llm.model}"

    per_section_summaries, overall_summary = llm.generate_summaries(content_for_llm)
    if overall_summary:
        report.ai_summary = f"<p style='color: gray' title='{tooltip}'>AI summary ✨</p>{overall_summary}"

    if per_section_summaries:
        summary_by_section_dict = {
            section_summary.name: section_summary.summary for section_summary in per_section_summaries.section_summaries
        }

        if genstats_summary := summary_by_section_dict.get("general_stats"):
            report.general_stats_ai_summary = (
                f"<p style='color: gray' title='{tooltip}'>AI summary ✨</p>{genstats_summary}"
            )

        for section in report.get_all_sections():
            ai_summary = summary_by_section_dict.get(section.id)
            if ai_summary:
                section.ai_summary = f"<p style='color: gray' title='{tooltip}'>AI summary ✨</p>{ai_summary}"
