import os
import sys
from textwrap import dedent
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union
import openai  # type: ignore
from openai.types.chat.chat_completion_message import ChatCompletionMessage
from openai.types.chat.chat_completion_message_param import ChatCompletionMessageParam
from pydantic.main import BaseModel  # type: ignore

from multiqc import config, report
from multiqc.plots.plotly.plot import Plot
from multiqc.types import Anchor, Section

from dotenv import load_dotenv  # type: ignore

if TYPE_CHECKING:
    import anthropic  # type: ignore

load_dotenv()


class SectionSummary(BaseModel):
    id: str
    summary: str


class SummaryResponse(BaseModel):
    section_summaries: List[SectionSummary]


OVERALL_SUMMARY_PROMPT = """\
You are an expert in bioinformatics, sequencing technologies, and genomics data analysis. You are given key findings from a MultiQC report from a bioinformatics workflow.

Your task is to analyse the data and give a short and concise summary of the section. Be concise. Do not waste words. Point to the key findings and results.

Do not add any headers.

Use HTML to format lists, paragraphs, and add color and style to text, but not for anything else. 

Data is presented below.
"""


class LLMClient:
    def __init__(self, model: str):
        self.model = model
        self.history: List = []

    def completion(
        self, messages: List[ChatCompletionMessageParam], **kwargs
    ) -> Union[ChatCompletionMessage, SummaryResponse]:
        raise NotImplementedError("Not implemented")


class OpenAIClient(LLMClient):
    def __init__(self, model: str, token: str):
        super().__init__(model)
        self.client = openai.OpenAI(api_key=token)

    # def generate_summaries(self, content: str) -> Tuple[Optional[SummaryResponse], Optional[str]]:
    #     self.history.append({"role": "system", "content": SECTION_SUMMARY_PROMPT})
    #     self.history.append({"role": "user", "content": content})
    #     per_section_completion = self.client.beta.chat.completions.parse(
    #         model=self.model,
    #         messages=self.history,
    #         temperature=0.0,
    #         response_format=SummaryResponse,
    #     )
    #     per_section_summaries = per_section_completion.choices[0].message.parsed

    #     # Now asking for an overall summary
    #     self.history.append({"role": "assistant", "content": per_section_completion.choices[0].message.content})
    #     self.history.append({"role": "user", "content": OVERALL_SUMMARY_PROMPT})
    #     response = self.client.chat.completions.create(
    #         model=self.model,
    #         messages=self.history,
    #         temperature=0.0,
    #     )
    #     overall_summary_message: ChatCompletionMessage = response.choices[0].message
    #     self.history.append(overall_summary_message)
    #     return per_section_summaries, overall_summary_message.content

    def completion(self, messages: List[ChatCompletionMessageParam], **kwargs) -> ChatCompletionMessage:
        response = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            temperature=0.0,
            **kwargs,
        )
        return response.choices[0].message


# class AnthropicClient(LLMClient):
#     def __init__(self, model: str, token: str):
#         try:
#             import anthropic  # type: ignore
#         except ImportError:
#             raise ImportError(
#                 "anthropic package is not installed, make sure to install MultiQC with `pip install multiqc[anthropic]`"
#             )

#         super().__init__(model)
#         self.client = anthropic.Anthropic(api_key=token)
#         self.history: List = []

#     def completion(self, content: str) -> Optional[str]:
#         self.history += {"role": "user", "content": content}
#         response = self.client.messages.create(
#             model=self.model,
#             messages=self.history,
#             temperature=0.0,
#         )
#         self.history.append(response.content[0])
#         return response.content[0].text
#         # self.history.append({"role": "user", "content": OVERALL_SUMMARY_PROMPT})
#         # response = self.client.messages.create(
#         #     model=self.model,
#         #     messages=self.history,
#         #     temperature=0.0,
#         # )
#         # return None, ai_message.text


def get_llm_client() -> Optional[OpenAIClient]:
    if not config.ai_summary:
        return None

    openai_token = os.environ.get("OPENAI_API_KEY")
    if openai_token:
        return OpenAIClient(
            model=config.ai_model,
            token=openai_token,
        )
    # anthropic_token = os.environ.get("ANTHROPIC_API_KEY")
    # if anthropic_token:
    #     return AnthropicClient(
    #         model=os.environ.get("ANTHROPIC_MODEL", "claude-3-5-sonnet-20240620"),
    #         token=anthropic_token,
    #     )

    return None


def _summarize_section(
    llm: LLMClient,
    content: str,
) -> Optional[str]:
    return llm.completion(
        messages=[
            {
                "role": "system",
                "content": dedent("""
                You are an expert in bioinformatics, sequencing technologies, and genomics data analysis. 
                You are given quality control data from a MultiQC report section generated in a bioinformatics workflow.
                Your task is to analyse this data and give a summary of this data. Use no more than one line per sample.
                Do not add any headers.
                Use HTML to format lists, paragraphs, and add color and style to text, but not for anything else. 
                Data is presented below.
                """).strip(),
            },
            {"role": "user", "content": content},
        ],
    )


TOOLTIP = f"This block is AI-generated. Take with a grain of salt. Model: {config.ai_model}"


def add_ai_summary_to_report():
    """
    Ask the LLM to summarize data in each report section with a table, a bar plot, and the general stats, all in parallel.
    On a secon iteration, pass all the data and all the summaries to the LLM and ask it to generate a final summary.
    Ask it generate a short version of the final summary.
    """

    openai_token = os.environ.get("OPENAI_API_KEY")
    if not openai_token:
        return

    llm = get_llm_client()
    if not llm:
        return

    content = ""
    if report.general_stats_plot:
        content += dedent(f"""
            **id**: general_stats
            **Title** MultiQC General Statistics (Overview of key QC metrics for each sample)

            {report.general_stats_plot.data_for_ai_prompt()}
            """)

    sections = report.get_all_sections()
    for section in sections:
        if section.plot_anchor and section.plot_anchor in report.plot_by_id:
            plot = report.plot_by_id[section.plot_anchor]
            if plot_content := plot.data_for_ai_prompt():
                content += dedent(f"""\
                    **id**: {section.anchor}
                    **Tool**: {section.module} ({section.module_info})
                    **Title** {section.name} {f"({section.description})" if section.description else ""}
                    {f"\n**Extra plot description**: {section.helptext}" if section.helptext else ""}

                    {plot_content}

                    ----------------------
                """)

    if not content:
        return

    per_section_response = llm.client.beta.chat.completions.parse(
        model=llm.model,
        messages=[
            {
                "role": "system",
                "content": dedent(f"""
                You are an expert in bioinformatics, sequencing technologies, and genomics data analysis. 

                You are given quality control data from a MultiQC report generated in a bioinformatics workflow.
                                  
                The data is split by section. ID of each section is given before the section contents in a line starting with **id**.

                Your task is to analyse this data and give a summary for each section, indexed by the section title.
                
                Use HTML to format lists, paragraphs, and add color and style to text, but not for anything else. 
                Do not add any headers.
                {content}
                """).strip(),
            }
        ],
        temperature=0.0,
        response_format=SummaryResponse,
    )

    if per_section_response.choices[0].message.parsed:
        summary_by_section_dict = {
            s.id: s.summary for s in per_section_response.choices[0].message.parsed.section_summaries
        }
        for section in report.get_all_sections():
            if summary := summary_by_section_dict.get(str(section.anchor)):
                section.ai_summary = f"<p style='color: gray' title='{TOOLTIP}'>AI summary ✨</p>{summary}"

        if gen_stats_summary := summary_by_section_dict.get("general_stats"):
            report.general_stats_ai_summary = (
                f"<p style='color: gray' title='{TOOLTIP}'>AI summary ✨</p>{gen_stats_summary}"
            )

    # Now asking for an overall summary
    messages = [
        {
            "role": "system",
            "content": dedent("""
            You are an expert in bioinformatics, sequencing technologies, and genomics data analysis. 
            You are given key findings from a MultiQC report generated by a bioinformatics workflow, per section.

            Your task is to analyse the data and give a short and concise overall summary. Do not list every section,
            but rather give a high-level summary. Only list the samples that have critical QC issues.
            Be very concise. Use colors to color code samples.

            Do not add any headers.

            Use HTML to format lists, paragraphs, and add color and style to text, but for nothing else.
            """).strip(),
        },
        {
            "role": "user",
            "content": per_section_response.choices[0].message.content,
        },
    ]
    completion = llm.client.chat.completions.create(
        model=llm.model,
        messages=messages,
        temperature=0.0,
    )
    report.ai_summary = (
        f"<p style='color: gray' title='{TOOLTIP}'>AI summary ✨</p>{completion.choices[0].message.content}"
    )

    # # per_section_summaries, overall_summary = llm.generate_summaries(content_for_llm)
    # if overall_summary:
    #     report.ai_summary = f"<p style='color: gray' title='{TOOLTIP}'>AI summary ✨</p>{overall_summary}"

    # if per_section_summaries:
    #     summary_by_section_dict = {
    #         section_summary.name: section_summary.summary for section_summary in per_section_summaries.section_summaries
    #     }

    #     if genstats_summary := summary_by_section_dict.get("general_stats"):
    #         report.general_stats_ai_summary = (
    #             f"<p style='color: gray' title='{tooltip}'>AI summary ✨</p>{genstats_summary}"
    #         )

    #     for section in report.get_all_sections():
    #         ai_summary = summary_by_section_dict.get(section.id)
    #         if ai_summary:
    #             section.ai_summary = f"<p style='color: gray' title='{tooltip}'>AI summary ✨</p>{ai_summary}"
