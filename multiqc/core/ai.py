import os
from typing import List, Optional
import openai
from openai.types.chat.chat_completion_message import ChatCompletionMessage

from multiqc import config
from multiqc.plots.table_object import DataTable
from multiqc.types import Section


class LLMClient:
    def __init__(self, model: str):
        self.model = model

    def chat(self, data_items: List[str]) -> str:
        raise NotImplementedError("Not implemented")


TABLE_SUMMARY_INSTRUCTION = """\
You are an expert in bioinformatics, sequencing technologies, and genomics data analysis. \
You are given quality control data from a bioinformatics workflow.

Your task is to eyeball the data and give a very short and helpful summary of the results \
to the customer. The summary should only warn about what a human would otherwise miss or leave unnoticed. \
Be specific: it's okay to simply say that all data is good, if no items are found to stand out as problematic.

Please use HTML formatting to make the summary visually attractive. Use colors to bring user's attention \
to really important pieces, but do not use a very bright color palette. Make highlighting very subtle.

Do no add a header or any other formatting apart from simple HTML.

The data is presented below."""


class OpenAIClient(LLMClient):
    def __init__(self, model: str, token: str):
        super().__init__(model)
        self.client = openai.OpenAI(api_key=token)
        self.history: List = []

    def chat(self, data_items: List[str]) -> Optional[str]:
        self.history += [{"role": "user", "content": data_item} for data_item in data_items]
        response = self.client.chat.completions.create(
            model=self.model,
            messages=self.history,
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

    def chat(self, data_items: List[str]) -> Optional[str]:
        self.history += [{"role": "user", "content": data_item} for data_item in data_items]

        response = self.client.messages.create(model=self.model, messages=self.history)
        ai_message = response.content[0]
        self.history.append(ai_message)
        return ai_message.text


def _get_llm_client() -> Optional[LLMClient]:
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


def table_llm_summary(
    dt: DataTable,
    table_title: str,
    report_section: Optional[Section] = None,
) -> Optional[str]:
    if not config.ai_summary:
        return None

    llm = _get_llm_client()
    if not llm:
        return None

    prompt = f"""
{TABLE_SUMMARY_INSTRUCTION}

#Title: {table_title}
"""
    if report_section:
        prompt += f"""\
#Description:
{report_section.description}
"""
    prompt += """\
#Data:

"""

    for table_section in dt.sections:
        for _, rows in table_section.rows_by_sgroup.items():
            row = rows[0]  # take only the first row in a group
            prompt += f"{dt.pconfig.col1_header}: {row.sample}\n{row.raw_data}\n\n"

    summary = llm.chat([prompt])
    return summary
