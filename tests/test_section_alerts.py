import pytest
from pydantic import ValidationError

from multiqc.base_module import BaseMultiqcModule
from multiqc.types import SectionAlert


def test_add_section_formats_alert_markdown_and_samples():
    module = BaseMultiqcModule(name="Test", anchor="test")

    module.add_section(
        name="Only alert",
        alerts={
            "message": "**2 samples** hidden from plot.",
            "level": "warning",
            "affected_samples": ["sample-a", "sample-b"],
        },
    )

    section = module.sections[0]
    assert section.print_section is True
    assert section.alerts[0].level == "warning"
    assert section.alerts[0].message == "<p><strong>2 samples</strong> hidden from plot.</p>"
    assert section.alerts[0].affected_samples == ["sample-a", "sample-b"]


def test_add_section_accepts_plain_string_alert():
    module = BaseMultiqcModule(name="Test", anchor="test")

    module.add_section(name="String alert", alerts="Plain **markdown** alert.")

    assert module.sections[0].alerts[0] == SectionAlert(
        message="<p>Plain <strong>markdown</strong> alert.</p>",
        level="info",
    )


def test_section_alert_rejects_unsafe_level():
    with pytest.raises(ValidationError):
        SectionAlert(message="Bad level", level='info" onclick="alert(1)')


def test_section_alert_rejects_non_bootstrap_level():
    with pytest.raises(ValidationError):
        SectionAlert(message="Bad level", level="important")


def test_add_section_skips_empty_alert_message_with_samples():
    module = BaseMultiqcModule(name="Test", anchor="test")

    module.add_section(
        name="Empty alert",
        alerts={"message": "", "affected_samples": ["sample-a"]},
    )

    assert module.sections[0].alerts == []
    assert module.sections[0].print_section is False
