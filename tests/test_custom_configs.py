"""
Tests for custom configuration options like custom_favicon and custom_logo.
"""

from pathlib import Path

import pytest

import multiqc
from multiqc.core.update_config import ClConfig


@pytest.fixture
def sample_data_file(tmp_path: Path) -> Path:
    """Create a simple custom content data file for generating reports."""
    data_file = tmp_path / "data_mqc.txt"
    data_file.write_text("sample1\t100\nsample2\t200\n")
    return data_file


class TestCustomFavicon:
    """Tests for custom_favicon configuration."""

    def test_custom_favicon_png_in_html(self, data_dir, sample_data_file: Path, tmp_path: Path) -> None:
        """Verify custom PNG favicon appears in generated HTML with correct MIME type."""
        favicon_path = data_dir / "custom_configs" / "favicon.png"
        config_file = tmp_path / "config.yaml"
        config_file.write_text(f'custom_favicon: "{favicon_path}"')

        multiqc.run(
            sample_data_file,
            cfg=ClConfig(
                force=True,
                output_dir=str(tmp_path),
                config_files=[config_file],
            ),
        )

        report_html = (tmp_path / "multiqc_report.html").read_text()
        expected_favicon_html = """\
<link
  rel="icon"
  type="image/png"
  href="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=="
/>"""
        assert expected_favicon_html in report_html

    def test_custom_favicon_svg_in_html(self, data_dir, sample_data_file: Path, tmp_path: Path) -> None:
        """Verify custom SVG favicon appears in generated HTML with correct MIME type."""
        favicon_path = data_dir / "custom_configs" / "favicon.svg"
        config_file = tmp_path / "config.yaml"
        config_file.write_text(f'custom_favicon: "{favicon_path}"')

        multiqc.run(
            sample_data_file,
            cfg=ClConfig(
                force=True,
                output_dir=str(tmp_path),
                config_files=[config_file],
            ),
        )

        report_html = (tmp_path / "multiqc_report.html").read_text()
        expected_favicon_html = """\
<link
  rel="icon"
  type="image/svg+xml"
  href="data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIzMiIgaGVpZ2h0PSIzMiIgdmlld0JveD0iMCAwIDMyIDMyIj4KICA8cmVjdCB3aWR0aD0iMzIiIGhlaWdodD0iMzIiIGZpbGw9IiM0Q0FGNTAiLz4KPC9zdmc+Cg=="
/>"""
        assert expected_favicon_html in report_html

    def test_default_favicon_when_not_configured(self, sample_data_file: Path, tmp_path: Path) -> None:
        """Verify default favicon is used when custom_favicon is not set."""
        multiqc.run(
            sample_data_file,
            cfg=ClConfig(
                force=True,
                output_dir=str(tmp_path),
            ),
        )

        report_html = (tmp_path / "multiqc_report.html").read_text()
        expected_favicon_html = """\
<link
  rel="icon"
  type="image/svg+xml"
  href="data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMjUxIiBoZWlnaHQ9IjI1MSIgdmlld0JveD0iMCAwIDI1MSAyNTEiIGZpbGw9Im5vbmUiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+CjxwYXRoIGQ9Ik00Ni42NjAxIDEyMC41NDVDNDkuMjgwMSA4MS4wOTcyIDgwLjk3MDEgNDkuNDg5NCAxMjAuNDYgNDcuMDA5NVYwLjk0MjYyN0M1NS41MDAxIDMuNTAyNDUgMy4yODAwOCA1NS42MTkgMC41ODAwNzggMTIwLjU0NUg0Ni42NjAxWiIgZmlsbD0iI0YxODA0NiIvPgo8cGF0aCBkPSJNMTIwLjE5IDIwNC44NDlDODAuNzQwMSAyMDIuMjI5IDQ5LjEzMDEgMTcwLjU0MiA0Ni42NTAxIDEzMS4wNTRIMC41ODAwNzhDMy4xNDAwOCAxOTYuMDEgNTUuMjYwMSAyNDguMjI2IDEyMC4xOSAyNTAuOTI2VjIwNC44NDlaIiBmaWxsPSIjRjE4MDQ2Ii8+CjxwYXRoIGQ9Ik0xMzAuOTY5IDQ3LjAxOTVDMTcwLjQxOSA0OS42Mzk0IDIwMi4wMjkgODEuMzI3MiAyMDQuNTA5IDEyMC44MTVIMjUwLjU3OUMyNDguMDE5IDU1Ljg1ODkgMTk1Ljg5OSAzLjY0MjQ1IDEzMC45NjkgMC45NDI2MjdWNDcuMDE5NVoiIGZpbGw9IiNGMTgwNDYiLz4KPHBhdGggZD0iTTI1MC41NzkgMjA0Ljg0OUMyMTEuMTI5IDIwMi4yMjkgMTc5LjUxOSAxNzAuNTQyIDE3Ny4wMzkgMTMxLjA1NEgxMzAuOTY5QzEzMy41MjkgMTk2LjAxIDE4NS42NDkgMjQ4LjIyNiAyNTAuNTc5IDI1MC45MjZWMjA0Ljg0OVoiIGZpbGw9IiNGMTgwNDYiLz4KPC9zdmc+Cg=="
/>"""
        assert expected_favicon_html in report_html


class TestCustomLogo:
    """Tests for custom_logo configuration with MIME type detection."""

    def test_custom_logo_png_in_html(self, data_dir, sample_data_file: Path, tmp_path: Path) -> None:
        """Verify custom PNG logo appears with correct MIME type."""
        logo_path = data_dir / "custom_configs" / "logo.png"
        config_file = tmp_path / "config.yaml"
        config_file.write_text(f'''\
custom_logo: "{logo_path}"
custom_logo_url: "https://customlogo.com"
''')

        multiqc.run(
            sample_data_file,
            cfg=ClConfig(
                force=True,
                output_dir=str(tmp_path),
                config_files=[config_file],
            ),
        )

        report_html = (tmp_path / "multiqc_report.html").read_text()
        expected_logo_html = """\
      <a href="https://customlogo.com" target="_blank">
      <img
        src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=="
        title=""
        class="custom_logo custom_logo_light"
        
      />
      
      </a>"""

        assert expected_logo_html in report_html

    def test_custom_logo_svg_in_html(self, data_dir, sample_data_file: Path, tmp_path: Path) -> None:
        """Verify custom SVG logo appears with correct MIME type."""
        logo_path = data_dir / "custom_configs" / "logo.svg"
        config_file = tmp_path / "config.yaml"
        config_file.write_text(f'''\
custom_logo: "{logo_path}"
custom_logo_url: "https://customlogo.com"
''')

        multiqc.run(
            sample_data_file,
            cfg=ClConfig(
                force=True,
                output_dir=str(tmp_path),
                config_files=[config_file],
            ),
        )

        report_html = (tmp_path / "multiqc_report.html").read_text()
        expected_logo_html = """\
      <a href="https://customlogo.com" target="_blank">
      <img
        src="data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxMDAiIGhlaWdodD0iMzAiIHZpZXdCb3g9IjAgMCAxMDAgMzAiPgogIDxyZWN0IHdpZHRoPSIxMDAiIGhlaWdodD0iMzAiIGZpbGw9IiMyMTk2RjMiLz4KICA8dGV4dCB4PSI1MCIgeT0iMjAiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZpbGw9IndoaXRlIiBmb250LXNpemU9IjEyIj5Mb2dvPC90ZXh0Pgo8L3N2Zz4K"
        title=""
        class="custom_logo custom_logo_light"
        
      />
      
      </a>"""

        assert expected_logo_html in report_html
