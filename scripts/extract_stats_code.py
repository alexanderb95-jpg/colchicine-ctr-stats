#!/usr/bin/env python3
"""Extract full statistical R code from CTR R Markdown into one .R file."""

from __future__ import annotations

import html
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RMD = ROOT / "CTR_Submission_Results_Protocol2101175.Rmd"
HELPERS = ROOT / "protocol2101175_redcap_helpers.R"
OUT = ROOT / "colchicine_ctr_statistical_code.R"
INDEX = ROOT / "index.html"
STAT_HTML = ROOT / "statistical-analysis.html"


def extract_full_code(rmd_text: str, helpers_text: str) -> str:
    pattern = re.compile(r"```\{r\s*([^,\}]*)([^\}]*)\}\s*\n(.*?)```", re.S)
    chunks = pattern.findall(rmd_text)
    parts = [
        "# STUDY-21-001175 — full statistical analysis code (extracted from CTR R Markdown)",
        "# Source: CTR_Submission_Results_Protocol2101175.Rmd",
        "# Data: ProtocolNo2101175PRM_DATA_2026-05-18_2135.csv",
        "# Helpers: protocol2101175_redcap_helpers.R",
        "",
        "# =============================================================================",
        "# protocol2101175_redcap_helpers.R",
        "# =============================================================================",
        helpers_text.rstrip(),
        "",
    ]
    for name, opts, body in chunks:
        name = name.strip() or "unnamed"
        parts.extend(
            [
                "# =============================================================================",
                f"# R chunk: {name}  ({opts.strip()})",
                "# =============================================================================",
                body.rstrip(),
                "",
            ]
        )
    return "\n".join(parts).rstrip() + "\n"


def inject_code_block(html_text: str, code_esc: str, pre_open: str) -> str:
    start = html_text.find(pre_open)
    if start < 0:
        return html_text
    code_start = html_text.find("<code>", start)
    code_end = html_text.find("</code>", code_start)
    if code_start < 0 or code_end < 0:
        return html_text
    return html_text[: code_start + len("<code>")] + code_esc + html_text[code_end:]


def main() -> None:
    rmd_text = RMD.read_text(encoding="utf-8")
    helpers_text = HELPERS.read_text(encoding="utf-8")
    full_code = extract_full_code(rmd_text, helpers_text)
    code_esc = html.escape(full_code)
    OUT.write_text(full_code, encoding="utf-8")

    if STAT_HTML.exists():
        t = STAT_HTML.read_text(encoding="utf-8")
        STAT_HTML.write_text(inject_code_block(t, code_esc, "<pre><code>"), encoding="utf-8")

    if INDEX.exists():
        t = INDEX.read_text(encoding="utf-8")
        marker = '<div id="statistical-code-appendix"'
        end_marker = "<!-- end statistical-code-appendix -->"
        if marker in t and end_marker in t:
            t = inject_code_block(t, code_esc, '<pre style="background:#f6f8fa')
            INDEX.write_text(t, encoding="utf-8")

    print(f"Wrote {OUT} ({len(full_code.splitlines())} lines)")


if __name__ == "__main__":
    main()
