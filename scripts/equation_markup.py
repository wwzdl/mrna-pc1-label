#!/usr/bin/env python3
"""Shared helpers for editable manuscript equations."""

from __future__ import annotations

import base64
from io import BytesIO
import re

from bs4 import BeautifulSoup, NavigableString, Tag
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls
from latex2mathml.converter import convert as latex_to_mathml
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import math_to_image
from mathml2omml import convert as mathml_to_omml


INLINE_MATH_PATTERN = re.compile(r"(?<!\\)\$(?!\$)(.+?)(?<!\\)\$")


def promote_inline_math(soup: BeautifulSoup) -> None:
    """Replace $...$ text outside code blocks with tagged inline-math nodes."""
    candidates = list(soup.find_all(string=lambda value: value and "$" in value))
    for node in candidates:
        if any(parent.name in {"code", "pre"} for parent in node.parents):
            continue
        text = str(node)
        matches = list(INLINE_MATH_PATTERN.finditer(text))
        if not matches:
            continue
        cursor = 0
        for match in matches:
            if match.start() > cursor:
                node.insert_before(NavigableString(text[cursor : match.start()]))
            span = soup.new_tag("span")
            span["class"] = "inline-math"
            span["data-latex"] = match.group(1).strip()
            node.insert_before(span)
            cursor = match.end()
        if cursor < len(text):
            node.insert_before(NavigableString(text[cursor:]))
        node.extract()


def equation_payload(tag: Tag) -> tuple[str, str] | None:
    """Return (LaTeX, number) for a fenced ```equation block."""
    if tag.name != "pre":
        return None
    code = tag.find("code", recursive=False)
    if code is None or "language-equation" not in code.get("class", []):
        return None
    lines = code.get_text().strip().splitlines()
    if not lines or not lines[0].lower().startswith("number:"):
        raise ValueError("equation blocks must begin with 'number: <label>'")
    number = lines[0].split(":", 1)[1].strip()
    latex = " ".join(line.strip() for line in lines[1:] if line.strip())
    if not number or not latex:
        raise ValueError("equation blocks require both a number and a LaTeX expression")
    return latex, number


def mathml(latex: str, display: str = "inline") -> str:
    """Convert LaTeX to MathML for HTML/PDF output."""
    return latex_to_mathml(latex, display=display)


def omml_element(latex: str):
    """Convert LaTeX to an editable Office Math (OMML) element."""
    omml = mathml_to_omml(mathml(latex, display="block"))
    if "xmlns:m=" not in omml:
        omml = omml.replace("<m:oMath>", f"<m:oMath {nsdecls('m')}>", 1)
    return parse_xml(omml)


def math_svg_data_uri(latex: str, font_size: float = 11) -> str:
    """Render LaTeX math as a vector SVG data URI for paginated PDF output."""
    buffer = BytesIO()
    math_to_image(
        f"${latex}$",
        buffer,
        format="svg",
        dpi=300,
        prop=FontProperties(size=font_size),
    )
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


def replace_math_for_html(soup: BeautifulSoup) -> None:
    """Replace equation fences and inline markers with MathML HTML nodes."""
    promote_inline_math(soup)

    for span in list(soup.select("span.inline-math")):
        latex = span.get("data-latex", "")
        image = soup.new_tag("img")
        image["class"] = "inline-equation-svg"
        image["alt"] = latex
        image["src"] = math_svg_data_uri(latex, font_size=11)
        span.replace_with(image)

    for pre in list(soup.find_all("pre")):
        payload = equation_payload(pre)
        if payload is None:
            continue
        latex, number = payload
        wrapper = soup.new_tag("div")
        wrapper["class"] = "display-equation"
        formula = soup.new_tag("div")
        formula["class"] = "equation-formula"
        image = soup.new_tag("img")
        image["class"] = "display-equation-svg"
        image["alt"] = latex
        image["src"] = math_svg_data_uri(latex, font_size=12)
        formula.append(image)
        label = soup.new_tag("div")
        label["class"] = "equation-number"
        label.string = f"({number})"
        wrapper.extend([formula, label])
        pre.replace_with(wrapper)
