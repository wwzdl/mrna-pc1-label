#!/usr/bin/env python3
"""Render manuscript Markdown to a paginated, submission-ready PDF."""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

import markdown
from weasyprint import HTML


CSS = r"""
@page {
  size: A4;
  margin: 22mm 19mm 22mm 19mm;
  @bottom-center {
    content: counter(page);
    font-family: "Liberation Serif", "Noto Serif CJK SC", serif;
    font-size: 9pt;
    color: #444;
  }
}

html {
  font-family: "Liberation Serif", "Noto Serif CJK SC", serif;
  font-size: 10.5pt;
  line-height: 1.42;
  color: #111;
}

body { margin: 0; }
h1, h2, h3, h4 {
  font-family: "Liberation Sans", "Noto Sans CJK SC", sans-serif;
  color: #111;
  page-break-after: avoid;
  break-after: avoid-page;
}
h1 { font-size: 16pt; line-height: 1.22; margin: 0 0 13pt; }
h2 { font-size: 13pt; margin: 15pt 0 7pt; }
h3 { font-size: 11.5pt; margin: 12pt 0 5pt; }
h4 { font-size: 10.5pt; margin: 10pt 0 4pt; }
p { margin: 0 0 7pt; orphans: 3; widows: 3; }
a { color: #164e8a; text-decoration: none; }
code {
  font-family: "Liberation Mono", "DejaVu Sans Mono", monospace;
  font-size: 0.88em;
  overflow-wrap: anywhere;
}
pre {
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  font-size: 8pt;
  background: #f5f5f5;
  border: 0.4pt solid #c8c8c8;
  padding: 6pt;
}
figure {
  margin: 11pt auto 12pt;
  break-inside: avoid-page;
  page-break-inside: avoid;
}
img {
  display: block;
  max-width: 100%;
  max-height: 205mm;
  width: auto;
  height: auto;
  margin: 0 auto 6pt;
}
figcaption { font-size: 9pt; line-height: 1.34; }
table {
  border-collapse: collapse;
  width: 100%;
  table-layout: fixed;
  margin: 6pt 0 11pt;
  font-size: 7.2pt;
  line-height: 1.22;
  break-inside: auto;
}
thead { display: table-header-group; }
tr { break-inside: avoid; page-break-inside: avoid; }
th, td {
  border-top: 0.45pt solid #777;
  border-bottom: 0.45pt solid #aaa;
  padding: 3pt 3.5pt;
  vertical-align: top;
  overflow-wrap: anywhere;
  word-break: normal;
}
th { font-weight: 700; background: #f1f3f5; }
blockquote {
  margin: 7pt 0 7pt 12pt;
  padding-left: 8pt;
  border-left: 2pt solid #888;
}
ul, ol { margin: 4pt 0 8pt 18pt; padding: 0; }
li { margin-bottom: 2pt; }
"""


FIGURE_PATTERN = re.compile(
    r"<p>\s*(<img\b[^>]*>)\s*</p>\s*"
    r"<p>\s*(<strong>(?:(?:Supplementary\s+)?Fig\.|补充图|图)\s*[^<]*</strong>.*?)</p>",
    flags=re.IGNORECASE | re.DOTALL,
)


def wrap_figures(fragment: str) -> str:
    return FIGURE_PATTERN.sub(r"<figure>\1<figcaption>\2</figcaption></figure>", fragment)


def render(input_path: Path, output_path: Path) -> None:
    source = input_path.read_text(encoding="utf-8")
    fragment = markdown.markdown(
        source,
        extensions=["tables", "fenced_code", "sane_lists"],
        output_format="html5",
    )
    fragment = wrap_figures(fragment)
    title_match = re.search(r"^#\s+(.+)$", source, flags=re.MULTILINE)
    title = title_match.group(1).strip() if title_match else input_path.stem
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>{CSS}</style>
</head>
<body>{fragment}</body>
</html>
"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    HTML(string=document, base_url=str(input_path.parent.resolve())).write_pdf(output_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input Markdown path")
    parser.add_argument("output", type=Path, help="Output PDF path")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    render(args.input.resolve(), args.output.resolve())
    print(f"Rendered {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
