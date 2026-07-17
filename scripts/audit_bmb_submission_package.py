#!/usr/bin/env python3
"""Audit the local BMB submission package for common hard-format issues."""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from PIL import Image

from equation_markup import omml_element


ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / "manuscript" / "bmb_submission"

MAIN_EN = SUB / "bmb_main_manuscript_en.md"
SUPP_EN = SUB / "bmb_supplementary_material_en.md"
MAIN_CN = SUB / "bmb_main_manuscript_cn.md"
SUPP_CN = SUB / "bmb_supplementary_material_cn.md"
TITLE_EN = SUB / "bmb_title_page_en.md"
COVER_EN = SUB / "bmb_cover_letter_en.md"
DECL_EN = SUB / "bmb_statements_and_declarations_en.md"
CITATION = ROOT / "CITATION.cff"
BIB = SUB / "references" / "bmb_references.bib"
RIS = SUB / "references" / "bmb_references_for_endnote.ris"

MAIN_DOCX = SUB / "bmb_main_manuscript_en.docx"
SUPP_DOCX = SUB / "bmb_supplementary_material_en.docx"
SUPP_PDF = SUB / "bmb_supplementary_material_en.pdf"

REQUIRED_MD = [
    MAIN_EN,
    SUPP_EN,
    TITLE_EN,
    COVER_EN,
    DECL_EN,
]

REQUIRED_DOCX = [
    SUB / "bmb_main_manuscript_en.docx",
    SUB / "bmb_supplementary_material_en.docx",
    SUB / "bmb_title_page_en.docx",
    SUB / "bmb_cover_letter_en.docx",
    SUB / "bmb_statements_and_declarations_en.docx",
]

SKIP_PDF = os.environ.get("BMB_SKIP_PDF", "0") == "1"
REQUIRED_PDF = [] if SKIP_PDF else [SUPP_PDF]

REQUIRED_STATEMENT_HEADINGS = [
    "Funding",
    "Author Contributions",
    "Data Availability",
    "Code Availability",
    "Competing Interests",
    "Ethics Approval",
]

SI_METADATA_MARKERS = [
    "**Journal:**",
    "**Article title:**",
    "**Authors:**",
    "**Affiliation:**",
    "**Corresponding author:**",
]

@dataclass
class Check:
    level: str
    message: str


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_row(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    matches = [row for row in rows if all(row.get(key) == value for key, value in criteria.items())]
    if len(matches) != 1:
        raise ValueError(f"Expected one row for {criteria}, found {len(matches)}")
    return matches[0]


def plain_text(markdown: str) -> str:
    text = re.sub(r"```.*?```", " ", markdown, flags=re.S)
    text = re.sub(r"!\[[^\]]*\]\([^)]+\)", " ", text)
    text = re.sub(r"\[[^\]]+\]\([^)]+\)", " ", text)
    text = re.sub(r"[#*_`>|]", " ", text)
    return " ".join(text.split())


def extract_section(text: str, heading: str) -> str:
    pattern = re.compile(rf"(?ms)^##\s+{re.escape(heading)}\s*$\n(.*?)(?=^##\s+|\Z)")
    match = pattern.search(text)
    return match.group(1).strip() if match else ""


def split_main_references(text: str) -> tuple[str, list[str]]:
    match = re.search(r"(?ms)^##\s+\d+\.\s+References\s*$", text)
    if not match:
        return text, []
    body = text[: match.start()]
    ref_text = text[match.end() :]
    refs = [line.strip() for line in ref_text.splitlines() if line.strip()]
    return body, refs


def split_supplement_references(text: str) -> tuple[str, list[str]]:
    match = re.search(r"(?ms)^##\s+S16\.\s+References\s*$", text)
    if not match:
        return text, []
    body = text[: match.start()]
    ref_text = text[match.end() :]
    refs = [line.strip() for line in ref_text.splitlines() if line.strip()]
    return body, refs


def reference_family_year(line: str) -> tuple[str, str] | None:
    match = re.match(r"([A-Z][A-Za-z'\\-]+)\b.*?\((19\d{2}|20\d{2})\)", line)
    if not match:
        return None
    return match.group(1), match.group(2)


def citation_family_years(body: str) -> set[tuple[str, str]]:
    cited: set[tuple[str, str]] = set()
    for parenthetical in re.findall(r"\(([^()]*?\b(?:19|20)\d{2}[a-z]?[^()]*)\)", body):
        for part in parenthetical.split(";"):
            years = re.findall(r"\b((?:19|20)\d{2})[a-z]?\b", part)
            first = re.search(r"\b([A-Z][A-Za-z'\\-]+)\b", part)
            if not first:
                continue
            for year in years:
                cited.add((first.group(1), year))
    return cited


def count_words(text: str) -> int:
    return len(re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)?", plain_text(text)))


def count_markdown_tables(text: str) -> int:
    """Count Markdown table separator rows, one per table."""
    count = 0
    for line in text.splitlines():
        if not line.lstrip().startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if cells and all(re.fullmatch(r":?-{2,}:?", cell) for cell in cells):
            count += 1
    return count


def docx_xml(path: Path, member: str) -> str:
    with zipfile.ZipFile(path) as zf:
        try:
            return zf.read(member).decode("utf-8", errors="ignore")
        except KeyError:
            return ""


def docx_visible_text(path: Path) -> str:
    namespace = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    root = ET.fromstring(docx_xml(path, "word/document.xml"))
    return "".join(node.text or "" for node in root.iter(f"{{{namespace}}}t"))


def docx_has_line_numbers(path: Path) -> bool:
    if not path.exists():
        return False
    return "w:lnNumType" in docx_xml(path, "word/document.xml")


def docx_has_page_field(path: Path) -> bool:
    if not path.exists():
        return False
    with zipfile.ZipFile(path) as zf:
        footer_names = [name for name in zf.namelist() if name.startswith("word/footer")]
        for name in footer_names:
            if "PAGE" in zf.read(name).decode("utf-8", errors="ignore"):
                return True
    return False


def wide_table_caption_issues(path: Path) -> list[int]:
    """Return 1-based wide-table indices whose caption is not the previous body item."""
    if not path.exists():
        return []
    namespace = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    root = ET.fromstring(docx_xml(path, "word/document.xml"))
    body = root.find(f"{{{namespace}}}body")
    if body is None:
        return []
    children = list(body)
    issues: list[int] = []
    wide_index = 0
    for index, child in enumerate(children):
        if child.tag != f"{{{namespace}}}tbl":
            continue
        grid = child.find(f"{{{namespace}}}tblGrid")
        column_count = len(grid) if grid is not None else 0
        if column_count < 9:
            continue
        wide_index += 1
        if index == 0:
            issues.append(wide_index)
            continue
        previous_text = "".join(
            node.text or ""
            for node in children[index - 1].iter(f"{{{namespace}}}t")
        ).strip()
        if not re.match(r"^(?:Supplementary Table S\d+|Table \d+)\.", previous_text):
            issues.append(wide_index)
    return issues


def markdown_equations(path: Path) -> list[tuple[str, str]]:
    matches = re.findall(
        r"```equation\s*\nnumber:\s*([^\s]+)\s*\n(.+?)```",
        read_text(path),
        flags=re.S,
    )
    return [
        (number, " ".join(line.strip() for line in latex.splitlines() if line.strip()))
        for number, latex in matches
    ]


def markdown_equation_numbers(path: Path) -> list[str]:
    return [number for number, _latex in markdown_equations(path)]


def numbered_equation_texts(path: Path) -> dict[str, str]:
    w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    m = "http://schemas.openxmlformats.org/officeDocument/2006/math"
    root = ET.fromstring(docx_xml(path, "word/document.xml"))
    equations: dict[str, str] = {}
    for table in root.iter(f"{{{w}}}tbl"):
        formula = table.find(f".//{{{m}}}oMath")
        row = table.find(f"{{{w}}}tr")
        if formula is None or row is None:
            continue
        cells = row.findall(f"{{{w}}}tc")
        if len(cells) != 3:
            continue
        label = "".join(node.text or "" for node in cells[-1].iter(f"{{{w}}}t"))
        match = re.fullmatch(r"\(((?:S)?\d+)\)", label)
        if match:
            equations[match.group(1)] = "".join(
                node.text or "" for node in formula.iter(f"{{{m}}}t")
            )
    return equations


def expected_equation_text(latex: str) -> str:
    m = "http://schemas.openxmlformats.org/officeDocument/2006/math"
    return "".join(node.text or "" for node in omml_element(latex).iter(f"{{{m}}}t"))


def equation_layout_issues(path: Path) -> list[str]:
    """Return numbered equations whose Word table geometry can clip in WPS."""
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))

    w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    m = "http://schemas.openxmlformats.org/officeDocument/2006/math"
    issues: list[str] = []
    for table in root.iter(f"{{{w}}}tbl"):
        if table.find(f".//{{{m}}}oMath") is None:
            continue
        row = table.find(f"{{{w}}}tr")
        if row is None:
            continue
        cells = row.findall(f"{{{w}}}tc")
        if len(cells) != 3:
            continue
        number = "".join(node.text or "" for node in cells[-1].iter(f"{{{w}}}t"))
        if not re.fullmatch(r"\((?:S)?\d+\)", number):
            continue

        grid = table.find(f"{{{w}}}tblGrid")
        grid_widths = (
            [
                int(node.get(f"{{{w}}}w", "0"))
                for node in grid.findall(f"{{{w}}}gridCol")
            ]
            if grid is not None
            else []
        )
        cell_widths = []
        for cell in cells:
            width = cell.find(f"./{{{w}}}tcPr/{{{w}}}tcW")
            cell_widths.append(int(width.get(f"{{{w}}}w", "0")) if width is not None else 0)

        table_width_node = table.find(f"./{{{w}}}tblPr/{{{w}}}tblW")
        table_width = (
            int(table_width_node.get(f"{{{w}}}w", "0"))
            if table_width_node is not None
            else 0
        )
        fixed = table.find(f"./{{{w}}}tblPr/{{{w}}}tblLayout")
        cant_split = row.find(f"./{{{w}}}trPr/{{{w}}}cantSplit")

        if (
            len(grid_widths) != 3
            or grid_widths != cell_widths
            or grid_widths[1] < 8000
            or table_width != sum(grid_widths)
            or fixed is None
            or fixed.get(f"{{{w}}}type") != "fixed"
            or cant_split is None
        ):
            issues.append(number)
    return issues


def image_refs(markdown_path: Path) -> list[Path]:
    refs = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", read_text(markdown_path))
    return [(markdown_path.parent / ref).resolve() for ref in refs]


def image_dpi(path: Path) -> tuple[float, float] | None:
    try:
        with Image.open(path) as im:
            dpi = im.info.get("dpi")
    except Exception:
        return None
    if not dpi:
        return None
    return float(dpi[0]), float(dpi[1])


def tiff_compression(path: Path) -> int | None:
    try:
        with Image.open(path) as im:
            return int(im.tag_v2.get(259))
    except Exception:
        return None


def check_files() -> list[Check]:
    checks: list[Check] = []
    for path in REQUIRED_MD + REQUIRED_DOCX + REQUIRED_PDF:
        if path.exists():
            checks.append(Check("PASS", f"required file exists: {path.relative_to(ROOT)}"))
        else:
            checks.append(Check("FAIL", f"required file missing: {path.relative_to(ROOT)}"))
    if SKIP_PDF:
        checks.append(Check("PASS", "supplementary PDF intentionally deferred (BMB_SKIP_PDF=1)"))
    return checks


def check_main_markdown() -> list[Check]:
    checks: list[Check] = []
    text = read_text(MAIN_EN)
    body, references = split_main_references(text)
    abstract = extract_section(text, "Abstract")
    words = count_words(abstract)
    if 150 <= words <= 250:
        checks.append(Check("PASS", f"abstract length is {words} words (BMB target 150-250)"))
    else:
        checks.append(Check("FAIL", f"abstract length is {words} words, outside 150-250"))

    # mRNA is treated as a standard field abbreviation and need not be expanded.
    required_abstract_definitions = {
        "PC1": "principal component (PC1)",
    }
    for abbreviation, definition in required_abstract_definitions.items():
        if abbreviation in abstract and definition not in abstract:
            checks.append(Check("WARN", f"abstract uses {abbreviation} without preferred definition `{definition}`"))
        else:
            checks.append(Check("PASS", f"abstract abbreviation definition ok: {abbreviation}"))

    match = re.search(r"^\*\*Keywords:\*\*\s*(.+)$", text, flags=re.M)
    keywords = [k.strip() for k in re.split(r";|,", match.group(1)) if k.strip()] if match else []
    if 4 <= len(keywords) <= 6:
        checks.append(Check("PASS", f"keyword count is {len(keywords)}"))
    else:
        checks.append(Check("FAIL", f"keyword count is {len(keywords)}, expected 4-6"))

    for heading in REQUIRED_STATEMENT_HEADINGS:
        if re.search(rf"^###\s+{re.escape(heading)}\s*$", text, flags=re.M):
            checks.append(Check("PASS", f"main manuscript statement present: {heading}"))
        else:
            checks.append(Check("FAIL", f"main manuscript statement missing: {heading}"))

    numeric_citations = re.findall(r"\[[0-9,\s;\-]+\]", text)
    if numeric_citations:
        checks.append(Check("FAIL", f"English main manuscript has numeric citation remnants: {numeric_citations[:5]}"))
    else:
        checks.append(Check("PASS", "English main manuscript has no numeric citation remnants"))

    if "Online Resource 1" in body:
        checks.append(Check("PASS", "main manuscript cites supplementary material as Online Resource 1"))
    else:
        checks.append(Check("WARN", "main manuscript does not cite supplementary material as Online Resource 1"))

    figure_numbers = [int(value) for value in re.findall(r"^\*\*Fig\.\s+(\d+)\*\*", text, flags=re.M)]
    if figure_numbers == list(range(1, 9)):
        checks.append(Check("PASS", "English main figure captions are contiguous: Fig. 1-8"))
    else:
        checks.append(Check("FAIL", f"English main figure captions are not Fig. 1-8: {figure_numbers}"))

    caption_lines = re.findall(r"^\*\*Fig\.\s+\d+\*\*.*$", text, flags=re.M)
    punctuated = [line for line in caption_lines if line.rstrip().endswith(tuple(".,;:!?"))]
    if punctuated:
        checks.append(Check("FAIL", "English main figure captions end with punctuation"))
    else:
        checks.append(Check("PASS", "English main figure captions follow BMB terminal-punctuation style"))

    table_numbers = [int(value) for value in re.findall(r"^\*\*Table\s+(\d+)\*\*", text, flags=re.M)]
    table_count = count_markdown_tables(text)
    if table_numbers == list(range(1, table_count + 1)):
        checks.append(Check("PASS", f"English main tables are uniquely captioned and contiguous: Table 1-{table_count}"))
    else:
        checks.append(
            Check(
                "FAIL",
                f"English main table captions do not match {table_count} Markdown tables: {table_numbers}",
            )
        )

    required_universes = [
        "16,444",
        "16,951",
        "15,788",
        "14,244",
        "13,265",
        "12,592",
        "13,532",
        "12,916",
        "11,107",
        "12,644",
        "12,307",
        "10,768",
        "10,682",
    ]
    missing_universes = [value for value in required_universes if value not in body]
    if missing_universes:
        checks.append(Check("FAIL", f"canonical analysis-universe counts missing from main text: {missing_universes}"))
    else:
        checks.append(Check("PASS", "canonical analysis-universe counts are present in the main text"))

    ref_keys = [reference_family_year(line) for line in references]
    parsed_refs = [(key, line) for key, line in zip(ref_keys, references) if key is not None]
    cited = citation_family_years(body)
    uncited = [line for key, line in parsed_refs if key not in cited]
    if uncited:
        checks.append(Check("FAIL", f"uncited references in English main manuscript: {uncited[:5]}"))
    else:
        checks.append(Check("PASS", "all parsed English main references are cited in the main text"))

    sort_keys = [key[0].lower() for key, _line in parsed_refs]
    if sort_keys == sorted(sort_keys):
        checks.append(Check("PASS", "English main references are alphabetized by first author"))
    else:
        checks.append(Check("FAIL", "English main references are not alphabetized by first author"))

    if "Generative AI tools were used" in body:
        checks.append(Check("PASS", "generative-AI use is disclosed in the Methods"))
    else:
        checks.append(Check("FAIL", "generative-AI use disclosure is missing from the Methods"))
    return checks


def check_supplement() -> list[Check]:
    checks: list[Check] = []
    text = read_text(SUPP_EN)
    body, references = split_supplement_references(text)
    for marker in SI_METADATA_MARKERS:
        if marker in text:
            checks.append(Check("PASS", f"SI metadata present: {marker.strip('*:')}"))
        else:
            checks.append(Check("FAIL", f"SI metadata missing: {marker}"))

    numeric_citations = re.findall(r"\[[0-9,\s;\-]+\]", text)
    if numeric_citations:
        checks.append(Check("FAIL", f"English SI has numeric citation remnants: {numeric_citations[:5]}"))
    else:
        checks.append(Check("PASS", "English SI has no numeric citation remnants"))

    figure_numbers = [int(value) for value in re.findall(r"^\*\*Supplementary Fig\. S(\d+)\*\*", text, flags=re.M)]
    if figure_numbers == list(range(1, 5)):
        checks.append(Check("PASS", "English SI figure captions are contiguous: Supplementary Fig. S1-S4"))
    else:
        checks.append(Check("FAIL", f"English SI figure captions are not Supplementary Fig. S1-S4: {figure_numbers}"))

    caption_lines = re.findall(r"^\*\*Supplementary Fig\.\s+S\d+\*\*.*$", text, flags=re.M)
    punctuated = [line for line in caption_lines if line.rstrip().endswith(tuple(".,;:!?"))]
    if punctuated:
        checks.append(Check("FAIL", "English SI figure captions end with punctuation"))
    else:
        checks.append(Check("PASS", "English SI figure captions follow BMB terminal-punctuation style"))

    table_numbers = [
        int(value)
        for value in re.findall(r"^###\s+Supplementary Table S(\d+)\.", text, flags=re.M)
    ]
    table_count = count_markdown_tables(text)
    if table_numbers == list(range(1, table_count + 1)):
        checks.append(
            Check(
                "PASS",
                f"English SI tables are uniquely captioned and contiguous: Supplementary Table S1-S{table_count}",
            )
        )
    else:
        checks.append(
            Check(
                "FAIL",
                f"English SI table captions do not match {table_count} Markdown tables: {table_numbers}",
            )
        )

    reader_facing_table_lines = "\n".join(
        line for line in text.splitlines() if line.lstrip().startswith("|")
    )
    machine_table_tokens = [
        "compact_all_plus_",
        "moesm3_no_gejman_",
        "moesm3_all_samples_",
        "Saluki_human_PC1",
        "prior_only_linear",
        "prior_plus_compact_residual",
        "pearson_mean_sd",
        "target_vs_Saluki_pearson",
        "all_one2one",
    ]
    remaining_machine_tokens = [
        token for token in machine_table_tokens if token in reader_facing_table_lines
    ]
    if remaining_machine_tokens:
        checks.append(
            Check(
                "FAIL",
                f"reader-facing SI tables retain machine identifiers: {remaining_machine_tokens}",
            )
        )
    else:
        checks.append(Check("PASS", "reader-facing SI tables use readable display labels"))

    ref_keys = [reference_family_year(line) for line in references]
    parsed_refs = [(key, line) for key, line in zip(ref_keys, references) if key is not None]
    cited = citation_family_years(body)
    uncited = [line for key, line in parsed_refs if key not in cited]
    if uncited:
        checks.append(Check("FAIL", f"uncited references in English SI: {uncited[:5]}"))
    else:
        checks.append(Check("PASS", "all parsed English SI references are cited in the supplement"))

    sort_keys = [key[0].lower() for key, _line in parsed_refs]
    if sort_keys == sorted(sort_keys):
        checks.append(Check("PASS", "English SI references are alphabetized by first author"))
    else:
        checks.append(Check("FAIL", "English SI references are not alphabetized by first author"))
    return checks


def check_docx() -> list[Check]:
    checks: list[Check] = []
    if docx_has_line_numbers(MAIN_DOCX):
        checks.append(Check("PASS", "main DOCX has continuous line-number setting"))
    else:
        checks.append(Check("FAIL", "main DOCX lacks line-number setting"))

    for path in (MAIN_DOCX, SUPP_DOCX):
        if docx_has_page_field(path):
            checks.append(Check("PASS", f"{path.name} has PAGE footer field"))
        else:
            checks.append(Check("FAIL", f"{path.name} lacks PAGE footer field"))

    caption_issues = wide_table_caption_issues(SUPP_DOCX)
    if caption_issues:
        checks.append(
            Check(
                "FAIL",
                f"wide supplementary tables are separated from their captions: {caption_issues}",
            )
        )
    else:
        checks.append(Check("PASS", "wide supplementary tables follow their captions without a section-break gap"))

    main_visible_text = docx_visible_text(MAIN_DOCX)
    title_visible_text = docx_visible_text(SUB / "bmb_title_page_en.docx")
    if "Keywords:*" in main_visible_text:
        checks.append(Check("FAIL", "main DOCX contains a stray asterisk after the Keywords label"))
    else:
        checks.append(Check("PASS", "main DOCX Keywords label has no stray author-marker asterisk"))
    if (
        "Authors: Wenzhuo Wang, Ying Shao" in main_visible_text
        and "Wenzhuo Wang, Ying Shao" in title_visible_text
        and "Ying Shao*" not in main_visible_text
        and "Ying Shao*" not in title_visible_text
    ):
        checks.append(Check("PASS", "single-affiliation author formatting is consistent in main DOCX and title page"))
    else:
        checks.append(Check("FAIL", "author markers are inconsistent between main DOCX and title page"))

    equation_specs = {
        MAIN_EN: (MAIN_DOCX, [str(index) for index in range(1, 10)]),
        MAIN_CN: (SUB / "bmb_main_manuscript_cn.docx", [str(index) for index in range(1, 10)]),
        SUPP_EN: (SUPP_DOCX, ["S1", "S2"]),
        SUPP_CN: (SUB / "bmb_supplementary_material_cn.docx", ["S1", "S2"]),
    }
    for markdown_path, (docx_path, expected_numbers) in equation_specs.items():
        source_numbers = markdown_equation_numbers(markdown_path)
        if source_numbers == expected_numbers:
            checks.append(
                Check("PASS", f"equation numbering is contiguous in {markdown_path.name}: {source_numbers}")
            )
        else:
            checks.append(
                Check(
                    "FAIL",
                    f"equation numbering mismatch in {markdown_path.name}: "
                    f"expected {expected_numbers}, found {source_numbers}",
                )
            )

        xml = docx_xml(docx_path, "word/document.xml")
        math_count = xml.count("<m:oMath")
        missing_labels = [number for number in source_numbers if f">({number})<" not in xml]
        if math_count >= len(source_numbers) and not missing_labels and "number:" not in xml:
            checks.append(
                Check(
                    "PASS",
                    f"editable OMML equations present in {docx_path.name}: "
                    f"{math_count} math objects, {len(source_numbers)} numbered",
                )
            )
        else:
            checks.append(
                Check(
                    "FAIL",
                    f"editable-equation audit failed for {docx_path.name}: "
                    f"math_objects={math_count}, missing_labels={missing_labels}",
                )
            )

        expected_texts = {
            number: expected_equation_text(latex)
            for number, latex in markdown_equations(markdown_path)
        }
        actual_texts = numbered_equation_texts(docx_path)
        content_mismatches = [
            number
            for number, expected_text in expected_texts.items()
            if actual_texts.get(number) != expected_text
        ]
        if content_mismatches:
            checks.append(
                Check(
                    "FAIL",
                    f"numbered-equation content mismatch in {docx_path.name}: {content_mismatches}",
                )
            )
        else:
            checks.append(
                Check(
                    "PASS",
                    f"all numbered equations retain their complete source content in {docx_path.name}",
                )
            )

        layout_issues = equation_layout_issues(docx_path)
        if layout_issues:
            checks.append(
                Check(
                    "FAIL",
                    f"WPS-safe equation layout failed in {docx_path.name}: {layout_issues}",
                )
            )
        else:
            checks.append(
                Check(
                    "PASS",
                    f"numbered equations have synchronized full-width WPS-safe layout in {docx_path.name}",
                )
            )
    return checks


def check_artifact_freshness() -> list[Check]:
    checks: list[Check] = []
    source_groups = {
        MAIN_DOCX: [MAIN_EN, *image_refs(MAIN_EN)],
        SUPP_DOCX: [SUPP_EN, *image_refs(SUPP_EN)],
        SUB / "bmb_main_manuscript_cn.docx": [SUB / "bmb_main_manuscript_cn.md", *image_refs(SUB / "bmb_main_manuscript_cn.md")],
        SUB / "bmb_supplementary_material_cn.docx": [SUB / "bmb_supplementary_material_cn.md", *image_refs(SUB / "bmb_supplementary_material_cn.md")],
        SUB / "bmb_title_page_en.docx": [TITLE_EN],
        SUB / "bmb_title_page_cn.docx": [SUB / "bmb_title_page_cn.md"],
        SUB / "bmb_cover_letter_en.docx": [COVER_EN],
        SUB / "bmb_cover_letter_cn.docx": [SUB / "bmb_cover_letter_cn.md"],
        SUB / "bmb_statements_and_declarations_en.docx": [DECL_EN],
        SUB / "bmb_statements_and_declarations_cn.docx": [SUB / "bmb_statements_and_declarations_cn.md"],
        SUB / "bmb_author_contributions_cn.docx": [SUB / "bmb_author_contributions_cn.md"],
        SUB / "bmb_editor_highlights_cn.docx": [SUB / "bmb_editor_highlights_cn.md"],
    }
    for artifact, sources in source_groups.items():
        if not artifact.exists():
            continue
        newest_source = max(source.stat().st_mtime for source in sources if source.exists())
        if artifact.stat().st_mtime >= newest_source:
            checks.append(Check("PASS", f"artifact is current relative to Markdown and figures: {artifact.name}"))
        else:
            checks.append(Check("FAIL", f"artifact is older than its Markdown or figure source: {artifact.name}"))

    if not SKIP_PDF and SUPP_PDF.exists() and SUPP_DOCX.exists():
        if SUPP_PDF.stat().st_mtime >= SUPP_DOCX.stat().st_mtime:
            checks.append(Check("PASS", "supplementary PDF is current relative to the DOCX"))
        else:
            checks.append(
                Check(
                    "WARN",
                    "supplementary PDF is older than the DOCX; regenerate only after the manuscript is frozen",
                )
            )

    lock_files = sorted(SUB.glob("~$*.docx"))
    if lock_files:
        checks.append(Check("WARN", f"temporary Office lock files remain: {[path.name for path in lock_files]}"))
    else:
        checks.append(Check("PASS", "no temporary Office lock files remain in the submission directory"))
    return checks


def check_metadata_and_reference_files() -> list[Check]:
    checks: list[Check] = []
    expected_affiliations = ("School of Science, Dalian Maritime University",)
    for path in (MAIN_EN, SUPP_EN, TITLE_EN):
        text = read_text(path)
        missing = [affiliation for affiliation in expected_affiliations if affiliation not in text]
        if missing:
            checks.append(Check("FAIL", f"formal affiliation wording missing from {path.name}: {missing}"))
        else:
            checks.append(Check("PASS", f"formal affiliation wording matches in {path.name}"))

    expected_title = (
        "Study-aware label auditing, cross-species transfer, and ortholog-informed "
        "target shrinkage for mammalian mRNA half-life prediction"
    )
    for path in (MAIN_EN, SUPP_EN, TITLE_EN, COVER_EN):
        if expected_title in read_text(path):
            checks.append(Check("PASS", f"manuscript title matches in {path.name}"))
        else:
            checks.append(Check("FAIL", f"manuscript title mismatch in {path.name}"))

    author_names = ("Wenzhuo Wang", "Ying Shao")
    for path in (MAIN_EN, SUPP_EN, TITLE_EN, DECL_EN):
        missing = [name for name in author_names if name not in read_text(path)]
        if missing:
            checks.append(Check("FAIL", f"author metadata missing from {path.name}: {missing}"))
        else:
            checks.append(Check("PASS", f"both author names are present in {path.name}"))

    ying_shao_orcid = "0000-0002-4056-5757"
    for path in (MAIN_EN, SUPP_EN, TITLE_EN, DECL_EN, CITATION):
        if ying_shao_orcid in read_text(path):
            checks.append(Check("PASS", f"Ying Shao ORCID matches in {path.name}"))
        else:
            checks.append(Check("FAIL", f"Ying Shao ORCID missing from {path.name}"))

    corresponding_author_email = "yshao@dlmu.edu.cn"
    for path in (MAIN_EN, SUPP_EN, TITLE_EN, COVER_EN, DECL_EN, CITATION):
        text = read_text(path)
        if any(email in text for email in ("myc@lnnu.edu.cn", "dlzhang@dicp.ac.cn")):
            checks.append(Check("FAIL", f"superseded corresponding-author email remains in {path.name}"))
        elif corresponding_author_email not in text:
            checks.append(Check("FAIL", f"Ying Shao email missing from {path.name}"))
        else:
            checks.append(Check("PASS", f"Ying Shao email matches in {path.name}"))

    funding_statement = "No specific funding was received for this work."
    for path in (MAIN_EN, TITLE_EN, DECL_EN):
        text = read_text(path)
        if funding_statement not in text:
            checks.append(Check("FAIL", f"no-funding statement missing from {path.name}"))
        elif any(grant in text for grant in ("22373101", "22203089", "22573043")):
            checks.append(Check("FAIL", f"superseded funding identifier remains in {path.name}"))
        else:
            checks.append(Check("PASS", f"no-funding statement matches in {path.name}"))

    repository_url = "https://github.com/wwzdl/mrna-pc1-label"
    for path in (MAIN_EN, SUPP_EN, TITLE_EN, COVER_EN, DECL_EN):
        if repository_url in read_text(path):
            checks.append(Check("PASS", f"repository URL matches in {path.name}"))
        else:
            checks.append(Check("FAIL", f"repository URL missing from {path.name}"))

    bib_text = read_text(BIB)
    ris_text = read_text(RIS)
    bib_ids = re.findall(r"^@\w+\{([^,]+),", bib_text, flags=re.M)
    ris_ids = re.findall(r"^ID  - (.+)$", ris_text, flags=re.M)
    if len(bib_ids) != len(set(bib_ids)):
        checks.append(Check("FAIL", "BibTeX library contains duplicate citation keys"))
    else:
        checks.append(Check("PASS", f"BibTeX library contains {len(bib_ids)} unique citation keys"))
    if len(ris_ids) != len(set(ris_ids)):
        checks.append(Check("FAIL", "RIS library contains duplicate citation IDs"))
    elif set(ris_ids) != set(bib_ids):
        missing_from_ris = sorted(set(bib_ids) - set(ris_ids))
        extra_in_ris = sorted(set(ris_ids) - set(bib_ids))
        checks.append(
            Check(
                "FAIL",
                f"BibTeX/RIS citation-ID mismatch: missing={missing_from_ris}, extra={extra_in_ris}",
            )
        )
    else:
        checks.append(Check("PASS", f"BibTeX and RIS libraries contain the same {len(bib_ids)} records"))

    main_ref_keys = {
        key
        for key in (reference_family_year(line) for line in split_main_references(read_text(MAIN_EN))[1])
        if key is not None
    }
    supp_ref_keys = {
        key
        for key in (reference_family_year(line) for line in split_supplement_references(read_text(SUPP_EN))[1])
        if key is not None
    }
    active_ref_count = len(main_ref_keys | supp_ref_keys)
    if len(bib_ids) == active_ref_count:
        checks.append(Check("PASS", f"reference-manager library size matches {active_ref_count} active manuscript/SI references"))
    else:
        checks.append(
            Check(
                "FAIL",
                f"reference-manager library has {len(bib_ids)} records but manuscripts use {active_ref_count} unique references",
            )
        )
    return checks


def check_quantitative_claims() -> list[Check]:
    checks: list[Check] = []
    manuscript = read_text(MAIN_EN)
    supplement = read_text(SUPP_EN)

    def require_tokens_in(label: str, tokens: list[str], text: str) -> None:
        missing = [token for token in tokens if token not in text]
        if missing:
            checks.append(Check("FAIL", f"quantitative claim mismatch for {label}: missing {missing}"))
        else:
            checks.append(Check("PASS", f"quantitative claim matches result table: {label}"))

    def require_tokens(label: str, tokens: list[str]) -> None:
        require_tokens_in(f"main text, {label}", tokens, manuscript)

    feature_rows = read_tsv(ROOT / "results" / "human_compact_regulatory_benchmark.tsv")
    feature_dimensions = {
        row["setting"]: int(row["n_features"])
        for row in feature_rows
    }
    expected_feature_dimensions = {
        "compact_base": 526,
        "compact_cwcs": 845,
        "compact_seqweaver": 1306,
        "compact_deepripe": 703,
        "compact_all": 1802,
    }
    if feature_dimensions == expected_feature_dimensions:
        checks.append(Check("PASS", "compact feature-block dimensions match the benchmark table"))
    else:
        checks.append(
            Check(
                "FAIL",
                f"compact feature-block dimensions mismatch: {feature_dimensions}",
            )
        )
    block_dimensions = {
        "CWCS": feature_dimensions["compact_cwcs"] - feature_dimensions["compact_base"],
        "SeqWeaver": feature_dimensions["compact_seqweaver"] - feature_dimensions["compact_base"],
        "DeepRiPe": feature_dimensions["compact_deepripe"] - feature_dimensions["compact_base"],
    }
    require_tokens(
        "compact_all feature definition",
        [
            "1,802-dimensional",
            "526 sequence features",
            f"{block_dimensions['CWCS']} CWCS",
            f"{block_dimensions['SeqWeaver']} SeqWeaver",
            f"{block_dimensions['DeepRiPe']} DeepRiPe",
            "not the name of a model reported by Saluki",
        ],
    )

    bootstrap = read_tsv(ROOT / "results" / "key_claim_bootstrap" / "paired_bootstrap_summary.tsv")
    for claim in (
        "Saluki agreement after Gejman exclusion",
        "Ortholog concordance after Gejman exclusion: all one-to-one orthologs",
        "Real prior versus shuffled prior",
    ):
        row = select_row(bootstrap, claim=claim)
        require_tokens(
            claim,
            [
                f"{float(row['difference_a_minus_b']):.4f}",
                f"{float(row['difference_ci95_low']):.4f}",
                f"{float(row['difference_ci95_high']):.4f}",
            ],
        )

    cross_target = read_tsv(
        ROOT / "results" / "ortholog_regularized_cross_target" / "summary_aggregate.tsv"
    )
    for training_target in (
        "human_no_gejman_pc1",
        "orthoreg_reconstructed_mouse_pc1_0.10",
    ):
        row = select_row(
            cross_target,
            training_target=training_target,
            evaluation_target="human_no_gejman_pc1",
        )
        require_tokens(
            f"cross-target {training_target}",
            [f"{float(row['pearson_mean']):.4f}", f"{float(row['pearson_sd']):.4f}"],
        )

    distances = read_tsv(ROOT / "results" / "ortholog_label_distance" / "distance_summary.tsv")
    distance = select_row(
        distances,
        comparison="orthoreg_lambda0p10_vs_human_no_gejman_pc1",
        subset="all_one2one",
    )
    require_tokens(
        "0.10 target geometry",
        [
            f"{float(distance['pearson']):.3f}",
            f"{float(distance['rmse']):.3f}",
            f"{float(distance['mae']):.3f}",
        ],
    )

    residual = read_tsv(ROOT / "results" / "prior_residual_analysis" / "summary.tsv")
    prior = select_row(residual, model="prior_only_linear")
    final = select_row(residual, model="prior_plus_compact_residual")
    require_tokens(
        "prior residual decomposition",
        [
            f"{float(prior['pearson']):.3f}",
            f"{float(prior['r2']):.3f}",
            f"{float(final['pearson']):.3f}",
            f"{float(final['r2']):.3f}",
            f"{100.0 * float(final['remaining_variance_explained']):.1f}%",
        ],
    )

    ablation = read_tsv(
        ROOT
        / "results"
        / "ortholog_prior_ablation_lambda0p1_10fold"
        / "summary_by_setting.tsv"
    )
    ablation_settings = (
        "human_only",
        "mask_only",
        "no_direct_reconstructed_mouse_pc1",
        "exact_target_and_saluki_priors",
        "primary_audit_and_saluki_priors",
    )
    ablation_tokens: list[str] = []
    for setting in ablation_settings:
        row = select_row(ablation, setting=setting)
        ablation_tokens.extend(
            [
                f"{float(row['pearson_vs_target_mean']):.4f}",
                f"{float(row['pearson_vs_target_sd']):.4f}",
            ]
        )
    strict_no_direct = select_row(
        ablation,
        setting="no_direct_reconstructed_mouse_pc1",
    )
    strict_exact = select_row(
        ablation,
        setting="exact_target_and_saluki_priors",
    )
    direct_increment = (
        float(strict_exact["pearson_vs_target_mean"])
        - float(strict_no_direct["pearson_vs_target_mean"])
    )
    ablation_tokens.append(f"{direct_increment:.4f}")
    require_tokens_in("main text, strict target-vector ablation", ablation_tokens, manuscript)
    require_tokens_in("supplement, strict target-vector ablation", ablation_tokens, supplement)

    ablation_metadata_path = (
        ROOT
        / "results"
        / "ortholog_prior_ablation_lambda0p1_10fold"
        / "analysis_metadata.json"
    )
    with ablation_metadata_path.open(encoding="utf-8") as handle:
        ablation_metadata = json.load(handle)
    expected_ablation_metadata = {
        "target_human_coverage_threshold": 10,
        "target_mouse_coverage_threshold": 5,
        "direct_prior_value_column": "target_reconstructed_mouse_pc1",
        "primary_target_prediction_setting": "primary_audit_and_saluki_priors",
    }
    metadata_mismatches = {
        key: (ablation_metadata.get(key), expected)
        for key, expected in expected_ablation_metadata.items()
        if ablation_metadata.get(key) != expected
    }
    if metadata_mismatches:
        checks.append(
            Check(
                "FAIL",
                f"strict target-vector ablation metadata mismatch: {metadata_mismatches}",
            )
        )
    else:
        checks.append(Check("PASS", "strict target-vector ablation metadata is explicit and current"))

    influence_comparators = read_tsv(
        ROOT / "results" / "study_influence_sensitivity" / "conditional_same_size_summary.tsv"
    )
    for metric in (
        "pc1_stability_pearson",
        "delta_saluki_pearson",
        "delta_ortholog_pearson",
    ):
        row = select_row(influence_comparators, metric=metric)
        display_precision = 3 if metric == "pc1_stability_pearson" else 4
        require_tokens(
            f"conditional same-size study-influence comparison: {metric}",
            [
                f"{float(row['observed_gejman']):.{display_precision}f}",
                f"{float(row['comparator_median']):.{display_precision}f}",
                f"{float(row['empirical_tail_proportion']):.3f}",
            ],
        )

    sensitivity = read_tsv(
        ROOT / "results" / "study_influence_sensitivity" / "pipeline_sensitivity_summary.tsv"
    )
    if len(sensitivity) == 12 and all(row["worst_left_out"] == "Gejman" for row in sensitivity):
        checks.append(
            Check(
                "PASS",
                "study-influence sensitivity covers 12 prespecified settings and ranks Gejman first in each",
            )
        )
    else:
        checks.append(
            Check(
                "FAIL",
                "study-influence sensitivity table is incomplete or does not support the stated 12-setting result",
            )
        )

    weighted = read_tsv(
        ROOT / "results" / "study_influence_sensitivity" / "study_weighted_summary.tsv"
    )
    for estimator in ("equal_study_weighted", "study_mean_collapsed"):
        row = select_row(weighted, estimator=estimator)
        if row["gejman_rank"] == "3":
            checks.append(Check("PASS", f"{estimator} Gejman rank matches manuscript: 3"))
        else:
            checks.append(
                Check(
                    "FAIL",
                    f"{estimator} Gejman rank is {row['gejman_rank']}, expected 3",
                )
            )

    provenance_path = ROOT / "results" / "human_xgboost_parameter_provenance.json"
    with provenance_path.open(encoding="utf-8") as handle:
        provenance = json.load(handle)
    expected_parameters = {
        "max_depth": 6,
        "learning_rate": 0.02,
        "min_child_weight": 4,
        "subsample": 0.9,
        "colsample_bytree": 0.75,
        "reg_alpha": 0.0,
        "reg_lambda": 1.0,
    }
    if (
        provenance.get("genes") == 12916
        and provenance.get("features") == 526
        and provenance.get("selected_parameters") == expected_parameters
        and "same compendium" in provenance.get("inference_boundary", "")
    ):
        checks.append(Check("PASS", "exploratory XGBoost parameter selection has machine-readable provenance"))
    else:
        checks.append(Check("FAIL", "XGBoost parameter-selection provenance is incomplete or inconsistent"))

    fig1_rows = read_tsv(
        SUB / "figures" / "data" / "Fig01_panel_data.tsv"
    )
    fig1_values = {
        (row["section"], row["item"], row["metric"]): float(row["value"])
        for row in fig1_rows
    }
    global_summary = read_tsv(
        ROOT / "results" / "tenfold_global_prior" / "summary_by_setting.tsv"
    )
    distance_rows = read_tsv(
        ROOT / "results" / "ortholog_label_distance" / "distance_summary.tsv"
    )
    target_rows = read_tsv(
        ROOT
        / "results"
        / "ortholog_regularized_label_10fold_main"
        / "summary_by_label.tsv"
    )
    conditional_rows = {
        row["metric"]: row
        for row in influence_comparators
    }
    global_rows = {
        row["prior_mode"]: row
        for row in global_summary
    }
    distance_row = select_row(
        distance_rows,
        comparison="orthoreg_lambda0p10_vs_human_no_gejman_pc1",
        subset="all_one2one",
    )
    target_row = select_row(
        target_rows,
        label="orthoreg_reconstructed_mouse_pc1_0.1",
    )
    expected_fig1 = {
        ("A", "study_audit", "saluki_delta_r"): float(
            conditional_rows["delta_saluki_pearson"]["observed_gejman"]
        ),
        ("A", "study_audit", "ortholog_delta_r"): float(
            conditional_rows["delta_ortholog_pearson"]["observed_gejman"]
        ),
        ("A", "study_audit", "geometry_tail_proportion"): float(
            conditional_rows["pc1_stability_pearson"]["empirical_tail_proportion"]
        ),
        ("B", "human_only", "pearson_mean"): float(global_rows["none"]["pearson_mean"]),
        ("B", "fixed_target", "genes"): float(global_rows["none"]["n"]),
        ("B", "prior_enhanced", "pearson_mean"): float(global_rows["real"]["pearson_mean"]),
        ("B", "prior_shuffled", "pearson_mean"): float(global_rows["shuffled"]["pearson_mean"]),
        ("C", "target_geometry", "human_pearson"): float(distance_row["pearson"]),
        ("C", "target_geometry", "shift_rmse"): float(distance_row["rmse"]),
        ("C", "target_geometry", "ortholog_pairs"): float(distance_row["n"]),
        ("C", "target_prediction", "genes"): float(target_row["n"]),
    }
    fig1_mismatches = {
        key: (fig1_values.get(key), expected)
        for key, expected in expected_fig1.items()
        if key not in fig1_values or abs(fig1_values[key] - expected) > 1e-12
    }
    if set(fig1_values) != set(expected_fig1):
        fig1_mismatches["row_keys"] = (sorted(fig1_values), sorted(expected_fig1))
    missing_fig1_sources = [
        row["source_file"]
        for row in fig1_rows
        if not (ROOT / row["source_file"]).exists()
    ]
    if fig1_mismatches or missing_fig1_sources:
        checks.append(
            Check(
                "FAIL",
                "Fig. 1 panel-data provenance mismatch: "
                f"values={fig1_mismatches}, missing_sources={missing_fig1_sources}",
            )
        )
    else:
        checks.append(Check("PASS", "Fig. 1 panel data match their declared result-table sources"))
    return checks


def check_figures() -> list[Check]:
    checks: list[Check] = []
    refs = image_refs(MAIN_EN) + image_refs(SUPP_EN)
    for path in refs:
        if not path.exists():
            checks.append(Check("FAIL", f"referenced figure missing: {path}"))
            continue
        checks.append(Check("PASS", f"referenced figure exists: {path.relative_to(SUB)}"))
        if path.suffix.lower() in {".png", ".tif", ".tiff"}:
            dpi = image_dpi(path)
            if dpi is None:
                checks.append(Check("WARN", f"no DPI metadata found for referenced figure: {path.relative_to(SUB)}"))
            elif min(dpi) >= 300:
                checks.append(Check("PASS", f"referenced figure DPI metadata ok: {path.relative_to(SUB)} {dpi[0]:.0f}x{dpi[1]:.0f}"))
            else:
                checks.append(Check("WARN", f"low DPI metadata for referenced figure: {path.relative_to(SUB)} {dpi[0]:.0f}x{dpi[1]:.0f}"))

    tiffs = sorted({path.with_suffix(".tiff") for path in refs})
    for path in tiffs:
        if not path.exists():
            checks.append(Check("WARN", f"standalone TIFF not found: {path.relative_to(SUB)}"))
            continue
        dpi = image_dpi(path)
        if dpi and min(dpi) >= 600:
            checks.append(Check("PASS", f"600 dpi TIFF available: {path.relative_to(SUB)}"))
        else:
            checks.append(Check("WARN", f"TIFF does not report 600 dpi: {path.relative_to(SUB)} dpi={dpi}"))
        compression = tiff_compression(path)
        if compression in {5, 8}:
            checks.append(Check("PASS", f"TIFF uses lossless compression: {path.relative_to(SUB)} tag={compression}"))
        else:
            checks.append(Check("WARN", f"TIFF is not LZW/ZIP compressed: {path.relative_to(SUB)} tag={compression}"))
    return checks


def main() -> int:
    checks: list[Check] = []
    checks.extend(check_files())
    checks.extend(check_main_markdown())
    checks.extend(check_supplement())
    checks.extend(check_docx())
    checks.extend(check_artifact_freshness())
    checks.extend(check_metadata_and_reference_files())
    checks.extend(check_quantitative_claims())
    checks.extend(check_figures())

    for check in checks:
        print(f"[{check.level}] {check.message}")

    failures = [c for c in checks if c.level == "FAIL"]
    warnings = [c for c in checks if c.level == "WARN"]
    print(f"\nSummary: {len(failures)} failure(s), {len(warnings)} warning(s), {len(checks)} checks.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
