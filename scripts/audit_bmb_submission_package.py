#!/usr/bin/env python3
"""Audit the local BMB submission package for common hard-format issues."""

from __future__ import annotations

import csv
import os
import re
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path

from PIL import Image


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


def markdown_equation_numbers(path: Path) -> list[str]:
    return re.findall(
        r"```equation\s*\nnumber:\s*([^\s]+)\s*\n.+?```",
        read_text(path),
        flags=re.S,
    )


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

    required_universes = ["12,916", "11,107", "12,307", "10,768"]
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

    def require_tokens(label: str, tokens: list[str]) -> None:
        missing = [token for token in tokens if token not in manuscript]
        if missing:
            checks.append(Check("FAIL", f"main-text quantitative claim mismatch for {label}: missing {missing}"))
        else:
            checks.append(Check("PASS", f"main-text quantitative claim matches result table: {label}"))

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

    influence_null = read_tsv(
        ROOT / "results" / "study_influence_sensitivity" / "size_matched_null_summary.tsv"
    )
    for metric in (
        "pc1_stability_pearson",
        "delta_saluki_pearson",
        "delta_ortholog_pearson",
    ):
        row = select_row(influence_null, metric=metric)
        display_precision = 3 if metric == "pc1_stability_pearson" else 4
        require_tokens(
            f"size-matched study-influence null: {metric}",
            [
                f"{float(row['observed_gejman']):.{display_precision}f}",
                f"{float(row['null_median']):.{display_precision}f}",
                f"{float(row['empirical_p_one_sided']):.3f}",
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
