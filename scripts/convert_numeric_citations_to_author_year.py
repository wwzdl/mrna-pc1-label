#!/usr/bin/env python3
"""Convert simple numbered manuscript citations to author-year style.

This helper is intentionally conservative. It uses the numbered reference list
already present in a markdown manuscript, maps each numbered item to the local
BibTeX library by DOI or title, replaces in-text numeric clusters such as
`[1, 2, 3]` with author-year clusters, and rewrites the reference list sorted by
first author. It is designed for the current BMB submission markdown files.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from render_bmb_markdown_with_references import (  # noqa: E402
    build_reference_line,
    build_year_suffix_map,
    citation_text,
    parse_bibtex_entries,
    reference_sort_key,
)


ROOT = SCRIPT_DIR.parents[0]
DEFAULT_BIB = ROOT / "manuscript" / "bmb_submission" / "references" / "bmb_references.bib"


def norm_doi(value: str) -> str:
    return value.lower().strip().rstrip(".").removeprefix("https://doi.org/")


def norm_title(value: str) -> str:
    value = re.sub(r"[*_`{}\\]", "", value.lower())
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return " ".join(value.split())


def split_reference_section(text: str) -> tuple[str, str, str]:
    matches = list(re.finditer(r"(?m)^(##\s+.*(?:References|参考文献).*)$", text, flags=re.I))
    if not matches:
        raise ValueError("No markdown reference heading found")
    match = matches[-1]
    return text[: match.start()], match.group(1), text[match.end() :]


def parse_numbered_references(ref_text: str) -> dict[int, str]:
    refs: dict[int, str] = {}
    current_num: int | None = None
    current_lines: list[str] = []
    for line in ref_text.splitlines():
        match = re.match(r"^\s*(\d+)\.\s+(.*)", line)
        if match:
            if current_num is not None:
                refs[current_num] = " ".join(current_lines).strip()
            current_num = int(match.group(1))
            current_lines = [match.group(2).strip()]
        elif current_num is not None and line.strip():
            current_lines.append(line.strip())
    if current_num is not None:
        refs[current_num] = " ".join(current_lines).strip()
    return refs


def raw_ref_doi(raw: str) -> str:
    match = re.search(r"doi:\s*(\S+)", raw, flags=re.I)
    return norm_doi(match.group(1)) if match else ""


def raw_ref_title_guess(raw: str) -> str:
    # Most current references follow: Authors. Title. *Venue*. ...
    match = re.search(r"\.\s+([^*.]+?)\.\s+\*", raw)
    if match:
        return norm_title(match.group(1))
    return ""


def first_family_from_raw(raw: str) -> str:
    raw = raw.strip()
    if "," in raw:
        return raw.split(",", 1)[0].strip()
    return raw.split()[0] if raw.split() else "Unknown"


def year_from_raw(raw: str) -> str:
    before_doi = raw.split("doi:", 1)[0]
    years = re.findall(r"\b(19\d{2}|20\d{2})\b", before_doi)
    return years[-1] if years else "n.d."


def build_maps(numbered_refs: dict[int, str], entries: dict[str, dict[str, str]]):
    by_doi = {
        norm_doi(entry.get("doi", "")): key
        for key, entry in entries.items()
        if entry.get("doi")
    }
    by_title = {
        norm_title(entry.get("title", "")): key
        for key, entry in entries.items()
        if entry.get("title")
    }
    number_to_key: dict[int, str] = {}
    fallback: dict[int, tuple[str, str]] = {}
    for number, raw in numbered_refs.items():
        key = ""
        doi = raw_ref_doi(raw)
        if doi:
            key = by_doi.get(doi, "")
        if not key:
            title = raw_ref_title_guess(raw)
            key = by_title.get(title, "")
        if key:
            number_to_key[number] = key
        else:
            fallback[number] = (first_family_from_raw(raw), year_from_raw(raw))
    return number_to_key, fallback


def expand_number_token(token: str) -> list[int]:
    token = token.strip()
    if "-" in token:
        left, right = [part.strip() for part in token.split("-", 1)]
        if left.isdigit() and right.isdigit():
            start, end = int(left), int(right)
            if start <= end:
                return list(range(start, end + 1))
    return [int(token)] if token.isdigit() else []


def replace_numeric_citations(
    body: str,
    number_to_key: dict[int, str],
    fallback: dict[int, tuple[str, str]],
    entries: dict[str, dict[str, str]],
    suffix_map: dict[str, str],
) -> str:
    def repl(match: re.Match[str]) -> str:
        numbers: list[int] = []
        for token in re.split(r",|;", match.group(1)):
            numbers.extend(expand_number_token(token))
        if not numbers:
            return match.group(0)
        cite_items: list[str] = []
        for number in numbers:
            key = number_to_key.get(number)
            if key and key in entries:
                cite_items.append(citation_text(key, entries[key], suffix_map))
            elif number in fallback:
                family, year = fallback[number]
                cite_items.append(f"{family} {year}")
        if not cite_items:
            return match.group(0)
        seen: set[str] = set()
        unique_items = []
        for item in cite_items:
            if item not in seen:
                seen.add(item)
                unique_items.append(item)
        return "(" + "; ".join(unique_items) + ")"

    return re.sub(r"\[([0-9,\s;\-]+)\]", repl, body)


def render_reference_section(
    heading: str,
    numbered_refs: dict[int, str],
    number_to_key: dict[int, str],
    fallback: dict[int, tuple[str, str]],
    entries: dict[str, dict[str, str]],
    suffix_map: dict[str, str],
) -> str:
    mapped_numbers = [n for n in numbered_refs if n in number_to_key]
    mapped_numbers = sorted(
        mapped_numbers,
        key=lambda n: reference_sort_key(number_to_key[n], entries[number_to_key[n]], suffix_map),
    )
    fallback_numbers = sorted(
        [n for n in numbered_refs if n in fallback],
        key=lambda n: (fallback[n][0].lower(), fallback[n][1], numbered_refs[n].lower()),
    )
    lines = [heading, ""]
    for number in mapped_numbers:
        key = number_to_key[number]
        lines.append(build_reference_line(key, entries[key], suffix_map))
    for number in fallback_numbers:
        family, year = fallback[number]
        raw = numbered_refs[number]
        lines.append(f"{family} ({year}) {raw}")
    return "\n".join(lines).rstrip() + "\n"


def convert_file(path: Path, bib: Path) -> tuple[int, int]:
    text = path.read_text(encoding="utf-8")
    body, heading, ref_text = split_reference_section(text)
    numbered_refs = parse_numbered_references(ref_text)
    entries = parse_bibtex_entries(bib)
    number_to_key, fallback = build_maps(numbered_refs, entries)
    ordered_keys = [number_to_key[n] for n in sorted(number_to_key)]
    suffix_map = build_year_suffix_map(ordered_keys, entries)
    converted_body = replace_numeric_citations(body, number_to_key, fallback, entries, suffix_map)
    converted_refs = render_reference_section(
        heading,
        numbered_refs,
        number_to_key,
        fallback,
        entries,
        suffix_map,
    )
    path.write_text(converted_body.rstrip() + "\n\n" + converted_refs, encoding="utf-8")
    return len(numbered_refs), len(number_to_key)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("markdown", type=Path)
    parser.add_argument("--bib", type=Path, default=DEFAULT_BIB)
    args = parser.parse_args()
    total, mapped = convert_file(args.markdown, args.bib)
    print(f"Converted {args.markdown}: {mapped}/{total} references mapped to BibTeX")


if __name__ == "__main__":
    main()
