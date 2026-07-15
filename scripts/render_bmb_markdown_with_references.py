#!/usr/bin/env python3
"""Render markdown with optional author-year citations from a BibTeX file.

Supported workflow:
1. Markdown contains simple pandoc-style citation clusters such as
   `[@agarwal2022; @chen2016]`.
2. The script replaces those clusters with author-year citations.
3. The script appends or replaces a `# References` / `## References` section.

If a source file contains no pandoc-style citation clusters, the script copies
the file unchanged. This keeps the helper safe for the current BMB workflow,
where some English drafts still use manually written reference sections.
"""

from __future__ import annotations

import argparse
import re
import shutil
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BIB = ROOT / "manuscript" / "bmb_submission" / "references" / "bmb_references.bib"
DEFAULT_JOBS = [
    (
        ROOT / "manuscript" / "bmb_submission" / "bmb_main_manuscript_en.md",
        ROOT / "manuscript" / "bmb_submission" / "bmb_main_manuscript_en.md",
    ),
    (
        ROOT / "manuscript" / "bmb_submission" / "bmb_supplementary_material_en.md",
        ROOT / "manuscript" / "bmb_submission" / "bmb_supplementary_material_en.md",
    ),
]


def parse_bibtex_entries(path: Path) -> dict[str, dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    entries: dict[str, dict[str, str]] = {}
    pattern = re.compile(r"@(\w+)\{([^,]+),")
    pos = 0
    while True:
        match = pattern.search(text, pos)
        if not match:
            break
        entry_type, key = match.group(1).lower(), match.group(2).strip()
        brace_depth = 0
        end = match.end()
        for idx in range(match.start(), len(text)):
            char = text[idx]
            if char == "{":
                brace_depth += 1
            elif char == "}":
                brace_depth -= 1
                if brace_depth == 0:
                    end = idx
                    break
        block = text[match.end() : end]
        fields: dict[str, str] = {"entry_type": entry_type, "key": key}
        field_pattern = re.compile(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", re.S)
        for field_match in field_pattern.finditer(block):
            field_name = field_match.group(1).strip().lower()
            field_value = " ".join(field_match.group(2).split())
            fields[field_name] = field_value
        entries[key] = fields
        pos = end + 1
    return entries


def split_authors(author_field: str) -> list[str]:
    return [part.strip() for part in author_field.split(" and ") if part.strip()]


def parse_person(name: str) -> tuple[str, str]:
    if "," in name:
        family, given = [part.strip() for part in name.split(",", 1)]
        return family, given
    parts = name.split()
    if not parts:
        return "", ""
    if len(parts) == 1:
        return parts[0], ""
    return parts[-1], " ".join(parts[:-1])


def initials(given: str) -> str:
    letters = []
    for token in re.split(r"[\s.-]+", given):
        token = token.strip()
        if token:
            letters.append(token[0].upper())
    return "".join(letters)


def citation_author_token(entry: dict[str, str]) -> str:
    authors = split_authors(entry.get("author", ""))
    if not authors:
        return "Unknown"
    families = [parse_person(name)[0] for name in authors]
    if len(families) == 1:
        return families[0]
    if len(families) == 2:
        return f"{families[0]} and {families[1]}"
    return f"{families[0]} et al."


def reference_author_list(entry: dict[str, str]) -> str:
    authors = split_authors(entry.get("author", ""))
    if not authors:
        return "Unknown"
    formatted = []
    for author in authors:
        family, given = parse_person(author)
        if author.lower() == "others":
            formatted.append("et al")
            continue
        if family:
            init = initials(given)
            formatted.append(f"{family} {init}".strip())
        else:
            formatted.append(author)
    return ", ".join(formatted)


def extract_ordered_keys(text: str) -> list[str]:
    keys: list[str] = []
    seen: set[str] = set()
    for match in re.finditer(r"\[(?:[^\]]*?@[^\]]*?)\]", text):
        inside = match.group(0)[1:-1]
        for part in inside.split(";"):
            part = part.strip()
            if "@" not in part:
                continue
            key = part.split("@", 1)[1].strip()
            key = key.split(",", 1)[0].strip()
            if key and key not in seen:
                seen.add(key)
                keys.append(key)
    return keys


def build_year_suffix_map(keys: list[str], entries: dict[str, dict[str, str]]) -> dict[str, str]:
    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for key in keys:
        entry = entries.get(key, {})
        authors = split_authors(entry.get("author", ""))
        families = tuple(parse_person(name)[0] for name in authors)
        year = entry.get("year", "n.d.")
        grouped[(str(families), year)].append(key)
    suffix_map: dict[str, str] = {}
    for _, dup_keys in grouped.items():
        if len(dup_keys) == 1:
            suffix_map[dup_keys[0]] = ""
            continue
        for idx, key in enumerate(sorted(dup_keys, key=lambda k: entries.get(k, {}).get("title", ""))):
            suffix_map[key] = chr(ord("a") + idx)
    return suffix_map


def citation_text(key: str, entry: dict[str, str], suffix_map: dict[str, str]) -> str:
    author_part = citation_author_token(entry)
    year = entry.get("year", "n.d.")
    suffix = suffix_map.get(key, "")
    return f"{author_part} {year}{suffix}"


def replace_citations(text: str, entries: dict[str, dict[str, str]], suffix_map: dict[str, str]) -> str:
    def repl(match: re.Match[str]) -> str:
        inside = match.group(0)[1:-1]
        cite_items: list[str] = []
        seen: set[str] = set()
        for part in inside.split(";"):
            part = part.strip()
            if "@" not in part:
                continue
            key = part.split("@", 1)[1].strip()
            key = key.split(",", 1)[0].strip()
            if key in entries and key not in seen:
                seen.add(key)
                cite_items.append(citation_text(key, entries[key], suffix_map))
        if not cite_items:
            return match.group(0)
        return "(" + "; ".join(cite_items) + ")"

    return re.sub(r"\[(?:[^\]]*?@[^\]]*?)\]", repl, text)


def reference_sort_key(key: str, entry: dict[str, str], suffix_map: dict[str, str]) -> tuple[str, str, str, str]:
    authors = split_authors(entry.get("author", ""))
    first_family = parse_person(authors[0])[0] if authors else "Unknown"
    coauthor_token = parse_person(authors[1])[0] if len(authors) > 1 else ""
    year = entry.get("year", "n.d.") + suffix_map.get(key, "")
    title = entry.get("title", "")
    return (first_family.lower(), coauthor_token.lower(), year, title.lower())


def build_reference_line(key: str, entry: dict[str, str], suffix_map: dict[str, str]) -> str:
    authors = reference_author_list(entry)
    title = entry.get("title", "Untitled").rstrip(".")
    venue = entry.get("journal") or entry.get("booktitle") or entry.get("publisher") or "Unknown venue"
    year = entry.get("year", "n.d.") + suffix_map.get(key, "")
    volume = entry.get("volume", "")
    number = entry.get("number", "")
    pages = entry.get("pages", "").replace("--", "-")
    doi = entry.get("doi", "")

    meta_parts = []
    if volume:
        meta_parts.append(f"{volume}({number})" if number else volume)
    elif number:
        meta_parts.append(f"({number})")
    if pages:
        meta_parts.append(pages)

    venue_part = venue
    if meta_parts:
        venue_part += " " + ":".join(meta_parts)

    line = f"{authors} ({year}) {title}. {venue_part}."
    if doi:
        doi = doi.removeprefix("https://doi.org/")
        line += f" https://doi.org/{doi}"
    return line


def inject_reference_section(text: str, reference_lines: list[str]) -> str:
    rendered_refs = "\n".join(reference_lines)
    heading_patterns = [
        r"(?m)^# References\s*$",
        r"(?m)^## References\s*$",
    ]
    for pattern in heading_patterns:
        match = re.search(pattern, text)
        if match:
            return text[: match.end()] + "\n\n" + rendered_refs + "\n"
    return text.rstrip() + "\n\n# References\n\n" + rendered_refs + "\n"


def render_text(text: str, entries: dict[str, dict[str, str]]) -> tuple[str, bool]:
    ordered_keys = extract_ordered_keys(text)
    if not ordered_keys:
        return text, False
    suffix_map = build_year_suffix_map(ordered_keys, entries)
    rendered_text = replace_citations(text, entries, suffix_map)
    sorted_keys = sorted(
        [key for key in ordered_keys if key in entries],
        key=lambda key: reference_sort_key(key, entries[key], suffix_map),
    )
    reference_lines = [build_reference_line(key, entries[key], suffix_map) for key in sorted_keys]
    rendered_text = inject_reference_section(rendered_text, reference_lines)
    return rendered_text, True


def render_file(source: Path, dest: Path, entries: dict[str, dict[str, str]]) -> str:
    text = source.read_text(encoding="utf-8")
    rendered_text, changed = render_text(text, entries)
    if not changed:
        if source.resolve() != dest.resolve():
            shutil.copyfile(source, dest)
        return "copied unchanged (no pandoc-style citation clusters found)"
    dest.write_text(rendered_text, encoding="utf-8")
    return "rendered author-year citations"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", nargs="?", type=Path, help="Source markdown file.")
    parser.add_argument("dest", nargs="?", type=Path, help="Output markdown file.")
    parser.add_argument(
        "--bib",
        type=Path,
        default=DEFAULT_BIB,
        help=f"BibTeX file to use (default: {DEFAULT_BIB}).",
    )
    parser.add_argument(
        "--all-bmb-en",
        action="store_true",
        help="Process the current English BMB main text and supplement in place.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.all_bmb_en:
        jobs = DEFAULT_JOBS
    else:
        if args.source is None or args.dest is None:
            raise SystemExit(
                "Provide SOURCE and DEST, or use --all-bmb-en to process the current "
                "English BMB main/supplement files in place."
            )
        jobs = [(args.source, args.dest)]

    if not args.bib.exists():
        raise SystemExit(f"BibTeX file not found: {args.bib}")

    entries = parse_bibtex_entries(args.bib)
    for source, dest in jobs:
        if not source.exists():
            raise SystemExit(f"Source markdown not found: {source}")
        status = render_file(source, dest, entries)
        try:
            rel = dest.relative_to(ROOT)
        except ValueError:
            rel = dest
        print(f"{rel}: {status}")


if __name__ == "__main__":
    main()
