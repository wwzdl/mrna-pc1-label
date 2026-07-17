from __future__ import annotations

import argparse
from pathlib import Path

import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode


RIS_TYPES = {
    "article": "JOUR",
    "inproceedings": "CPAPER",
    "conference": "CPAPER",
    "book": "BOOK",
    "inbook": "CHAP",
    "incollection": "CHAP",
    "phdthesis": "THES",
    "mastersthesis": "THES",
    "techreport": "RPRT",
}


def _page_fields(value: str) -> list[tuple[str, str]]:
    normalized = value.replace("---", "--")
    if "--" not in normalized:
        return [("SP", normalized)]
    start, end = normalized.split("--", maxsplit=1)
    return [("SP", start), ("EP", end)]


def _ris_lines(entry: dict[str, str]) -> list[str]:
    entry_type = entry.get("ENTRYTYPE", "misc").lower()
    lines = [f"TY  - {RIS_TYPES.get(entry_type, 'GEN')}" ]

    for author in entry.get("author", "").split(" and "):
        author = author.strip()
        if author:
            lines.append(f"AU  - {author}")

    field_map = (
        ("title", "TI"),
        ("journal", "JO"),
        ("booktitle", "T2"),
        ("year", "PY"),
        ("volume", "VL"),
        ("number", "IS"),
        ("publisher", "PB"),
        ("doi", "DO"),
        ("url", "UR"),
    )
    for bib_field, ris_field in field_map:
        value = entry.get(bib_field, "").strip()
        if value:
            lines.append(f"{ris_field}  - {value}")

    pages = entry.get("pages", "").strip()
    if pages:
        lines.extend(f"{field}  - {value}" for field, value in _page_fields(pages))

    lines.extend((f"ID  - {entry['ID']}", "ER  -", ""))
    return lines


def export_bib_to_ris(input_path: Path, output_path: Path) -> int:
    parser = BibTexParser(common_strings=True)
    parser.customization = convert_to_unicode
    with input_path.open(encoding="utf-8") as handle:
        database = bibtexparser.load(handle, parser=parser)

    ids = [entry["ID"] for entry in database.entries]
    if len(ids) != len(set(ids)):
        duplicates = sorted({key for key in ids if ids.count(key) > 1})
        raise ValueError(f"duplicate BibTeX keys: {', '.join(duplicates)}")

    entries = sorted(
        database.entries,
        key=lambda entry: (entry.get("author", "").casefold(), entry.get("year", ""), entry["ID"]),
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "\n".join(line for entry in entries for line in _ris_lines(entry)),
        encoding="utf-8",
    )
    return len(entries)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export the project BibTeX library as EndNote-compatible RIS.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    count = export_bib_to_ris(args.input, args.output)
    print(f"[done] exported {count} unique references to {args.output}")


if __name__ == "__main__":
    main()
