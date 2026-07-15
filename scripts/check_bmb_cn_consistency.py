#!/usr/bin/env python3
"""Validate the active BMB manuscript chains.

Checks:
1. Main and supplementary markdown files reference only existing images.
2. Figure captions are contiguous in the selected language chain.
3. Referenced main vs supplementary images are not byte-identical duplicates.

Use `--lang cn`, `--lang en`, or `--lang both` to select the manuscript chain.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = ROOT / "manuscript" / "bmb_submission"

IMAGE_RE = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")

LANG_CONFIG = {
    "cn": {
        "main_md": MANUSCRIPT_DIR / "bmb_main_manuscript_cn.md",
        "supp_md": MANUSCRIPT_DIR / "bmb_supplementary_material_cn.md",
        "main_caption_re": re.compile(r"\*\*图\s+(\d+)\*\*"),
        "supp_caption_re": re.compile(r"\*\*补充图\s+S(\d+)\*\*"),
        "main_label": "图",
        "supp_label": "补充图 S",
    },
    "en": {
        "main_md": MANUSCRIPT_DIR / "bmb_main_manuscript_en.md",
        "supp_md": MANUSCRIPT_DIR / "bmb_supplementary_material_en.md",
        "main_caption_re": re.compile(r"\*\*Fig\.\s+(\d+)\*\*", re.I),
        "supp_caption_re": re.compile(r"\*\*Supplementary Fig\.\s+S(\d+)\*\*", re.I),
        "main_label": "Fig.",
        "supp_label": "Supplementary Fig. S",
    },
}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def extract_images(md_path: Path) -> list[Path]:
    text = read_text(md_path)
    return [(md_path.parent / rel).resolve() for rel in IMAGE_RE.findall(text)]


def extract_caption_numbers(md_path: Path, pattern: re.Pattern[str]) -> list[int]:
    return [int(x) for x in pattern.findall(read_text(md_path))]


def ensure_exists(paths: list[Path], label: str) -> list[str]:
    errors: list[str] = []
    for path in paths:
        if not path.exists():
            errors.append(f"[missing] {label}: {path}")
    return errors


def ensure_contiguous(values: list[int], prefix: str) -> list[str]:
    if not values:
        return [f"[numbering] {prefix}: no captions found"]
    expected = list(range(1, len(values) + 1))
    if values != expected:
        return [f"[numbering] {prefix}: expected {expected}, found {values}"]
    return []


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def find_cross_duplicates(main_images: list[Path], supp_images: list[Path]) -> list[str]:
    errors: list[str] = []
    supp_by_hash: dict[str, Path] = {md5(path): path for path in supp_images}
    for main_path in main_images:
        digest = md5(main_path)
        supp_path = supp_by_hash.get(digest)
        if supp_path is not None:
            errors.append(
                f"[duplicate] main vs supplementary: {main_path.name} == {supp_path.name}"
            )
    return errors


def summarize(
    lang: str,
    main_images: list[Path],
    supp_images: list[Path],
    main_caps: list[int],
    supp_caps: list[int],
    main_label: str,
    supp_label: str,
) -> None:
    print(f"[ok] main image refs: {len(main_images)}")
    print(f"[ok] supplementary image refs: {len(supp_images)}")
    print(f"[ok] {lang} main figure captions: {main_label} 1-{len(main_caps)}")
    supp_range = f"{supp_label}1-S{len(supp_caps)}"
    print(f"[ok] {lang} supplementary figure captions: {supp_range}")


def validate_lang(lang: str) -> list[str]:
    cfg = LANG_CONFIG[lang]
    main_md = cfg["main_md"]
    supp_md = cfg["supp_md"]

    errors: list[str] = []
    for required in (main_md, supp_md):
        if not required.exists():
            errors.append(f"[missing] {lang}: required markdown not found: {required}")
            return errors

    main_images = extract_images(main_md)
    supp_images = extract_images(supp_md)
    main_caps = extract_caption_numbers(main_md, cfg["main_caption_re"])
    supp_caps = extract_caption_numbers(supp_md, cfg["supp_caption_re"])

    errors.extend(ensure_exists(main_images, "main"))
    errors.extend(ensure_exists(supp_images, "supplement"))
    errors.extend(ensure_contiguous(main_caps, f"{lang} main figures"))
    errors.extend(ensure_contiguous(supp_caps, f"{lang} supplementary figures"))
    errors.extend(find_cross_duplicates(main_images, supp_images))

    if not errors:
        summarize(
            lang,
            main_images,
            supp_images,
            main_caps,
            supp_caps,
            cfg["main_label"],
            cfg["supp_label"],
        )
        print(
            "[ok] no byte-identical main/supplement duplicate images are referenced "
            f"in the {lang} manuscript chain"
        )
    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lang",
        choices=("cn", "en", "both"),
        default="cn",
        help="Which manuscript chain to validate (default: cn).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    langs = ("cn", "en") if args.lang == "both" else (args.lang,)
    errors: list[str] = []
    for lang in langs:
        errors.extend(validate_lang(lang))

    if errors:
        for err in errors:
            print(err, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
