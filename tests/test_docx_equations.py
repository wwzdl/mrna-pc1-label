from __future__ import annotations

import sys
import zipfile
from pathlib import Path

from lxml import etree


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from render_markdown_to_docx import render_markdown  # noqa: E402


def test_inline_and_numbered_equations_are_editable_omml(tmp_path: Path) -> None:
    source = """# Equation test

Inline $T_{\\lambda}$ remains editable.

```equation
number: 1
T_{\\lambda} = \\frac{R_{\\mathrm{final}}^2}{1-R_{\\mathrm{prior}}^2}
```
"""
    output = tmp_path / "equation_test.docx"
    render_markdown(source, output, tmp_path)

    with zipfile.ZipFile(output) as archive:
        xml = archive.read("word/document.xml").decode("utf-8")

    assert xml.count("<m:oMath") >= 2
    assert "<m:f>" in xml
    assert "<m:sSub>" in xml
    assert "<m:sSup>" in xml or "<m:sSubSup>" in xml
    assert "<w:t>(1)</w:t>" in xml
    assert "number: 1" not in xml

    root = etree.fromstring(xml.encode("utf-8"))
    ns = {
        "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
        "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    }
    equation_table = root.xpath("//w:tbl[.//m:oMath]", namespaces=ns)[0]
    grid_widths = [
        int(node.get(f"{{{ns['w']}}}w"))
        for node in equation_table.xpath("./w:tblGrid/w:gridCol", namespaces=ns)
    ]
    cell_widths = [
        int(node.get(f"{{{ns['w']}}}w"))
        for node in equation_table.xpath("./w:tr/w:tc/w:tcPr/w:tcW", namespaces=ns)
    ]
    table_width = int(
        equation_table.xpath("./w:tblPr/w:tblW", namespaces=ns)[0].get(
            f"{{{ns['w']}}}w"
        )
    )

    assert grid_widths == cell_widths == [691, 8482, 691]
    assert table_width == sum(grid_widths)
    assert equation_table.xpath("./w:tblPr/w:tblLayout[@w:type='fixed']", namespaces=ns)
    assert equation_table.xpath("./w:tr/w:trPr/w:cantSplit", namespaces=ns)
