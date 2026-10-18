"""Re-render the 117 cached zw-diagram outputs in Chang--Chen's worksheet.

The source notebook is the Mathematica worksheet distributed by the authors.
This renderer extracts its 117 cached GraphBox objects and independently
rebuilds SVG and a paginated TikZ/PDF source. It verifies extraction counts,
labels, edge endpoints, and the manuscript's label partition; it does not rerun
the authors' Mathematica Algorithm I.
"""
from __future__ import annotations

import json
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "SingSeqCCver202303.nb"
OUT = ROOT / "diagrams_117"
POSITIONS = (
    (1.0, 0.0),
    (0.5, 0.8660254037844386),
    (-0.5, 0.8660254037844386),
    (-1.0, 0.0),
    (-0.5, -0.8660254037844386),
    (0.5, -0.8660254037844386),
)
COLOR_RGB = {
    "0.2, 0.7, 1": (51, 179, 255),
    "1, 0.3, 0.3": (255, 77, 77),
    "0.5, 0.3, 0.9": (128, 77, 230),
}
COLOR_NAMES = {
    "0.2, 0.7, 1": "blue",
    "1, 0.3, 0.3": "red",
    "0.5, 0.3, 0.9": "purple",
}


def graphics_boxes(text: str) -> list[str]:
    """Return balanced GraphicsBox expressions, respecting quoted strings."""
    blocks = []
    cursor = 0
    token = "GraphicsBox["
    while True:
        start = text.find(token, cursor)
        if start < 0:
            return blocks
        i = start + len(token)
        depth, quoted, escaped = 1, False, False
        while i < len(text) and depth:
            char = text[i]
            if quoted:
                if escaped:
                    escaped = False
                elif char == "\\":
                    escaped = True
                elif char == '"':
                    quoted = False
            else:
                if char == '"':
                    quoted = True
                elif char == "[":
                    depth += 1
                elif char == "]":
                    depth -= 1
            i += 1
        if depth:
            raise ValueError(f"Unclosed GraphicsBox at byte {start}")
        blocks.append(text[start:i])
        cursor = i


def color_key(raw: str) -> str:
    return ", ".join(str(float(x.strip())).rstrip("0").rstrip(".")
                     if "." in x else x.strip()
                     for x in raw.split(","))


def vertex_index(x: float, y: float) -> int:
    distances = [(x - px) ** 2 + (y - py) ** 2 for px, py in POSITIONS]
    i = min(range(6), key=distances.__getitem__)
    if distances[i] > 1e-12:
        raise ValueError(f"Unexpected vertex coordinate {(x, y)}")
    return i + 1


def parse_graph(block: str, expected_number: int) -> dict:
    labels = {int(n) for n in re.findall(r'InsetBox\["(\d+)", \{1\.1, -1\}\]', block)}
    if labels != {expected_number}:
        raise ValueError(f"Graphic #{expected_number} has number labels {labels}")

    disk_matches = list(re.finditer(
        r'(RGBColor\[([^\]]+)\]|GrayLevel\[([^\]]+)\]),\s*'
        r'EdgeForm\[GrayLevel\[0\]\],\s*DiskBox\[\{([^}]+)\},\s*([^]]+)\]', block
    ))
    if len(disk_matches) != 6:
        raise ValueError(f"Graphic #{expected_number}: expected 6 node disks, got {len(disk_matches)}")
    vertices = []
    for match in disk_matches:
        if match.group(2) is not None:
            key = color_key(match.group(2))
            circled = True
        else:
            shade = float(match.group(3))
            key = "white" if shade == 1 else f"gray-{shade:g}"
            circled = shade != 1
        xy = [float(v.strip()) for v in match.group(4).split(",")]
        idx = vertex_index(*xy)
        vertices.append({"vertex": idx, "color": key, "circled": circled})
    vertices.sort(key=lambda row: row["vertex"])
    if [v["vertex"] for v in vertices] != list(range(1, 7)):
        raise ValueError(f"Graphic #{expected_number}: duplicate/missing vertex labels")

    edge_source = block[:block.find("DiskBox[")]
    edges = []
    for match in re.finditer(r'RGBColor\[([^\]]+)\],\s*Thickness\[[^]]+\],\s*LineBox\[\{\{([^}]+)\},\s*\{([^}]+)\}\}\]', edge_source):
        key = color_key(match.group(1))
        if key not in COLOR_RGB:
            raise ValueError(f"Unknown edge color {key} in graphic #{expected_number}")
        a = [float(v.strip()) for v in match.group(2).split(",")]
        b = [float(v.strip()) for v in match.group(3).split(",")]
        i, j = sorted((vertex_index(*a), vertex_index(*b)))
        edges.append({"vertices": [i, j], "color": key, "palette": COLOR_NAMES[key]})
    if len(edges) != len({tuple(e["vertices"]) for e in edges}):
        raise ValueError(f"Graphic #{expected_number}: duplicate edge")
    return {"number": expected_number, "vertices": vertices, "edges": sorted(edges, key=lambda e: e["vertices"])}


def rgb_css(key: str) -> str:
    r, g, b = COLOR_RGB[key]
    return f"#{r:02x}{g:02x}{b:02x}"


def render_svg(graph: dict) -> str:
    width = height = 480
    cx, cy, scale = width / 2, height / 2 + 10, 150
    coords = {i + 1: (cx + x * scale, cy - y * scale) for i, (x, y) in enumerate(POSITIONS)}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{cx}" y="32" text-anchor="middle" font-family="Arial,sans-serif" font-size="22" font-weight="bold">No. {graph["number"]}</text>']
    for edge in graph["edges"]:
        (x1, y1), (x2, y2) = (coords[v] for v in edge["vertices"])
        parts.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="{rgb_css(edge["color"])}" stroke-width="5" stroke-linecap="round"/>')
    for vertex in graph["vertices"]:
        x, y = coords[vertex["vertex"]]
        fill = "white" if vertex["color"] == "white" else rgb_css(vertex["color"])
        parts.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="37" fill="{fill}" stroke="#222" stroke-width="2"/>')
        parts.append(f'<text x="{x:.3f}" y="{y + 7:.3f}" text-anchor="middle" font-family="Arial,sans-serif" font-size="21" fill="#111">{vertex["vertex"]}</text>')
    parts += ['<text x="240" y="460" text-anchor="middle" font-family="Arial,sans-serif" font-size="12" fill="#555">vertex and edge colors retained from the source worksheet</text>', '</svg>']
    return "\n".join(parts) + "\n"


def tex_rgb(key: str) -> str:
    r, g, b = COLOR_RGB[key]
    return f"{{rgb,255:red,{r};green,{g};blue,{b}}}"


def render_tikz(graph: dict) -> str:
    coords = {i + 1: (x * 1.45, y * 1.45) for i, (x, y) in enumerate(POSITIONS)}
    out = [r"\begin{minipage}[c][0.39\textheight][c]{0.315\textwidth}", r"\centering", f"{{\\large\\bfseries No. {graph['number']}}}", r"\vspace{2mm}", r"\begin{tikzpicture}[line cap=round,line join=round]"]
    for edge in graph["edges"]:
        i, j = edge["vertices"]
        out.append(f"\\draw[line width=1.05pt,color={tex_rgb(edge['color'])}] ({coords[i][0]:.6f},{coords[i][1]:.6f}) -- ({coords[j][0]:.6f},{coords[j][1]:.6f});")
    for vertex in graph["vertices"]:
        i = vertex["vertex"]
        x, y = coords[i]
        fill = "white" if vertex["color"] == "white" else tex_rgb(vertex["color"])
        out.append(f"\\filldraw[draw=black,fill={fill}] ({x:.6f},{y:.6f}) circle (0.31cm) node[font=\\small,text=black] {{{i}}};")
    out += [r"\end{tikzpicture}", r"\end{minipage}"]
    return "\n".join(out)


def write_catalog_tex(graphs: list[dict]) -> None:
    tex = [r"\documentclass[10pt]{article}", r"\usepackage[margin=1.1cm]{geometry}", r"\usepackage{tikz}", r"\usepackage{xcolor}", r"\pagestyle{plain}", r"\setlength{\parindent}{0pt}", r"\begin{document}"]
    for i, graph in enumerate(graphs):
        tex.append(render_tikz(graph))
        if i % 6 == 5 and i != len(graphs) - 1:
            tex.append(r"\newpage")
        elif i % 3 == 2 and i != len(graphs) - 1:
            tex.append(r"\par\vspace{2mm}")
    tex.append(r"\end{document}")
    (OUT / "diagram_catalog_117.tex").write_text("\n".join(tex) + "\n", encoding="utf-8")


def write_gallery(graphs: list[dict]) -> None:
    cards = []
    for graph in graphs:
        number = graph["number"]
        name = f"diagram_{number:03d}.svg"
        cards.append(
            f'<article><h2>No. {number}</h2><a href="{name}"><img src="{name}" alt="No. {number} zw-diagram"></a>'
            f'<p><a href="{name}" download>Open / save SVG</a></p></article>'
        )
    html = """<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>N=6: 117 zw-diagrams</title>
<style>body{font-family:system-ui,sans-serif;margin:2rem;color:#222}header{max-width:70rem;margin:auto}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:1rem;max-width:90rem;margin:2rem auto}article{border:1px solid #ddd;border-radius:8px;padding:.8rem;text-align:center}h2{font-size:1rem;margin:.2rem}img{width:200px;height:200px;object-fit:contain}p{font-size:.85rem;margin:.3rem}small{color:#555}</style>
<header><h1>117 个 N=6 zw 图</h1><p>按作者 Mathematica 工作表中缓存的图形输出重绘。该图册没有重新运行 Algorithm I。</p><small>来源：SingSeqCCver202303.nb；机器可读数据见 diagram_data_117.json。</small></header>
<main>""" + "\n".join(cards) + "</main></html>\n"
    (OUT / "index.html").write_text(html, encoding="utf-8")


def main() -> None:
    if not SOURCE.exists():
        raise SystemExit(f"Missing official source worksheet: {SOURCE}")
    boxes = graphics_boxes(SOURCE.read_text(encoding="utf-8", errors="replace"))
    if len(boxes) != 117:
        raise SystemExit(f"Expected exactly 117 cached diagram graphics; found {len(boxes)}")
    graphs = [parse_graph(block, number) for number, block in enumerate(boxes, 1)]
    OUT.mkdir(exist_ok=True)
    for graph in graphs:
        path = OUT / f"diagram_{graph['number']:03d}.svg"
        path.write_text(render_svg(graph), encoding="utf-8")
    (OUT / "diagram_data_117.json").write_text(json.dumps({
        "source": "SingSeqCCver202303.nb (author-distributed worksheet)",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_graphics_count": len(boxes),
        "label_check": "all cached graphics were ordered and numbered 1 through 117",
        "scope": "re-rendering of cached author outputs, not an independent rerun of Algorithm I",
        "color_key": {"edge colors": COLOR_NAMES, "vertices": "white means no vertex circle; nonwhite source fill means a circled vertex"},
        "diagrams": graphs,
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_catalog_tex(graphs)
    write_gallery(graphs)
    assert len(list(OUT.glob("diagram_*.svg"))) == 117
    assert {g["number"] for g in graphs} == set(range(1, 118))
    print(f"Extracted and rendered {len(graphs)} numbered diagrams.")
    print(f"Individual SVGs and machine-readable graph data: {OUT}")
    print(f"Paginated PDF source: {OUT / 'diagram_catalog_117.tex'}")
    print("Scope: exact cached output redraw; classification algorithm not independently rerun.")


if __name__ == "__main__":
    main()
