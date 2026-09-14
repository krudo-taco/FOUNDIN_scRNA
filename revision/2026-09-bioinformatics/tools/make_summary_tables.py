#!/usr/bin/env python3
"""Create compact public GSE243639 summary table and SVG figure."""
import argparse
import csv
from pathlib import Path


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def write_svg(rows, out_path):
    tested = [r for r in rows if r.get("p_omnibus") not in ("", "nan")]
    width, height = 900, 320
    margin = 70
    max_value = max(float(r["p_omnibus"]) for r in tested) if tested else 1.0
    bar_w = (width - 2 * margin) / max(len(tested), 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="70" y="35" font-family="Arial" font-size="18">GSE243639 pilot omnibus P values</text>',
    ]
    for i, row in enumerate(tested):
        p = float(row["p_omnibus"])
        h = 220 * p / max_value
        x = margin + i * bar_w + 4
        y = height - margin - h
        label = f'{row["context"]}/{row["method"]}/{row["test"]}'
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w - 8:.1f}" height="{h:.1f}" fill="#4f7cac"/>')
        parts.append(f'<text x="{x:.1f}" y="{height - 45}" transform="rotate(45 {x:.1f},{height - 45})" font-family="Arial" font-size="10">{label}</text>')
        parts.append(f'<text x="{x:.1f}" y="{y - 4:.1f}" font-family="Arial" font-size="10">{p:.3g}</text>')
    parts.append(f'<line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="black"/>')
    parts.append("</svg>")
    out_path.write_text("\n".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-root", required=True)
    args = parser.parse_args()
    work = Path(args.work_root)
    results = work / "results"
    summary = work / "summary"
    summary.mkdir(parents=True, exist_ok=True)
    panel_rows = read_csv(results / "gse243639_pilot_panel_tests.csv")
    out_csv = summary / "gse243639_panel_tests_summary.csv"
    with open(out_csv, "w", newline="") as handle:
        fields = ["context", "method", "model", "test", "genes_tested", "p_omnibus", "p_holm_all_pilot_tests"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in panel_rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    write_svg(panel_rows, summary / "gse243639_panel_tests.svg")
    print(f"Wrote {out_csv}")
    print(f"Wrote {summary / 'gse243639_panel_tests.svg'}")


if __name__ == "__main__":
    main()

