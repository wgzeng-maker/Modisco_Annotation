#!/usr/bin/env python
r"""
upgrade_modisco_report.py
=========================
Upgrade a TF-MoDISco HTML report using information from the matching
modisco-lite HDF5 file, and optionally write a standalone PDF.

Why PDF?
--------
The original MoDISco HTML report often references motif/logo images by relative
paths. If the HTML is moved without its image folder, images disappear. The PDF
writer resolves those image paths while the report is still next to its image
folder, then embeds the images into the PDF file.

Examples
--------
Create a standalone PDF:

    python upgrade_modisco_report.py \
      --report profile_report/motifs.html \
      --h5 GC_modisco_profile_v2.h5 \
      --meme JASPAR2024_CORE_non-redundant_pfms_meme.txt \
      --out-pdf motifs_upgraded.pdf

Also write upgraded HTML:

    python upgrade_modisco_report.py \
      --report profile_report/motifs.html \
      --h5 GC_modisco_profile_v2.h5 \
      --meme JASPAR2024_CORE_non-redundant_pfms_meme.txt \
      --out-html motifs_upgraded.html \
      --embed-html-images

Requirements
------------
Core: beautifulsoup4, h5py, hdf5plugin, numpy, matplotlib, logomaker
Optional: none
"""

import argparse
import base64
import csv
import io
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse


BASES = ["A", "C", "G", "T"]
ADDED_HEADERS = [
    "PFM forward",
    "PFM reverse",
    "matched TFs",
    "#subpatterns",
    "functional annotation",
]


def _optional_hdf5plugin():
    try:
        import hdf5plugin  # noqa: F401  -- must be imported before h5py when needed
    except ModuleNotFoundError:
        pass


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--report", required=True,
                        help="Original MoDISco report HTML, e.g. profile_report/motifs.html")
    parser.add_argument("--h5", required=True,
                        help="modisco-lite output .h5")
    parser.add_argument("--meme", required=True,
                        help="JASPAR MEME file used for the original report")
    parser.add_argument("--out", default=None,
                        help="Backward-compatible output path. .pdf writes PDF; other suffix writes HTML.")
    parser.add_argument("--out-html", default=None,
                        help="Optional upgraded HTML output path")
    parser.add_argument("--out-pdf", default=None,
                        help="Standalone PDF output path. Default if no --out/--out-html is provided.")
    parser.add_argument("--embed-html-images", action="store_true",
                        help="Convert relative HTML <img> paths to base64 data URIs in --out-html")
    parser.add_argument("--trim-frac", type=float, default=0.30,
                        help="Trim logos to columns whose CWM signal >= this fraction of max")
    parser.add_argument("--notes", default=None,
                        help="Optional CSV with columns pattern,annotation")
    parser.add_argument("--max-original-images", type=int, default=4,
                        help="Max original report images to include per PDF pattern page")
    parser.add_argument("--max-pages", type=int, default=None,
                        help="Debug/testing limit for number of PDF pages")
    return parser.parse_args(argv)


def parse_meme_names(path):
    """Map JASPAR motif ID -> TF name from MEME header lines."""
    id2name = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("MOTIF"):
                parts = line.split()
                if len(parts) >= 2:
                    id2name[parts[1]] = " ".join(parts[2:]) if len(parts) > 2 else parts[1]
    return id2name


def load_notes(path):
    """Load optional pattern -> annotation mapping."""
    if not path:
        return {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = {name.lower(): name for name in (reader.fieldnames or [])}
        if "pattern" not in fields or "annotation" not in fields:
            raise SystemExit("--notes CSV must contain columns named 'pattern' and 'annotation'")
        return {
            row[fields["pattern"]].strip(): row[fields["annotation"]]
            for row in reader
        }


def normalize_pattern_key(key):
    """Return (group, pattern_name) from common report key formats."""
    key = key.strip()
    if "." in key:
        group, name = key.split(".", 1)
        return group, name
    if "/" in key:
        group, name = key.split("/", 1)
        if group in {"pos", "neg"}:
            group = f"{group}_patterns"
        return group, name
    raise ValueError(f"Could not parse pattern key {key!r}")


def parse_report(report_path):
    """Parse report HTML and return soup/table metadata."""
    from bs4 import BeautifulSoup

    report_path = Path(report_path)
    soup = BeautifulSoup(report_path.read_text(errors="replace"), "html.parser")
    table = soup.find("table")
    if table is None:
        raise SystemExit(f"No <table> found in {report_path}")
    thead = table.find("thead")
    tbody = table.find("tbody")
    if thead is None or tbody is None:
        raise SystemExit("Expected report table to contain <thead> and <tbody>")
    header_row = thead.find("tr")
    headers = [th.get_text(strip=True) for th in header_row.find_all("th")]
    if not headers:
        raise SystemExit("Could not parse report table headers")

    match_cols = [
        idx for idx, header in enumerate(headers)
        if header.lower() in {"match0", "match1", "match2"}
    ]
    if not match_cols:
        raise SystemExit("Could not find match0/match1/match2 columns in report")

    rows = []
    for tr in tbody.find_all("tr"):
        cells = tr.find_all("td")
        if not cells:
            continue
        key = cells[0].get_text(strip=True)
        if not key:
            continue
        match_ids = []
        for idx in match_cols:
            if idx < len(cells):
                value = cells[idx].get_text(strip=True)
                if value and value.lower() != "nan":
                    match_ids.append(value)
        rows.append({
            "key": key,
            "cells": cells,
            "tr": tr,
            "match_ids": match_ids,
            "image_srcs": [img.get("src") for img in tr.find_all("img") if img.get("src")],
        })

    return {
        "soup": soup,
        "table": table,
        "headers": headers,
        "header_row": header_row,
        "rows": rows,
        "report_path": report_path,
    }


def n_subpatterns(pattern_group):
    """Count subpattern groups in a modisco-lite pattern group."""
    return len([k for k in pattern_group.keys() if k.startswith("subpattern_")])


def pfm_to_ic(pfm):
    """Scale a probability PFM to information-content heights."""
    import numpy as np

    p = np.clip(pfm, 1e-9, 1.0)
    p = p / p.sum(axis=1, keepdims=True)
    entropy = -(p * np.log2(p)).sum(axis=1)
    info = 2.0 - entropy
    return p * info[:, None]


def revcomp_pfm(pfm):
    """Reverse-complement a PFM with columns [A,C,G,T]."""
    return pfm[::-1, ::-1]


def trimmed_pfm(pattern_group, trim_frac):
    """Trim PFM to the contribution-supported core."""
    import numpy as np

    pfm = pattern_group["sequence"][:]
    cwm = pattern_group["contrib_scores"][:]
    if pfm.shape[0] == 4 and pfm.shape[1] != 4:
        pfm = pfm.T
    if cwm.shape[0] == 4 and cwm.shape[1] != 4:
        cwm = cwm.T
    score = np.abs(cwm).sum(axis=1)
    if len(score) == 0 or score.max() <= 0:
        return pfm
    keep = np.where(score >= trim_frac * score.max())[0]
    lo, hi = (keep.min(), keep.max() + 1) if len(keep) else (0, len(score))
    return pfm[lo:hi]


def logo_png_bytes(matrix):
    """Render a logo matrix to PNG bytes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import logomaker

    df = pd.DataFrame(matrix, columns=BASES)
    fig, ax = plt.subplots(figsize=(3.4, 0.95))
    logomaker.Logo(df, ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return buf.getvalue()


def data_uri(png_bytes):
    encoded = base64.b64encode(png_bytes).decode()
    return f"data:image/png;base64,{encoded}"


def resolve_image_bytes(src, report_path):
    """Resolve HTML image src to bytes, supporting data URIs and local relative paths."""
    if not src:
        return None, "empty src"
    if src.startswith("data:image"):
        try:
            _, payload = src.split(",", 1)
            return base64.b64decode(payload), None
        except Exception as exc:  # noqa: BLE001
            return None, f"could not decode data URI: {exc}"

    parsed = urlparse(src)
    if parsed.scheme in {"http", "https"}:
        return None, "remote image URL not embedded automatically"
    if parsed.scheme == "file":
        candidate = Path(unquote(parsed.path))
    else:
        clean = unquote(parsed.path)
        candidate = (Path(report_path).parent / clean).resolve()

    if not candidate.exists():
        return None, f"missing image file: {candidate}"
    try:
        return candidate.read_bytes(), None
    except OSError as exc:
        return None, f"could not read image file {candidate}: {exc}"


def image_array_from_bytes(image_bytes):
    """Read PNG/JPEG/SVG-compatible image bytes through matplotlib when possible."""
    import matplotlib.image as mpimg

    try:
        return mpimg.imread(io.BytesIO(image_bytes))
    except Exception:  # noqa: BLE001
        return None


def pattern_summary(row, id2name, annotations, h5_file, trim_frac):
    """Build all derived information for a report row."""
    key = row["key"]
    group, pattern_name = normalize_pattern_key(key)
    pattern_group = h5_file[group][pattern_name]
    pfm = trimmed_pfm(pattern_group, trim_frac)
    fwd_png = logo_png_bytes(pfm_to_ic(pfm))
    rev_png = logo_png_bytes(pfm_to_ic(revcomp_pfm(pfm)))
    tf_matches = [
        (motif_id, id2name.get(motif_id, "?"))
        for motif_id in row["match_ids"]
    ]
    return {
        "key": key,
        "group": group,
        "pattern_name": pattern_name,
        "tf_matches": tf_matches,
        "n_subpatterns": n_subpatterns(pattern_group),
        "annotation": annotations.get(key, annotations.get(f"{group}.{pattern_name}", "")),
        "fwd_png": fwd_png,
        "rev_png": rev_png,
    }


def write_upgraded_html(report, summaries, out_html, embed_images=False):
    """Write HTML with appended annotation columns; optionally inline original images."""
    from bs4 import BeautifulSoup

    soup = report["soup"]
    headers = report["headers"]
    header_row = report["header_row"]
    existing = set(headers)
    already_upgraded = all(header in existing for header in ADDED_HEADERS)
    if already_upgraded:
        print("[info] report already contains upgraded columns; not appending them again")
    else:
        for header in ADDED_HEADERS:
            th = soup.new_tag("th")
            th.string = header
            header_row.append(th)

        by_key = {summary["key"]: summary for summary in summaries}
        for row in report["rows"]:
            summary = by_key[row["key"]]
            tr = row["tr"]
            tf_html = "<br>".join(
                f"{motif_id} &rarr; <b>{tf_name}</b>"
                for motif_id, tf_name in summary["tf_matches"]
            )
            values = [
                f'<img src="{data_uri(summary["fwd_png"])}" width="240">',
                f'<img src="{data_uri(summary["rev_png"])}" width="240">',
                tf_html,
                str(summary["n_subpatterns"]),
                summary["annotation"],
            ]
            for value in values:
                td = soup.new_tag("td")
                td.append(BeautifulSoup(value, "html.parser"))
                tr.append(td)

    if embed_images:
        for img in soup.find_all("img"):
            src = img.get("src")
            if not src or src.startswith("data:image"):
                continue
            image_bytes, warning = resolve_image_bytes(src, report["report_path"])
            if image_bytes is None:
                print(f"[warn] {warning}", file=sys.stderr)
                continue
            img["src"] = data_uri(image_bytes)

    Path(out_html).write_text(str(soup))
    print(f"[done] wrote HTML: {out_html}")


def draw_image_panel(fig, rect, image_bytes, title, missing_text=None):
    """Draw one image slot on the PDF page."""
    import matplotlib.pyplot as plt

    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.set_title(title, fontsize=8)
    if image_bytes is None:
        ax.text(0.5, 0.5, missing_text or "image unavailable",
                ha="center", va="center", fontsize=8, wrap=True)
        return
    arr = image_array_from_bytes(image_bytes)
    if arr is None:
        ax.text(0.5, 0.5, "image format unreadable",
                ha="center", va="center", fontsize=8, wrap=True)
        return
    ax.imshow(arr)


def write_pdf(report, summaries, out_pdf, max_original_images=4, max_pages=None):
    """Write a standalone PDF, embedding all images drawn on each page."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    by_key = {summary["key"]: summary for summary in summaries}
    warnings = []
    with PdfPages(out_pdf) as pdf:
        for page_idx, row in enumerate(report["rows"]):
            if max_pages is not None and page_idx >= max_pages:
                break
            summary = by_key[row["key"]]
            fig = plt.figure(figsize=(11, 8.5))
            fig.patch.set_facecolor("white")

            tf_text = ", ".join(
                f"{motif_id} → {tf_name}"
                for motif_id, tf_name in summary["tf_matches"]
            ) or "No motif matches parsed"
            header = (
                f"{summary['key']}\n"
                f"Matched TFs: {tf_text}\n"
                f"Subpatterns: {summary['n_subpatterns']}"
            )
            if summary["annotation"]:
                header += f"\nAnnotation: {summary['annotation']}"
            fig.text(0.05, 0.95, header, ha="left", va="top", fontsize=11)

            slots = [
                (0.05, 0.58, 0.42, 0.25),
                (0.53, 0.58, 0.42, 0.25),
                (0.05, 0.28, 0.42, 0.25),
                (0.53, 0.28, 0.42, 0.25),
                (0.05, 0.05, 0.42, 0.16),
                (0.53, 0.05, 0.42, 0.16),
            ]

            original_srcs = row["image_srcs"][:max_original_images]
            for idx, src in enumerate(original_srcs[:4]):
                image_bytes, warning = resolve_image_bytes(src, report["report_path"])
                if warning:
                    warnings.append(f"{summary['key']}: {warning}")
                draw_image_panel(
                    fig, slots[idx], image_bytes,
                    title=f"Original report image {idx + 1}",
                    missing_text=warning,
                )

            draw_image_panel(fig, slots[4], summary["fwd_png"], "PFM forward")
            draw_image_panel(fig, slots[5], summary["rev_png"], "PFM reverse")

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"[done] wrote PDF: {out_pdf}")
    if warnings:
        print(f"[warn] {len(warnings)} report image(s) could not be embedded", file=sys.stderr)
        for warning in warnings[:20]:
            print(f"[warn] {warning}", file=sys.stderr)


def infer_outputs(args):
    """Resolve output paths while keeping --out backward compatible."""
    out_html = args.out_html
    out_pdf = args.out_pdf
    if args.out:
        suffix = Path(args.out).suffix.lower()
        if suffix == ".pdf":
            out_pdf = args.out
        else:
            out_html = args.out
    if not out_html and not out_pdf:
        report = Path(args.report)
        out_pdf = str(report.with_name(f"{report.stem}_upgraded.pdf"))
    return out_html, out_pdf


def main(argv=None):
    args = parse_args(argv)
    out_html, out_pdf = infer_outputs(args)

    _optional_hdf5plugin()
    import h5py

    id2name = parse_meme_names(args.meme)
    annotations = load_notes(args.notes)
    report = parse_report(args.report)

    summaries = []
    with h5py.File(args.h5, "r") as h5_file:
        for row in report["rows"]:
            summaries.append(pattern_summary(row, id2name, annotations, h5_file, args.trim_frac))

    print(f"[info] parsed {len(report['rows'])} report rows")
    print(f"[info] parsed {len(id2name)} MEME motif names")
    if out_html:
        write_upgraded_html(report, summaries, out_html, embed_images=args.embed_html_images)
    if out_pdf:
        write_pdf(
            report,
            summaries,
            out_pdf,
            max_original_images=args.max_original_images,
            max_pages=args.max_pages,
        )


if __name__ == "__main__":
    main()
