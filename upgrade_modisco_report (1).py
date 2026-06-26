#!/usr/bin/env python
"""
upgrade_modisco_report.py

Append annotation columns onto an existing TF-MoDISco report (motifs.html),
keeping every original column intact. Adds, at the right edge:

    PFM forward | PFM reverse | matched TFs | #subpatterns | functional annotation

  - PFM logos are drawn from the modisco .h5 'sequence' arrays (base64-embedded,
    so they are self-contained and portable).
  - "matched TFs" reuses the match IDs already present in your report and resolves
    them to TF names using the JASPAR .meme header lines (e.g. MA0671.2 -> NFIX).
    No re-running of TOMTOM.
  - "#subpatterns" is read from the .h5.
  - "functional annotation" is left blank (hook for downstream ontology work; see
    --notes to populate it from a CSV).

Example:
    python upgrade_modisco_report.py \\
        --report  GC_contribs/profile_report/motifs.html \\
        --h5      modisco_results.h5 \\
        --meme    JASPAR2024_CORE_non-redundant_pfms_meme.txt \\
        --out     motifs_upgraded.html

Dependencies: numpy, pandas, h5py, logomaker, matplotlib, beautifulsoup4
"""

import argparse
import base64
import io
import sys

import h5py
import numpy as np
import pandas as pd
import logomaker
import matplotlib
matplotlib.use("Agg")  # headless: no display needed
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

BASES = ["A", "C", "G", "T"]


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Append PFM logos, TF names, and subpattern counts to an existing modisco report.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--report", required=True,
                   help="Existing modisco report HTML (e.g. profile_report/motifs.html).")
    p.add_argument("--h5", required=True,
                   help="modisco results .h5 (provides PFM 'sequence' arrays and subpatterns).")
    p.add_argument("--meme", required=True,
                   help="JASPAR .meme file used for the original report (for ID -> TF name).")
    p.add_argument("--out", default=None,
                   help="Output HTML path. Default: <report stem>_upgraded.html")
    p.add_argument("--trim-frac", type=float, default=0.30,
                   help="Trim the PFM logo to columns whose CWM signal >= this fraction of the max.")
    p.add_argument("--notes", default=None,
                   help="Optional CSV with columns 'pattern,annotation' to fill the last column.")
    return p.parse_args(argv)


def parse_meme_names(path):
    """Map JASPAR matrix ID -> TF name from 'MOTIF <id> <name>' header lines."""
    id2name = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("MOTIF"):
                parts = line.split()
                id2name[parts[1]] = " ".join(parts[2:]) if len(parts) > 2 else parts[1]
    return id2name


def load_notes(path):
    """Optional pattern -> annotation map for the last column."""
    if not path:
        return {}
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    pcol = cols.get("pattern")
    acol = cols.get("annotation")
    if not pcol or not acol:
        sys.exit("--notes CSV must have 'pattern' and 'annotation' columns.")
    return {str(r[pcol]).strip(): str(r[acol]) for _, r in df.iterrows()}


# --------------------------------------------------------------------------- #
# Logos
# --------------------------------------------------------------------------- #
def to_df(mat):
    return pd.DataFrame(mat, columns=BASES)


def pfm_to_ic(pfm):
    """Scale a probability PFM to information-content (bits) for a standard logo."""
    p = np.clip(pfm, 1e-9, 1.0)
    p = p / p.sum(axis=1, keepdims=True)
    ic = 2.0 - (-(p * np.log2(p)).sum(axis=1))
    return p * ic[:, None]


def revcomp_pfm(pfm):
    """Reverse-complement: columns are [A,C,G,T], so reverse rows AND columns."""
    return pfm[::-1, ::-1]


def logo_b64(df):
    fig, ax = plt.subplots(figsize=(3.4, 0.95))
    logomaker.Logo(df, ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def img_tag(b64):
    return f'<img src="data:image/png;base64,{b64}" width="240">'


# --------------------------------------------------------------------------- #
# h5 helpers
# --------------------------------------------------------------------------- #
def n_subpatterns(patt):
    subs = [k for k in patt.keys() if k.startswith("subpattern")]
    if not subs and "subpatterns" in patt:
        subs = list(patt["subpatterns"].keys())
    return len(subs)


def trimmed_pfm(patt, trim_frac):
    pfm = patt["sequence"][:]
    cwm = patt["contrib_scores"][:]
    score = np.abs(cwm).sum(axis=1)
    keep = np.where(score >= trim_frac * score.max())[0]
    lo, hi = (keep.min(), keep.max() + 1) if len(keep) else (0, len(score))
    return pfm[lo:hi]


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)
    out_path = args.out
    if out_path is None:
        stem = args.report.rsplit(".", 1)[0]
        out_path = f"{stem}_upgraded.html"

    id2name = parse_meme_names(args.meme)
    notes = load_notes(args.notes)

    with open(args.report) as f:
        soup = BeautifulSoup(f.read(), "html.parser")
    table = soup.find("table")
    if table is None:
        sys.exit(f"No <table> found in {args.report}")
    thead_tr = table.find("thead").find("tr")
    headers = [th.get_text(strip=True) for th in thead_tr.find_all("th")]
    match_cols = [headers.index(c) for c in ("match0", "match1", "match2") if c in headers]
    if not match_cols:
        sys.exit("Could not find any match0/1/2 columns in the report header.")

    print(f"[verify] headers: {headers}")
    print(f"[verify] match-ID columns at: {match_cols}")
    print(f"[verify] meme names: {len(id2name)} (MA0671.2 -> {id2name.get('MA0671.2', '??')})")

    for h in ("PFM forward", "PFM reverse", "matched TFs", "#subpatterns", "functional annotation"):
        th = soup.new_tag("th")
        th.string = h
        thead_tr.append(th)

    rows = table.find("tbody").find_all("tr")
    with h5py.File(args.h5, "r") as hf:
        for ri, tr in enumerate(rows):
            tds = tr.find_all("td")
            key = tds[0].get_text(strip=True)          # e.g. pos_patterns.pattern_0
            group, name = key.split(".")
            patt = hf[group][name]

            pfm_t = trimmed_pfm(patt, args.trim_frac)
            fwd = img_tag(logo_b64(to_df(pfm_to_ic(pfm_t))))
            rev = img_tag(logo_b64(to_df(pfm_to_ic(revcomp_pfm(pfm_t)))))

            tf_bits = []
            for ci in match_cols:
                mid = tds[ci].get_text(strip=True)
                if mid and mid.lower() != "nan":
                    tf_bits.append(f"{mid} &rarr; <b>{id2name.get(mid, '?')}</b>")
            tf_html = "<br>".join(tf_bits)

            for content in (fwd, rev, tf_html, str(n_subpatterns(patt)), notes.get(key, "")):
                td = soup.new_tag("td")
                td.append(BeautifulSoup(content, "html.parser"))
                tr.append(td)

            if ri == 0:
                print(f"[verify] first row TFs: {tf_bits}")

    with open(out_path, "w") as f:
        f.write(str(soup))
    print(f"[done] wrote {out_path} ({len(rows)} patterns)")


if __name__ == "__main__":
    main()
