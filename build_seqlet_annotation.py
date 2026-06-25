"""
build_seqlet_annotation.py
==========================
Map every clustered MoDISco seqlet back to its genome coordinates, tagged with
its pattern and subpattern. Produces the master annotation table that all
downstream analysis (co-occurrence, gene overlap, cross-model comparison) builds
on.

Coordinate logic (verified against the mm10 genome FASTA)
---------------------------------------------------------
A seqlet's start/end are relative to the MoDISco motif-discovery window, which
is the central `-w` crop of the full ChromBPNet input window. To recover genome
coordinates:

    summit            = peak_chrom_start + summit_offset      (narrowPeak col 9)
    window_genome_start = summit - input_len // 2             (window centered on summit)
    seqlet_genome_start = window_genome_start + crop_offset + seqlet.start
    crop_offset         = (input_len - w) // 2

Strand is '-' if the seqlet's is_revcomp flag is set, else '+'. For revcomp
seqlets the stored sequence is the reverse-complement of the genome.

ALWAYS verify the coordinate frame on a new dataset with --verify (decodes one
seqlet and matches it against the genome FASTA) before trusting the table.

Usage
-----
    python build_seqlet_annotation.py \\
        --modisco GC_modisco_profile_v2.h5 \\
        --bed     GC_mm10.interpreted_regions.bed \\
        --output  GC_profile_seqlet_annotation \\
        --window  1000 \\
        --input-len 2114 \\
        --genome  mm10.fa        # optional, enables --verify

Requirements: h5py, hdf5plugin, numpy, pandas  (pyfaidx only for --verify)
"""

import argparse
import hdf5plugin   # noqa: F401  -- before h5py
import h5py
import numpy as np
import pandas as pd


def _decode(onehot_2d):
    """(width, 4) one-hot -> ACGT string."""
    return "".join("ACGT"[i] for i in onehot_2d.argmax(axis=1))


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def load_peak_windows(bed_path, input_len):
    """From a narrowPeak-style BED, compute each window's genomic start.

    Returns (chrom array, window_start array), indexed by example_idx (row order).
    """
    bed = pd.read_csv(bed_path, sep="\t", header=None)
    if bed.shape[1] < 10:
        raise ValueError(
            f"Expected a 10-column narrowPeak BED; got {bed.shape[1]} columns. "
            "Column 9 (summit offset) is required.")
    chrom = bed[0].values
    summit = bed[1].values + bed[9].values          # genomic summit
    window_start = summit - (input_len // 2)         # window centered on summit
    return chrom, window_start


def build_annotation(modisco_h5, chrom, window_start, crop_offset):
    """Walk every seqlet in every pattern/subpattern -> annotation DataFrame."""
    rows = []
    with h5py.File(modisco_h5, "r") as f:
        for arm in ["pos_patterns", "neg_patterns"]:
            if arm not in f:
                continue
            for pat in f[arm]:
                parent = f[arm][pat]
                subs = [k for k in parent if k.startswith("subpattern_")]
                # If a pattern somehow has no subpatterns, fall back to its
                # top-level seqlets so nothing is silently dropped.
                groups = subs if subs else ["__self__"]
                for sub in groups:
                    node = parent if sub == "__self__" else parent[sub]
                    sl = node["seqlets"]
                    ex = sl["example_idx"][:]
                    s = sl["start"][:]
                    e = sl["end"][:]
                    rc = sl["is_revcomp"][:]
                    g_start = window_start[ex] + crop_offset + s
                    g_end = window_start[ex] + crop_offset + e
                    for i in range(len(ex)):
                        rows.append((
                            chrom[ex[i]], int(g_start[i]), int(g_end[i]),
                            "-" if rc[i] else "+",
                            f"{arm.split('_')[0]}/{pat}",
                            "" if sub == "__self__" else sub,
                            int(ex[i]),
                        ))
    return pd.DataFrame(rows, columns=[
        "chrom", "start", "end", "strand",
        "pattern", "subpattern", "example_idx"])


def verify_frame(modisco_h5, genome_fa, chrom, window_start, crop_offset):
    """Decode one seqlet and match it against the genome FASTA. Returns bool."""
    from pyfaidx import Fasta
    genome = Fasta(genome_fa)
    with h5py.File(modisco_h5, "r") as f:
        # use the first pattern/subpattern available
        arm = "pos_patterns" if "pos_patterns" in f else "neg_patterns"
        pat = list(f[arm].keys())[0]
        parent = f[arm][pat]
        sub = next((k for k in parent if k.startswith("subpattern_")), None)
        node = parent[sub] if sub else parent
        sl = node["seqlets"]
        ex = int(sl["example_idx"][0]); s = int(sl["start"][0]); e = int(sl["end"][0])
        rc = bool(sl["is_revcomp"][0])
        stored = sl["sequence"][0]
    if stored.shape[0] == 4:
        stored = stored.T
    stored_str = _decode(stored)

    g_s = window_start[ex] + crop_offset + s
    g_e = window_start[ex] + crop_offset + e
    fetched = str(genome[chrom[ex]][g_s:g_e]).upper()
    expected = _revcomp(fetched) if rc else fetched

    ok = stored_str == expected
    print(f"  verify: example_idx={ex} revcomp={rc}")
    print(f"    stored : {stored_str}")
    print(f"    genome : {expected}")
    print(f"    MATCH  : {ok}")
    return ok


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modisco", required=True, help="modisco-lite output .h5")
    ap.add_argument("--bed", required=True,
                    help="narrowPeak-style interpreted_regions.bed")
    ap.add_argument("--output", required=True,
                    help="output prefix (writes .bed and .csv)")
    ap.add_argument("--window", type=int, default=1000,
                    help="MoDISco -w value (default 1000)")
    ap.add_argument("--input-len", type=int, default=2114,
                    help="full ChromBPNet input length (default 2114)")
    ap.add_argument("--genome", default=None,
                    help="genome FASTA; if given, verifies the coordinate frame")
    args = ap.parse_args()

    crop_offset = (args.input_len - args.window) // 2
    print(f"crop offset = ({args.input_len} - {args.window}) // 2 = {crop_offset}")

    chrom, window_start = load_peak_windows(args.bed, args.input_len)
    print(f"loaded {len(chrom):,} peak windows")

    if args.genome:
        ok = verify_frame(args.modisco, args.genome, chrom, window_start, crop_offset)
        if not ok:
            raise SystemExit("Coordinate frame verification FAILED -- stopping. "
                             "Check --window / --input-len / BED order.")

    ann = build_annotation(args.modisco, chrom, window_start, crop_offset)
    print(f"annotated {len(ann):,} seqlets")

    if (ann["start"] < 0).any():
        print(f"WARNING: {(ann['start'] < 0).sum()} seqlets have negative starts "
              "(near a window edge / chromosome start)")

    # Write a BED (pattern|subpattern packed into name) and a CSV with headers.
    bed_out = ann.copy()
    bed_out["name"] = bed_out["pattern"] + "|" + bed_out["subpattern"]
    bed_out["score"] = 0
    bed_out = bed_out[["chrom", "start", "end", "name", "score", "strand",
                       "pattern", "subpattern", "example_idx"]]
    bed_out.to_csv(args.output + ".bed", sep="\t", header=False, index=False)
    ann.to_csv(args.output + ".csv", index=False)
    print(f"saved {args.output}.bed and {args.output}.csv")

    # Quick coverage summary
    n_peaks = ann["example_idx"].nunique()
    print(f"\npeaks with >=1 annotated seqlet: {n_peaks:,}")


if __name__ == "__main__":
    main()
