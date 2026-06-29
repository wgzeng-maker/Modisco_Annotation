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

Requirements: h5py, numpy, pandas
Optional: hdf5plugin for compressed production HDF5 files; pyfaidx for --verify
"""

import argparse
from contextlib import contextmanager

try:
    import hdf5plugin   # noqa: F401  -- imported before h5py when available
    HDF5PLUGIN_AVAILABLE = True
except ModuleNotFoundError:
    HDF5PLUGIN_AVAILABLE = False

import h5py
import numpy as np
import pandas as pd


@contextmanager
def open_h5(path):
    """Open an HDF5 file, with a clear message for compressed modisco-lite files."""
    try:
        with h5py.File(path, "r") as f:
            yield f
    except OSError as exc:
        msg = str(exc).lower()
        likely_filter_error = "filter" in msg or "plugin" in msg or "blosc" in msg
        if likely_filter_error and not HDF5PLUGIN_AVAILABLE:
            raise OSError(
                f"Could not open {path!r}. This file may use HDF5 compression "
                "filters that require the optional package 'hdf5plugin'. Install "
                "it in the active environment, then rerun this command."
            ) from exc
        raise


def _decode(onehot_2d):
    """(width, 4) one-hot -> ACGT string."""
    return "".join("ACGT"[i] for i in onehot_2d.argmax(axis=1))


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _natural_key(name):
    """Sort pattern_2 before pattern_10 while keeping unknown names stable."""
    prefix, _, suffix = name.rpartition("_")
    if suffix.isdigit():
        return prefix, int(suffix)
    return name, -1


def _read_n_seqlets(seqlets_group):
    """Read modisco-lite seqlet count robustly across scalar/1-element layouts."""
    if "n_seqlets" not in seqlets_group:
        return len(seqlets_group["example_idx"])
    value = seqlets_group["n_seqlets"][()]
    return int(np.asarray(value).reshape(-1)[0])


def _normalize_seqlet_matrix(onehot_2d):
    """Return a seqlet sequence as (width, 4), accepting (4, width) if needed."""
    arr = np.asarray(onehot_2d)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D one-hot seqlet matrix; got shape {arr.shape}")
    if arr.shape[-1] == 4:
        return arr
    if arr.shape[0] == 4:
        return arr.T
    raise ValueError(f"Expected one dimension of seqlet matrix to be 4; got {arr.shape}")


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
    with open_h5(modisco_h5) as f:
        for arm in ["pos_patterns", "neg_patterns"]:
            if arm not in f:
                continue
            for pat in sorted(f[arm], key=_natural_key):
                parent = f[arm][pat]
                subs = sorted(
                    [k for k in parent if k.startswith("subpattern_")],
                    key=_natural_key,
                )
                # If a pattern somehow has no subpatterns, fall back to its
                # top-level seqlets so nothing is silently dropped.
                groups = subs if subs else ["__self__"]
                if subs and "seqlets" in parent:
                    parent_n = _read_n_seqlets(parent["seqlets"])
                    sub_n = sum(_read_n_seqlets(parent[sub]["seqlets"]) for sub in subs)
                    if sub_n != parent_n:
                        raise ValueError(
                            f"{arm}/{pat}: subpatterns contain {sub_n} seqlets, "
                            f"but parent pattern reports {parent_n}. Refusing to "
                            "write an incomplete or duplicated annotation table."
                        )
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


def _verification_candidates(h5_file, max_checks):
    """Pick representative seqlets across arms/patterns/subpatterns/strands."""
    candidates = []
    for arm in ["pos_patterns", "neg_patterns"]:
        if arm not in h5_file:
            continue
        for pat in sorted(h5_file[arm], key=_natural_key):
            parent = h5_file[arm][pat]
            subs = sorted(
                [k for k in parent if k.startswith("subpattern_")],
                key=_natural_key,
            )
            groups = subs if subs else ["__self__"]
            for sub in groups:
                node = parent if sub == "__self__" else parent[sub]
                sl = node["seqlets"]
                n = len(sl["example_idx"])
                if n == 0:
                    continue
                rc = np.asarray(sl["is_revcomp"][:], dtype=bool)
                chosen = []
                for want_rc in [False, True]:
                    idx = np.flatnonzero(rc == want_rc)
                    if len(idx):
                        chosen.append(int(idx[0]))
                if not chosen:
                    chosen.append(0)
                for idx in dict.fromkeys(chosen):
                    candidates.append((arm, pat, sub, sl, idx))
                    if len(candidates) >= max_checks:
                        return candidates
    return candidates


def verify_frame(modisco_h5, genome_fa, chrom, window_start, crop_offset, max_checks=20):
    """Decode sampled seqlets and match them against the genome FASTA. Returns bool."""
    from pyfaidx import Fasta
    genome = Fasta(genome_fa)
    failures = []
    checked = 0
    with open_h5(modisco_h5) as f:
        candidates = _verification_candidates(f, max_checks=max_checks)
        if not candidates:
            raise ValueError("No seqlets found to verify.")
        for arm, pat, sub, sl, idx in candidates:
            ex = int(sl["example_idx"][idx])
            s = int(sl["start"][idx])
            e = int(sl["end"][idx])
            rc = bool(sl["is_revcomp"][idx])
            stored = _normalize_seqlet_matrix(sl["sequence"][idx])
            stored_str = _decode(stored)

            g_s = int(window_start[ex] + crop_offset + s)
            g_e = int(window_start[ex] + crop_offset + e)
            fetched = str(genome[str(chrom[ex])][g_s:g_e]).upper()
            expected = _revcomp(fetched) if rc else fetched
            ok = stored_str == expected
            checked += 1
            label = f"{arm}/{pat}/{sub}[{idx}]"
            print(f"  verify: {label} example_idx={ex} revcomp={rc} {g_s}-{g_e} MATCH={ok}")
            if not ok:
                failures.append((label, stored_str, expected))

    print(f"  verified {checked} seqlets against {genome_fa}")
    if failures:
        print("  failed examples:")
        for label, stored_str, expected in failures[:5]:
            print(f"    {label}")
            print(f"      stored : {stored_str}")
            print(f"      genome : {expected}")
        return False
    return True


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
    ap.add_argument("--verify-max-seqlets", type=int, default=20,
                    help="maximum sampled seqlets to check with --genome (default 20)")
    args = ap.parse_args()

    crop_offset = (args.input_len - args.window) // 2
    if crop_offset < 0:
        raise SystemExit("--input-len must be greater than or equal to --window")
    print(f"crop offset = ({args.input_len} - {args.window}) // 2 = {crop_offset}")

    chrom, window_start = load_peak_windows(args.bed, args.input_len)
    print(f"loaded {len(chrom):,} peak windows")

    if args.genome:
        ok = verify_frame(
            args.modisco,
            args.genome,
            chrom,
            window_start,
            crop_offset,
            max_checks=args.verify_max_seqlets,
        )
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
