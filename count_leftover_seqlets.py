"""
count_leftover_seqlets.py
=========================
Reproduce the MoDISco seqlet "funnel" for a ChromBPNet run:

    candidates identified  ->  survived pos/neg threshold  ->  in final patterns

This quantifies how many identified seqlets never enter clustering (because of
the max_seqlets_per_metacluster cap) and how many end up in final patterns.
See LIMITATION_seqlet_coverage.md for why this matters.

REQUIRES modiscolite (run inside the kundajelab/chrombpnet Docker image, or any
env where `import modiscolite` works):

    docker run --rm -v /path/to/data:/work kundajelab/chrombpnet:latest \\
        python3 count_leftover_seqlets.py \\
            --scores /work/contribs/MODEL.profile_scores.h5 \\
            --modisco /work/contribs/MODEL_modisco_profile.h5

The two inputs:
  --scores   : the contribution-score .h5 fed to MoDISco (has raw/seq + shap/seq)
  --modisco  : the modisco-lite output .h5 (to count seqlets in final patterns)

These default parameters match modisco-lite's TFMoDISco() defaults. If you ran
MoDISco with non-default values, override them with the matching flags.
"""

import argparse
import numpy as np


def count_candidate_seqlets(scores_h5, window_size, flank, target_fdr,
                            min_passing_windows_frac, max_passing_windows_frac,
                            weak_threshold):
    """Run ONLY seqlet identification + the pos/neg threshold split.

    Returns (n_candidates, n_pos, n_neg). No clustering is performed.
    """
    import hdf5plugin   # noqa: F401  -- before h5py
    import h5py
    from modiscolite.extract_seqlets import extract_seqlets
    from modiscolite import core

    with h5py.File(scores_h5, "r") as f:
        one_hot = f["raw"]["seq"][:]
        hyp = f["shap"]["seq"][:]

    # modiscolite expects arrays as (n, length, 4). ChromBPNet writes
    # (n, 4, length), so swap the last two axes if needed.
    if one_hot.shape[1] == 4:
        one_hot = one_hot.transpose(0, 2, 1)
        hyp = hyp.transpose(0, 2, 1)

    contrib = one_hot * hyp                      # projected contribution scores
    suppress = int(0.5 * window_size) + flank

    # --- Stage 1: identify candidate seqlets ---
    seqlet_coords, threshold = extract_seqlets(
        attribution_scores=contrib.sum(axis=2),  # (n, length)
        window_size=window_size, flank=flank, suppress=suppress,
        target_fdr=target_fdr,
        min_passing_windows_frac=min_passing_windows_frac,
        max_passing_windows_frac=max_passing_windows_frac,
        weak_threshold_for_counting_sign=weak_threshold)
    n_candidates = len(seqlet_coords)

    # --- Stage 2: build seqlet objects + apply the pos/neg threshold split ---
    track_set = core.TrackSet(one_hot=one_hot, contrib_scores=contrib,
                              hypothetical_contribs=hyp)
    seqlets = track_set.create_seqlets(seqlet_coords)

    n_pos = n_neg = 0
    for seqlet in seqlets:
        fl = int(0.5 * (len(seqlet) - window_size))
        attr = np.sum(seqlet.contrib_scores[fl:-fl])   # core only, all 4 bases
        if attr > threshold:
            n_pos += 1
        elif attr < -threshold:
            n_neg += 1

    return n_candidates, n_pos, n_neg


def count_pattern_seqlets(modisco_h5):
    """Count seqlets assigned to final patterns in a modisco-lite output."""
    import hdf5plugin   # noqa: F401
    import h5py

    total = 0
    with h5py.File(modisco_h5, "r") as f:
        for arm in ["pos_patterns", "neg_patterns"]:
            if arm not in f:
                continue
            for p in f[arm]:
                total += int(f[arm][p]["seqlets"]["n_seqlets"][()][0])
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True,
                    help="contribution-score .h5 fed to MoDISco")
    ap.add_argument("--modisco", required=True,
                    help="modisco-lite output .h5")
    # Parameters matching modisco-lite TFMoDISco() defaults
    ap.add_argument("--window_size", type=int, default=21)
    ap.add_argument("--flank", type=int, default=10)
    ap.add_argument("--target_fdr", type=float, default=0.2)
    ap.add_argument("--min_passing_windows_frac", type=float, default=0.03)
    ap.add_argument("--max_passing_windows_frac", type=float, default=0.2)
    ap.add_argument("--weak_threshold", type=float, default=0.8)
    args = ap.parse_args()

    n_cand, n_pos, n_neg = count_candidate_seqlets(
        args.scores, args.window_size, args.flank, args.target_fdr,
        args.min_passing_windows_frac, args.max_passing_windows_frac,
        args.weak_threshold)
    n_survived = n_pos + n_neg
    n_patterns = count_pattern_seqlets(args.modisco)

    print("\n===== MoDISco seqlet funnel =====")
    print(f"  candidates identified     : {n_cand:>12,}")
    print(f"  survived pos/neg threshold: {n_survived:>12,}  "
          f"(pos {n_pos:,} / neg {n_neg:,})")
    print(f"  assigned to final patterns: {n_patterns:>12,}")
    if n_cand:
        print(f"  -> final / candidates     : {100*n_patterns/n_cand:>11.2f}%")
    print("\nMost of the drop is the 20,000-per-metacluster cap, not the")
    print("importance threshold. See LIMITATION_seqlet_coverage.md.")


if __name__ == "__main__":
    main()
