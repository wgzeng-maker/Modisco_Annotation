"""
count_leftover_seqlets.py
=========================
Reproduce the measurable parts of the MoDISco seqlet "funnel" for a
ChromBPNet run:

    candidates identified -> survived pos/neg threshold -> in final patterns

This script runs seqlet identification and the positive/negative threshold
split, then counts how many seqlets were saved in final MoDISco patterns. It
does not re-run clustering. If you provide the `-n/--max-seqlets-per-metacluster`
value from the real MoDISco run, it also reports the maximum number of threshold
survivors that could have entered clustering after that cap.

REQUIRES modiscolite. The safest environment is the same container used for
ChromBPNet/MoDISco:

    docker run --rm -v /path/to/data:/work kundajelab/chrombpnet:latest \\
        python3 count_leftover_seqlets.py \\
            --scores /work/contribs/MODEL.profile_scores.h5 \\
            --modisco /work/contribs/MODEL_modisco_profile.h5 \\
            --window 1000 \\
            -n 100000

The three key inputs:
  --scores   : contribution-score .h5 fed to MoDISco (raw/seq + shap/seq)
  --modisco  : modisco-lite output .h5 (to count final pattern seqlets)
  --window   : the -w value used in the real run. CRITICAL: using the wrong
               window inflates the candidate count.

These default parameters match modisco-lite's TFMoDISco() defaults. If you ran
MoDISco with non-default values, override them with the matching flags.
"""

import argparse
import json


def _optional_hdf5plugin():
    try:
        import hdf5plugin  # noqa: F401
    except ModuleNotFoundError:
        pass


def count_candidate_seqlets(
    scores_h5,
    window,
    window_size,
    flank,
    target_fdr,
    min_passing_windows_frac,
    max_passing_windows_frac,
    weak_threshold,
):
    """Run only seqlet identification and the pos/neg threshold split.

    `window` is the MoDISco -w value: the central crop applied before motif
    discovery. It must match the real run.

    Returns (n_candidates, n_pos, n_neg, threshold). No clustering is performed.
    """
    import numpy as np

    _optional_hdf5plugin()
    import h5py
    from modiscolite import core
    from modiscolite.extract_seqlets import extract_seqlets

    with h5py.File(scores_h5, "r") as f:
        one_hot = f["raw"]["seq"][:]
        hyp = f["shap"]["seq"][:]

    # modiscolite expects (n, length, 4). ChromBPNet commonly writes
    # (n, 4, length), so swap the last two axes if needed.
    if one_hot.shape[1] == 4:
        one_hot = one_hot.transpose(0, 2, 1)
        hyp = hyp.transpose(0, 2, 1)

    full_len = one_hot.shape[1]
    if window > full_len:
        raise ValueError(f"--window ({window}) cannot exceed input length ({full_len})")
    if window < full_len:
        offset = (full_len - window) // 2
        one_hot = one_hot[:, offset:offset + window, :]
        hyp = hyp[:, offset:offset + window, :]

    contrib = one_hot * hyp
    suppress = int(0.5 * window_size) + flank

    seqlet_coords, threshold = extract_seqlets(
        attribution_scores=contrib.sum(axis=2),
        window_size=window_size,
        flank=flank,
        suppress=suppress,
        target_fdr=target_fdr,
        min_passing_windows_frac=min_passing_windows_frac,
        max_passing_windows_frac=max_passing_windows_frac,
        weak_threshold_for_counting_sign=weak_threshold,
    )
    n_candidates = len(seqlet_coords)

    track_set = core.TrackSet(
        one_hot=one_hot,
        contrib_scores=contrib,
        hypothetical_contribs=hyp,
    )
    seqlets = track_set.create_seqlets(seqlet_coords)

    n_pos = n_neg = 0
    for seqlet in seqlets:
        fl = int(0.5 * (len(seqlet) - window_size))
        core_scores = seqlet.contrib_scores if fl == 0 else seqlet.contrib_scores[fl:-fl]
        attr = np.sum(core_scores)
        if attr > threshold:
            n_pos += 1
        elif attr < -threshold:
            n_neg += 1

    return n_candidates, n_pos, n_neg, float(threshold)


def count_pattern_seqlets(modisco_h5):
    """Count seqlets assigned to final patterns in a modisco-lite output."""
    import numpy as np

    _optional_hdf5plugin()
    import h5py

    total = 0
    per_pattern = {}
    with h5py.File(modisco_h5, "r") as f:
        for arm in ["pos_patterns", "neg_patterns"]:
            if arm not in f:
                continue
            for p in f[arm]:
                sl = f[arm][p]["seqlets"]
                if "n_seqlets" in sl:
                    value = sl["n_seqlets"][()]
                    n = int(np.asarray(value).reshape(-1)[0])
                else:
                    n = len(sl["example_idx"])
                label = f"{arm.split('_')[0]}/{p}"
                per_pattern[label] = n
                total += n
    return total, per_pattern


def estimate_clustering_entrants(n_pos, n_neg, max_seqlets_per_metacluster):
    """Estimate the maximum threshold survivors allowed past the metacluster cap."""
    if max_seqlets_per_metacluster is None:
        return None
    return min(n_pos, max_seqlets_per_metacluster) + min(
        n_neg, max_seqlets_per_metacluster
    )


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--scores", required=True,
                        help="contribution-score .h5 fed to MoDISco")
    parser.add_argument("--modisco", required=True,
                        help="modisco-lite output .h5")
    parser.add_argument("--window", type=int, default=1000,
                        help="MoDISco -w central crop. Must match the real run.")
    parser.add_argument("-n", "--max-seqlets-per-metacluster", type=int, default=None,
                        help="MoDISco -n value from the real run. Used only to "
                             "estimate the max clustering input after the cap.")

    # Parameters matching modisco-lite TFMoDISco() defaults.
    parser.add_argument("--window-size", type=int, default=21)
    parser.add_argument("--flank", type=int, default=10)
    parser.add_argument("--target-fdr", type=float, default=0.2)
    parser.add_argument("--min-passing-windows-frac", type=float, default=0.03)
    parser.add_argument("--max-passing-windows-frac", type=float, default=0.2)
    parser.add_argument("--weak-threshold", type=float, default=0.8)
    parser.add_argument("--json", action="store_true",
                        help="print machine-readable JSON instead of text")
    return parser


def main(argv=None):
    args = build_arg_parser().parse_args(argv)

    n_cand, n_pos, n_neg, threshold = count_candidate_seqlets(
        args.scores,
        args.window,
        args.window_size,
        args.flank,
        args.target_fdr,
        args.min_passing_windows_frac,
        args.max_passing_windows_frac,
        args.weak_threshold,
    )
    n_survived = n_pos + n_neg
    n_patterns, per_pattern = count_pattern_seqlets(args.modisco)
    n_entered_est = estimate_clustering_entrants(
        n_pos,
        n_neg,
        args.max_seqlets_per_metacluster,
    )

    result = {
        "candidates_identified": n_cand,
        "threshold": threshold,
        "survived_pos_neg_threshold": n_survived,
        "positive_threshold_survivors": n_pos,
        "negative_threshold_survivors": n_neg,
        "max_seqlets_per_metacluster": args.max_seqlets_per_metacluster,
        "estimated_clustering_entrants_after_cap": n_entered_est,
        "assigned_to_final_patterns": n_patterns,
        "final_over_candidates_percent": (100 * n_patterns / n_cand) if n_cand else None,
        "per_pattern": per_pattern,
    }

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    print("\n===== MoDISco seqlet funnel =====")
    print(f"  candidates identified     : {n_cand:>12,}")
    print(f"  threshold value           : {threshold:>12.6g}")
    print(f"  survived pos/neg threshold: {n_survived:>12,}  "
          f"(pos {n_pos:,} / neg {n_neg:,})")
    if n_entered_est is not None:
        print(f"  estimated clustering input: {n_entered_est:>12,}  "
              f"(after -n {args.max_seqlets_per_metacluster:,} cap)")
    print(f"  assigned to final patterns: {n_patterns:>12,}")
    if n_cand:
        print(f"  -> final / candidates     : {100*n_patterns/n_cand:>11.2f}%")
    print("\nNote: this script does not rerun clustering. The cap line is an")
    print("estimate from pos/neg threshold survivors and the supplied -n value.")


if __name__ == "__main__":
    main()
