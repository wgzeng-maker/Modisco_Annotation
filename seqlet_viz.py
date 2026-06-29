r"""
seqlet_viz.py
=============
Visualization toolkit for TF-MoDISco seqlets from ChromBPNet models.

This script treats each seqlet like a point in a single-cell-style embedding:
the same seqlets can be compared in sequence space or attribution space. The
goal is to reveal when DNA sequence similarity and model-attribution similarity
agree or disagree within a MoDISco pattern.

Examples
--------
Atlas view across all patterns:

    python seqlet_viz.py atlas \
      --modisco GC_modisco_profile_v2.h5 \
      --output umap_all_patterns.png

Subpatterns inside one pattern:

    python seqlet_viz.py subpattern \
      --modisco GC_modisco_profile_v2.h5 \
      --arm pos_patterns \
      --pattern pattern_0 \
      --output umap_pattern0_subpatterns.png

Sequence-vs-attribution two-panel view:

    python seqlet_viz.py design1 \
      --modisco GC_modisco_profile_v2.h5 \
      --pattern pattern_0 \
      --output design1_pattern0.png

Interpretable center-similarity axes:

    python seqlet_viz.py design2 \
      --modisco GC_modisco_profile_v2.h5 \
      --pattern pattern_0 \
      --output design2_pattern0.png

Important scaling note
----------------------
The UMAP plots use a full pairwise Continuous Jaccard similarity matrix. This
is O(N^2) in memory and time. Use `--max-seqlets` for one-pattern plots and
`--per-pattern-cap` for atlas plots. The script refuses very large matrices
unless you raise `--max-pairwise-entries` deliberately.

Requirements
------------
    pip install h5py hdf5plugin numpy matplotlib umap-learn
"""

import argparse
from contextlib import contextmanager


@contextmanager
def open_h5(path):
    """Open an HDF5 file, with a clear message for compressed production files."""
    try:
        import hdf5plugin  # noqa: F401  -- imported before h5py when available
        hdf5plugin_available = True
    except ModuleNotFoundError:
        hdf5plugin_available = False
    import h5py

    try:
        with h5py.File(path, "r") as f:
            yield f
    except OSError as exc:
        msg = str(exc).lower()
        likely_filter_error = "filter" in msg or "plugin" in msg or "blosc" in msg
        if likely_filter_error and not hdf5plugin_available:
            raise OSError(
                f"Could not open {path!r}. This file may use compressed HDF5 "
                "filters that require 'hdf5plugin'. Install it, then rerun."
            ) from exc
        raise


def _natural_key(name):
    prefix, _, suffix = name.rpartition("_")
    if suffix.isdigit():
        return prefix, int(suffix)
    return name, -1


def _np():
    import numpy as np
    return np


def _plotting():
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    return cm, plt


# =====================================================================
# Continuous Jaccard similarity  (from kundajelab/tfmodisco)
# =====================================================================
def l1_norm(M):
    """L1-normalize each row so its absolute values sum to 1."""
    np = _np()
    return M / np.maximum(np.sum(np.abs(M), axis=1)[:, None], 1e-7)


def cj_row(a_row, mat):
    """Continuous Jaccard similarity of one row vs every row of `mat`.

    Intersection uses the sign-aware minimum of magnitudes; union uses the
    maximum of magnitudes. Result is in [-1, 1] (1 = identical shape and sign).
    """
    np = _np()
    union = np.sum(np.maximum(np.abs(a_row[None, :]), np.abs(mat)), axis=1)
    inter = np.sum(
        np.minimum(np.abs(a_row[None, :]), np.abs(mat))
        * np.sign(a_row[None, :])
        * np.sign(mat),
        axis=1,
    )
    return inter.astype(np.float32) / np.maximum(union, 1e-7)


def cj_simmat(M, dtype=None):
    """Full pairwise Continuous Jaccard similarity matrix for rows of M."""
    np = _np()
    if dtype is None:
        dtype = np.float32
    Mn = l1_norm(M.astype(dtype, copy=False))
    S = np.zeros((len(Mn), len(Mn)), dtype=dtype)
    for i in range(len(Mn)):
        S[i] = cj_row(Mn[i], Mn)
    return S


def _check_pairwise_size(n, max_pairwise_entries):
    entries = n * n
    if entries > max_pairwise_entries:
        approx_gb = entries * 4 / 1e9
        raise ValueError(
            f"{n:,} seqlets require {entries:,} pairwise entries "
            f"(~{approx_gb:.2f} GB per float32 matrix). Reduce --max-seqlets "
            "or --per-pattern-cap, or raise --max-pairwise-entries deliberately."
        )


def umap_from_sim(S, random_state=0):
    """Run UMAP on a precomputed similarity matrix."""
    np = _np()
    import umap

    D = (1.0 - S).astype(np.float32, copy=False)
    np.fill_diagonal(D, 0.0)
    return umap.UMAP(metric="precomputed", random_state=random_state).fit_transform(D)


# =====================================================================
# Data loading
# =====================================================================
def _normalize_seqlet_array(arr):
    """Return seqlet arrays as (n_seqlets, width, 4)."""
    np = _np()
    arr = np.asarray(arr)
    if arr.ndim != 3:
        raise ValueError(f"Expected a 3D seqlet array; got shape {arr.shape}")
    if arr.shape[-1] == 4:
        return arr
    if arr.shape[1] == 4:
        return arr.transpose(0, 2, 1)
    raise ValueError(f"Expected one seqlet dimension to be 4; got shape {arr.shape}")


def _flatten(arr):
    """(n_seqlets, width, 4) or (n_seqlets, 4, width) -> feature matrix."""
    np = _np()
    arr = _normalize_seqlet_array(arr)
    return arr.reshape(len(arr), -1).astype(np.float32, copy=False)


def _read_take(dataset, take):
    """Read only selected HDF5 rows. h5py fancy indexes must be sorted."""
    np = _np()
    take = np.asarray(take, dtype=np.int64)
    order = np.argsort(take)
    sorted_take = take[order]
    data = dataset[sorted_take]
    undo = np.argsort(order)
    return data[undo]


def load_pattern(modisco_h5, pattern, arm="pos_patterns", max_seqlets=3000, seed=0):
    """Load one pattern's seqlets, grouped by subpattern labels."""
    np = _np()
    rng = np.random.default_rng(seed)
    seqs, attrs, labels = [], [], []
    with open_h5(modisco_h5) as f:
        if arm not in f:
            raise KeyError(f"{arm!r} not found in {modisco_h5}")
        if pattern not in f[arm]:
            raise KeyError(f"{arm}/{pattern} not found in {modisco_h5}")

        parent = f[arm][pattern]
        subs = sorted([k for k in parent if k.startswith("subpattern_")], key=_natural_key)
        groups = subs if subs else ["__self__"]
        for name in groups:
            node = parent if name == "__self__" else parent[name]
            sl = node["seqlets"]
            n = sl["sequence"].shape[0]
            seqs.append(sl["sequence"][:])
            attrs.append(sl["contrib_scores"][:])
            label = "no_subpattern" if name == "__self__" else name
            labels += [label] * n

    labels = np.asarray(labels)
    seq_arr = np.vstack(seqs)
    attr_arr = np.vstack(attrs)
    if len(labels) > max_seqlets:
        idx = rng.choice(len(labels), max_seqlets, replace=False)
        seq_arr = seq_arr[idx]
        attr_arr = attr_arr[idx]
        labels = labels[idx]
        print(f"subsampled to {max_seqlets:,} seqlets")

    seq = _flatten(seq_arr)
    attr = _flatten(attr_arr)
    print(f"{arm}/{pattern}: {len(labels):,} seqlets, {len(set(labels))} labels")
    return {
        "seq": seq,
        "attr": attr,
        "labels": labels,
        "title": f"{arm.split('_')[0]}/{pattern}",
    }


def load_all_patterns(modisco_h5, per_pattern_cap=75, seed=0):
    """Pool seqlets across every pattern, with a cap per parent pattern."""
    np = _np()
    rng = np.random.default_rng(seed)
    seqs, attrs, labels = [], [], []
    with open_h5(modisco_h5) as f:
        for arm in ["pos_patterns", "neg_patterns"]:
            if arm not in f:
                continue
            for p in sorted(f[arm].keys(), key=_natural_key):
                sl = f[arm][p]["seqlets"]
                n = sl["sequence"].shape[0]
                take = rng.choice(n, min(n, per_pattern_cap), replace=False)
                seqs.append(_read_take(sl["sequence"], take))
                attrs.append(_read_take(sl["contrib_scores"], take))
                labels += [f"{arm.split('_')[0]}/{p}"] * len(take)

    if not labels:
        raise ValueError(f"No patterns found in {modisco_h5}")
    seq = _flatten(np.vstack(seqs))
    attr = _flatten(np.vstack(attrs))
    labels = np.asarray(labels)
    print(f"pooled {len(labels):,} seqlets across {len(set(labels))} patterns")
    return {"seq": seq, "attr": attr, "labels": labels, "title": "all patterns"}


# =====================================================================
# Plotting
# =====================================================================
def _color_map(labels):
    """Stable color per label; handles pos/pattern_0 and neg/pattern_0 ties."""
    np = _np()
    cm, _ = _plotting()

    def keyfn(label):
        prefix, _, suffix = str(label).rpartition("_")
        return (int(suffix), str(label)) if suffix.isdigit() else (10**9, str(label))

    uniq = sorted(set(labels), key=keyfn)
    palette = cm.tab20(np.linspace(0, 1, max(len(uniq), 1)))
    return uniq, {u: c for u, c in zip(uniq, palette)}


def _maybe_show(show):
    _, plt = _plotting()
    if show:
        plt.show()
    else:
        plt.close()


def plot_umap(
    bundle,
    view="attr",
    label_min=150,
    figsize=(11, 9),
    save=None,
    show=False,
    random_state=0,
    max_pairwise_entries=25_000_000,
):
    """One UMAP, colored by label."""
    _, plt = _plotting()
    M = bundle[view]
    _check_pairwise_size(len(M), max_pairwise_entries)
    emb = umap_from_sim(cj_simmat(M), random_state=random_state)
    labels = bundle["labels"]
    uniq, colors = _color_map(labels)

    fig, ax = plt.subplots(figsize=figsize)
    for u in uniq:
        m = labels == u
        ax.scatter(
            emb[m, 0],
            emb[m, 1],
            s=4,
            alpha=0.5,
            color=colors[u],
            label=u if m.sum() > label_min else None,
        )
    space = "attribution" if view == "attr" else "sequence"
    ax.set_title(f"{bundle['title']} — seqlet UMAP ({space} space)")
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.legend(markerscale=4, fontsize=7, ncol=2, loc="best")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print("saved", save)
    _maybe_show(show)


def plot_all_patterns_umap(bundle, **kwargs):
    return plot_umap(bundle, **kwargs)


def plot_subpattern_umap(bundle, **kwargs):
    return plot_umap(bundle, **kwargs)


def plot_design1(
    bundle,
    label_min=100,
    figsize=(16, 7),
    save=None,
    show=False,
    random_state=0,
    max_pairwise_entries=25_000_000,
):
    """Sequence-space UMAP and attribution-space UMAP, side by side."""
    _, plt = _plotting()
    _check_pairwise_size(len(bundle["seq"]), max_pairwise_entries)
    emb_seq = umap_from_sim(cj_simmat(bundle["seq"]), random_state=random_state)
    emb_attr = umap_from_sim(cj_simmat(bundle["attr"]), random_state=random_state)
    labels = bundle["labels"]
    uniq, colors = _color_map(labels)

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    for emb, ax, title in [
        (emb_seq, axes[0], "SEQUENCE space"),
        (emb_attr, axes[1], "ATTRIBUTION space"),
    ]:
        for u in uniq:
            m = labels == u
            ax.scatter(
                emb[m, 0],
                emb[m, 1],
                s=4,
                alpha=0.5,
                color=colors[u],
                label=u if m.sum() > label_min else None,
            )
        ax.set_title(title)
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
    axes[1].legend(markerscale=4, fontsize=7, loc="best")
    plt.suptitle(f"{bundle['title']} — same seqlets, two views")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print("saved", save)
    _maybe_show(show)


def plot_design2(bundle, label_min=100, figsize=(9, 9), save=None, show=False):
    """Plot sequence-center similarity vs attribution-center similarity."""
    np = _np()
    _, plt = _plotting()

    def sim_to_center(M):
        Mn = l1_norm(M)
        center = l1_norm(M.mean(axis=0, keepdims=True))[0]
        return cj_row(center, Mn)

    x = sim_to_center(bundle["seq"])
    y = sim_to_center(bundle["attr"])
    labels = bundle["labels"]
    uniq, colors = _color_map(labels)

    fig, ax = plt.subplots(figsize=figsize)
    for u in uniq:
        m = labels == u
        ax.scatter(
            x[m],
            y[m],
            s=6,
            alpha=0.5,
            color=colors[u],
            label=u if m.sum() > label_min else None,
        )
    ax.axvline(np.median(x), color="grey", ls="--", lw=0.8)
    ax.axhline(np.median(y), color="grey", ls="--", lw=0.8)
    ax.set_xlabel("sequence similarity to cluster center")
    ax.set_ylabel("attribution similarity to cluster center")
    ax.set_title(f"{bundle['title']} — each dot = one seqlet")
    ax.text(0.98, 0.98, "canonical\n(seq+ attr+)", ha="right", va="top",
            transform=ax.transAxes, fontsize=9)
    ax.text(0.02, 0.98, "different seq,\nsimilar attribution", ha="left", va="top",
            transform=ax.transAxes, fontsize=9)
    ax.text(0.98, 0.02, "similar seq,\natypical attribution", ha="right", va="bottom",
            transform=ax.transAxes, fontsize=9)
    ax.legend(markerscale=3, fontsize=7, loc="center left")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print("saved", save)
    _maybe_show(show)


def _add_common_plot_args(parser):
    parser.add_argument("--output", "--save", dest="output", required=True,
                        help="output PNG path")
    parser.add_argument("--label-min", type=int, default=100,
                        help="only show legend entries above this size")
    parser.add_argument("--seed", type=int, default=0,
                        help="random seed for subsampling and UMAP")
    parser.add_argument("--show", action="store_true",
                        help="also display an interactive plot window")


def _add_pattern_args(parser):
    parser.add_argument("--modisco", required=True, help="modisco-lite output .h5")
    parser.add_argument("--arm", default="pos_patterns",
                        choices=["pos_patterns", "neg_patterns"])
    parser.add_argument("--pattern", required=True, help='pattern name, e.g. "pattern_0"')
    parser.add_argument("--max-seqlets", type=int, default=3000,
                        help="subsample cap for one-pattern plots")


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    atlas = sub.add_parser("atlas", help="UMAP across all parent patterns")
    atlas.add_argument("--modisco", required=True, help="modisco-lite output .h5")
    atlas.add_argument("--view", default="attr", choices=["attr", "seq"])
    atlas.add_argument("--per-pattern-cap", type=int, default=75,
                       help="max seqlets sampled from each parent pattern")
    atlas.add_argument("--max-pairwise-entries", type=int, default=25_000_000)
    _add_common_plot_args(atlas)

    subpattern = sub.add_parser("subpattern", help="UMAP within one pattern")
    _add_pattern_args(subpattern)
    subpattern.add_argument("--view", default="attr", choices=["attr", "seq"])
    subpattern.add_argument("--max-pairwise-entries", type=int, default=25_000_000)
    _add_common_plot_args(subpattern)

    design1 = sub.add_parser("design1", help="sequence UMAP vs attribution UMAP")
    _add_pattern_args(design1)
    design1.add_argument("--max-pairwise-entries", type=int, default=25_000_000)
    _add_common_plot_args(design1)

    design2 = sub.add_parser("design2", help="center-similarity scatter")
    _add_pattern_args(design2)
    _add_common_plot_args(design2)

    return parser


def main(argv=None):
    args = build_arg_parser().parse_args(argv)

    if args.command == "atlas":
        bundle = load_all_patterns(args.modisco, args.per_pattern_cap, seed=args.seed)
        plot_all_patterns_umap(
            bundle,
            view=args.view,
            label_min=args.label_min,
            save=args.output,
            show=args.show,
            random_state=args.seed,
            max_pairwise_entries=args.max_pairwise_entries,
        )
    elif args.command == "subpattern":
        bundle = load_pattern(args.modisco, args.pattern, args.arm, args.max_seqlets, args.seed)
        plot_subpattern_umap(
            bundle,
            view=args.view,
            label_min=args.label_min,
            save=args.output,
            show=args.show,
            random_state=args.seed,
            max_pairwise_entries=args.max_pairwise_entries,
        )
    elif args.command == "design1":
        bundle = load_pattern(args.modisco, args.pattern, args.arm, args.max_seqlets, args.seed)
        plot_design1(
            bundle,
            label_min=args.label_min,
            save=args.output,
            show=args.show,
            random_state=args.seed,
            max_pairwise_entries=args.max_pairwise_entries,
        )
    elif args.command == "design2":
        bundle = load_pattern(args.modisco, args.pattern, args.arm, args.max_seqlets, args.seed)
        plot_design2(bundle, label_min=args.label_min, save=args.output, show=args.show)
    else:
        raise SystemExit(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
