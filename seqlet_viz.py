"""
seqlet_viz.py
=============
Visualization toolkit for TF-MoDISco seqlets from ChromBPNet models.

Treats each seqlet like a "cell" in single-cell analysis: every seqlet is a
point, embedded and clustered by either its SEQUENCE (what base is present) or
its ATTRIBUTION profile (how much the model weights each base). The gap between
those two views is the object of interest -- "same sequence, different grammar".

Four visualizations
-------------------
1. plot_all_patterns_umap   : UMAP of every seqlet across all patterns,
                              colored by parent pattern (the "atlas" view).
2. plot_subpattern_umap     : UMAP of one pattern's seqlets, colored by the
                              subpatterns MoDISco found inside it.
3. plot_design1             : two UMAPs side by side for one pattern --
                              sequence space vs attribution space.
4. plot_design2             : one scatter for one pattern with interpretable
                              axes (sequence-sim vs attribution-sim to center).

Similarity metric
-----------------
Uses the Continuous Jaccard similarity from the original tfmodisco source
(sign-aware, magnitude-robust via L1 normalization). This is the same metric
MoDISco uses internally, so embeddings reflect how the algorithm actually
"sees" seqlet similarity -- more faithful than plain Euclidean distance.

Requirements
-----------
    pip install h5py hdf5plugin numpy umap-learn matplotlib
    (import name is `umap`, but the pip package is `umap-learn`)

Author: (your name)
"""

import hdf5plugin            # noqa: F401  -- MUST be imported before h5py
import h5py
import numpy as np
import umap
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# =====================================================================
# Continuous Jaccard similarity  (from kundajelab/tfmodisco)
# =====================================================================
def l1_norm(M):
    """L1-normalize each row so its absolute values sum to 1."""
    return M / np.maximum(np.sum(np.abs(M), axis=1)[:, None], 1e-7)


def cj_row(a_row, mat):
    """Continuous Jaccard similarity of one row vs every row of `mat`.

    Intersection uses the sign-aware minimum of magnitudes; union uses the
    maximum of magnitudes. Result is in [-1, 1] (1 = identical shape & sign).
    """
    union = np.sum(np.maximum(np.abs(a_row[None, :]), np.abs(mat)), axis=1)
    inter = np.sum(np.minimum(np.abs(a_row[None, :]), np.abs(mat))
                   * np.sign(a_row[None, :]) * np.sign(mat), axis=1)
    return inter.astype(float) / np.maximum(union, 1e-7)


def cj_simmat(M):
    """Full pairwise Continuous Jaccard similarity matrix for rows of M.

    WARNING: O(N^2) in memory and time. For >~5000 rows, subsample first
    (the loader functions below expose a `max_seqlets` / `per_pattern_cap`).
    """
    Mn = l1_norm(M)
    S = np.zeros((len(Mn), len(Mn)))
    for i in range(len(Mn)):
        S[i] = cj_row(Mn[i], Mn)
    return S


def umap_from_sim(S, random_state=0):
    """Run UMAP on a precomputed similarity matrix.

    UMAP needs distances, so we convert with (1 - similarity). The diagonal is
    forced to 0 (a point is identical to itself).
    """
    D = 1.0 - S
    np.fill_diagonal(D, 0.0)
    return umap.UMAP(metric="precomputed",
                     random_state=random_state).fit_transform(D)


# =====================================================================
# Data loading
# =====================================================================
def _flatten(arr):
    """(n_seqlets, width, 4) -> (n_seqlets, width*4) feature matrix."""
    return arr.reshape(len(arr), -1)


def load_pattern(modisco_h5, pattern, arm="pos_patterns",
                 max_seqlets=4000, seed=0):
    """Load one pattern's seqlets, grouped by its subpatterns.

    Parameters
    ----------
    modisco_h5 : str   path to a modisco-lite output .h5
    pattern    : str   e.g. "pattern_0"
    arm        : str   "pos_patterns" or "neg_patterns"
    max_seqlets: int   subsample cap (keeps the O(N^2) matrix affordable)

    Returns
    -------
    dict with keys: seq (N, W*4), attr (N, W*4), labels (N,), title (str).
    `labels` are the subpattern names, e.g. "subpattern_3".
    """
    rng = np.random.default_rng(seed)
    seqs, attrs, labels = [], [], []
    with h5py.File(modisco_h5, "r") as f:
        parent = f[arm][pattern]
        subs = sorted([k for k in parent if k.startswith("subpattern_")],
                      key=lambda s: int(s.split("_")[1]))
        if not subs:
            raise ValueError(
                f"{arm}/{pattern} has no subpatterns; "
                "use load_all_patterns or color by pattern instead.")
        for name in subs:
            sl = parent[name]["seqlets"]
            seqs.append(sl["sequence"][:])
            attrs.append(sl["contrib_scores"][:])
            labels += [name] * sl["sequence"].shape[0]

    seq = _flatten(np.vstack(seqs))
    attr = _flatten(np.vstack(attrs))
    labels = np.array(labels)

    if len(labels) > max_seqlets:
        idx = rng.choice(len(labels), max_seqlets, replace=False)
        seq, attr, labels = seq[idx], attr[idx], labels[idx]
        print(f"subsampled to {max_seqlets} seqlets")

    print(f"{arm}/{pattern}: {len(labels)} seqlets, "
          f"{len(set(labels))} subpatterns")
    return {"seq": seq, "attr": attr, "labels": labels,
            "title": f"{arm.split('_')[0]}/{pattern}"}


def load_all_patterns(modisco_h5, per_pattern_cap=300, seed=0):
    """Pool seqlets across EVERY pattern (both pos and neg).

    Subsamples each pattern to `per_pattern_cap` so large patterns don't
    dominate and the total stays small enough for an O(N^2) similarity matrix.

    Returns
    -------
    dict with keys: seq, attr, labels, title.
    `labels` are parent-pattern names, e.g. "pos/pattern_3".
    """
    rng = np.random.default_rng(seed)
    seqs, attrs, labels = [], [], []
    with h5py.File(modisco_h5, "r") as f:
        for arm in ["pos_patterns", "neg_patterns"]:
            if arm not in f:
                continue
            for p in sorted(f[arm].keys(), key=lambda s: int(s.split("_")[1])):
                sl = f[arm][p]["seqlets"]
                n = sl["sequence"].shape[0]
                take = rng.choice(n, min(n, per_pattern_cap), replace=False)
                take.sort()                       # h5 fancy-index needs sorted
                seqs.append(sl["sequence"][:][take])
                attrs.append(sl["contrib_scores"][:][take])
                labels += [f"{arm.split('_')[0]}/{p}"] * len(take)

    seq = _flatten(np.vstack(seqs))
    attr = _flatten(np.vstack(attrs))
    labels = np.array(labels)
    print(f"pooled {len(labels)} seqlets across {len(set(labels))} patterns")
    return {"seq": seq, "attr": attr, "labels": labels,
            "title": "all patterns"}


# =====================================================================
# Color helper
# =====================================================================
def _color_map(labels):
    """Stable color per label; sorts subpattern_/pattern_ names numerically."""
    def keyfn(s):
        # works for "subpattern_3" and "pos/pattern_3"
        return int(s.split("_")[-1])
    uniq = sorted(set(labels), key=keyfn)
    palette = cm.tab20(np.linspace(0, 1, max(len(uniq), 1)))
    return uniq, {u: c for u, c in zip(uniq, palette)}


# =====================================================================
# VISUALIZATION 1 + 2 : single-view UMAP (atlas or per-pattern)
# =====================================================================
def plot_umap(bundle, view="attr", label_min=150, figsize=(11, 9), save=None):
    """One UMAP, colored by label. Backbone for both the all-patterns atlas
    and the per-pattern subpattern view -- which one you get depends on the
    bundle you pass in (load_all_patterns vs load_pattern).

    view : "attr" (attribution space) or "seq" (sequence space).
    label_min : only label groups with more than this many points (keeps the
                legend readable when there are many groups).
    """
    M = bundle[view]
    emb = umap_from_sim(cj_simmat(M))
    labels = bundle["labels"]
    uniq, colors = _color_map(labels)

    fig, ax = plt.subplots(figsize=figsize)
    for u in uniq:
        m = labels == u
        ax.scatter(emb[m, 0], emb[m, 1], s=4, alpha=0.5, color=colors[u],
                   label=u if m.sum() > label_min else None)
    space = "attribution" if view == "attr" else "sequence"
    ax.set_title(f"{bundle['title']} — seqlet UMAP ({space} space)")
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.legend(markerscale=4, fontsize=7, ncol=2, loc="best")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print("saved", save)
    plt.show()


def plot_all_patterns_umap(bundle, **kwargs):
    """Convenience wrapper: the 'atlas' view (one dot per seqlet, colored by
    parent pattern). Pass a bundle from load_all_patterns()."""
    return plot_umap(bundle, **kwargs)


def plot_subpattern_umap(bundle, **kwargs):
    """Convenience wrapper: one pattern's seqlets colored by subpattern.
    Pass a bundle from load_pattern()."""
    return plot_umap(bundle, **kwargs)


# =====================================================================
# VISUALIZATION 3 : Design 1 -- two UMAPs side by side
# =====================================================================
def plot_design1(bundle, label_min=100, figsize=(16, 7), save=None):
    """Sequence-space UMAP and attribution-space UMAP, side by side, same
    points and colors. A group that is TIGHT in one panel but SPREAD in the
    other reveals a sequence/attribution mismatch.
    """
    emb_seq = umap_from_sim(cj_simmat(bundle["seq"]))
    emb_attr = umap_from_sim(cj_simmat(bundle["attr"]))
    labels = bundle["labels"]
    uniq, colors = _color_map(labels)

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    for emb, ax, title in [(emb_seq, axes[0], "SEQUENCE space"),
                           (emb_attr, axes[1], "ATTRIBUTION space")]:
        for u in uniq:
            m = labels == u
            ax.scatter(emb[m, 0], emb[m, 1], s=4, alpha=0.5, color=colors[u],
                       label=u if m.sum() > label_min else None)
        ax.set_title(title)
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
    axes[1].legend(markerscale=4, fontsize=7, loc="best")
    plt.suptitle(f"{bundle['title']} — same seqlets, two views")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print("saved", save)
    plt.show()


# =====================================================================
# VISUALIZATION 4 : Design 2 -- interpretable axes
# =====================================================================
def plot_design2(bundle, label_min=100, figsize=(9, 9), save=None):
    """Each dot = one seqlet. Axes are real quantities (unlike UMAP):
        x = sequence similarity to the cluster center
        y = attribution similarity to the cluster center
    The four corners are interpretable; median lines split them.
    """
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
        ax.scatter(x[m], y[m], s=6, alpha=0.5, color=colors[u],
                   label=u if m.sum() > label_min else None)
    ax.axvline(np.median(x), color="grey", ls="--", lw=0.8)
    ax.axhline(np.median(y), color="grey", ls="--", lw=0.8)
    ax.set_xlabel("sequence similarity to cluster center")
    ax.set_ylabel("attribution similarity to cluster center")
    ax.set_title(f"{bundle['title']} — each dot = one seqlet")
    ax.text(0.98, 0.98, "canonical\n(seq+ attr+)", ha="right", va="top",
            transform=ax.transAxes, fontsize=9)
    ax.text(0.02, 0.98, "different seq,\nsame grammar", ha="left", va="top",
            transform=ax.transAxes, fontsize=9)
    ax.text(0.98, 0.02, "same seq,\nmodel ignores it", ha="right", va="bottom",
            transform=ax.transAxes, fontsize=9)
    ax.legend(markerscale=3, fontsize=7, loc="center left")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print("saved", save)
    plt.show()


# =====================================================================
# Example usage (only runs if you execute this file directly)
# =====================================================================
if __name__ == "__main__":
    MODISCO = "/home/jupyter/GC_chrombpnet/GC_contribs/GC_modisco_profile_v2.h5"

    # 1. Atlas: all seqlets, all patterns
    atlas = load_all_patterns(MODISCO, per_pattern_cap=300)
    plot_all_patterns_umap(atlas, save="umap_all_patterns.png")

    # 2. One pattern, colored by subpattern
    p0 = load_pattern(MODISCO, "pattern_0")
    plot_subpattern_umap(p0, save="umap_pattern0_subpatterns.png")

    # 3. Design 1: sequence vs attribution, side by side
    plot_design1(p0, save="design1_pattern0.png")

    # 4. Design 2: interpretable axes
    plot_design2(p0, save="design2_pattern0.png")
