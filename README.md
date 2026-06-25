# Modisco_ClusterAnalysis

A downstream analysis pipeline built on Jacob Schreiber's [Modisco-lite](https://github.com/jmschrei/tfmodisco-lite).

Modisco-lite groups seqlets into clusters (called **Patterns**) and matches each one to a known motif from the public motif database JASPAR. It also splits each Pattern into sub-clusters (called **Subpatterns**). This repository analyzes and visualizes the heterogeneity within these Patterns and Subpatterns.

## Visualizing the clusters

Each point is one seqlet. Seqlets are embedded by similarity, so points that sit close together are alike.

**UMAP of all Patterns** — different motifs form separate islands:

<img width="1950" height="1650" alt="all_patterns_seqlet_umap" src="https://github.com/user-attachments/assets/9f1141e5-e951-42fb-9bfd-ba17eb70b737" />

**UMAP of the Subpatterns inside one Pattern** — reveals the finer structure Modisco found within a single cluster:

<img width="1500" height="1350" alt="pattern0_seqlet_umap" src="https://github.com/user-attachments/assets/9a4785c6-9254-430c-86cc-feffccdec3df" />


## What drives heterogeneity within a Pattern?

A Pattern can vary in two independent ways: its **sequence** (which bases are
present) and its **attribution** (how much the model weights each base). Each
seqlet therefore has two feature vectors:

- a **sequence vector** — one-hot encoded (1 for the base present, 0 otherwise)
- an **attribution vector** — the model's per-base contribution scores

Both plots below compare these two views. To compare seqlets we use
**Continuous Jaccard similarity**, the same metric Modisco uses internally — it
is sign-aware and robust to overall magnitude, so it captures the *shape* of a
seqlet rather than its raw scale.

### Plot type 1: two views, side by side

Two UMAP embeddings of the same seqlets — one built from the sequence vectors
(left), one from the attribution vectors (right). Each point is the same seqlet
in both panels, colored identically.

- **Left** — how seqlets group by sequence.
- **Right** — how seqlets group by attribution.

When a group is tight on one side but spread on the other, sequence and
attribution disagree — for example, *same sequence, different grammar*.

> Note: UMAP preserves which points are neighbors, but distances between
> separated blobs are not meaningful — read groupings, not gaps.

<img width="2400" height="1050" alt="design1_seq_vs_attr_umap" src="https://github.com/user-attachments/assets/60138751-0b99-4fce-a2c5-dfd8071defee" />

### Plot type 2: both views in one plot

This asks the same question on a single pair of axes that have direct meaning
(unlike UMAP axes, which are abstract). We first compute the **cluster center**
for each view — the average of all seqlet vectors in the Pattern, giving one
"typical sequence" and one "typical attribution profile." Then, for every
seqlet, we measure its Continuous Jaccard similarity to that center:

- **X axis** — sequence similarity to the cluster center (how typical this
  seqlet's *sequence* is).
- **Y axis** — attribution similarity to the cluster center (how typical this
  seqlet's *attribution* is).

Each seqlet becomes one point, and the corners are interpretable:

- **top-right** — typical in both: the canonical motif.
- **bottom-right** — typical sequence, atypical attribution: *same sequence,
  the model treats it differently*.
- **top-left** — atypical sequence, typical attribution: *different sequence,
  same grammar*.
- **bottom-left** — atypical in both.

The dashed lines mark the median of each axis, splitting the plot into these
four regions.

<img width="1350" height="1350" alt="design2_seq_vs_attr_axes" src="https://github.com/user-attachments/assets/42362391-d0bc-405a-9f8c-1e4117ff712b" />

## Limitation: most identified seqlets never enter clustering

A key limitation worth understanding before interpreting any MoDISco pattern set
— **the patterns are built from only a small fraction of the seqlets MoDISco
actually identifies.**

### What happens

MoDISco runs in stages, and seqlets are lost at each one:

1. **Seqlet identification** scans every position and flags windows whose
   importance passes a (deliberately lenient) FDR threshold.
2. **Metacluster cap** — each metacluster is capped at `max_seqlets_per_metacluster`
   (default **20,000**). Only the *strongest-attribution* seqlets enter
   clustering; the rest are discarded before clustering even begins.
3. **Clustering** drops more (noise filtering, clusters too small to survive).

This cap is not mentioned in the MoDISco README or technical note — it exists
only as a default in the function signature. It is a reasonable engineering
choice (clustering millions of seqlets is infeasible) and it keeps the
*strongest* seqlets, not random ones. But a user reading the documentation would
not know that most of their signal was set aside.

### Concrete numbers from one run

A ChromBPNet profile-head model on mouse cerebellar granule cells
(158,710 peak regions):

| Stage | Seqlets | Note |
|---|---|---|
| Candidate seqlets identified | ~3,635,000 | ~23 per 2 kb region |
| Survived pos/neg threshold | ~3,578,000 | threshold drops only ~1.5% |
| **Assigned to final patterns** | **62,767** | **~1.7% of candidates** |

The vast majority of the drop is **not** the importance threshold — it is the
20,000-per-metacluster cap. ~98% of identified seqlets pass the threshold, but
only the top 20,000 per sign enter clustering.

### Why this matters

The ~23 candidate seqlets per region are mostly weak-importance windows, not
strong motif instances, so much of what is dropped is genuinely low-signal.
**But not all of it.** Among the millions of sub-threshold-strength seqlets there
are very likely *real but weaker* regulatory elements — rarer motifs, weaker
binding sites, cell-type-specific grammar with modest attribution — that are
discarded simply for not being in the top 20,000 by strength.

In other words: **MoDISco patterns describe the strongest, most common motifs
well, but say nothing about the long tail of weaker signal.**

### Future direction

A robust method is needed to analyze the left-over seqlets — the ones that pass
the importance threshold but never enter clustering. Possible directions:

- cluster the discarded seqlets separately, in batches, rather than capping;
- match weak seqlets against the *strong* patterns already found, to rescue weak
  instances of known motifs;
- raise or remove the cap with a more scalable clustering backend.

Until then, any biological claim from a MoDISco run should be read as a statement
about the **dominant** motifs, explicitly not about the full set of important
sequence in the genome.

*(The script `count_leftover_seqlets.py` in this repo reproduces the funnel
numbers above for any modisco-lite run.)*

