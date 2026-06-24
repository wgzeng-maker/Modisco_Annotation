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
