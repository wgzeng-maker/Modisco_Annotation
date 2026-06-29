# Modisco_ClusterAnalysis

A downstream analysis pipeline built on Jacob Schreiber's [Modisco-lite](https://github.com/jmschrei/tfmodisco-lite).

Modisco-lite groups seqlets into clusters (called **Patterns**) and matches each one to a known motif from the public motif database JASPAR. It also splits each Pattern into sub-clusters (called **Subpatterns**). This repository analyzes and visualizes the heterogeneity within these Patterns and Subpatterns.

## Current development status

This repository is still in research-toolkit form. The first tested command is
the back-annotation script, which maps clustered MoDISco seqlets back to genome
coordinates:

```bash
python build_seqlet_annotation.py \
  --modisco GC_modisco_profile_v2.h5 \
  --bed GC_mm10.interpreted_regions.bed \
  --output GC_profile_seqlet_annotation \
  --window 1000 \
  --input-len 2114 \
  --genome mm10.fa
```

Coordinates are BED-style: 0-based, half-open intervals. When `--genome` is
provided, the script samples seqlets across patterns/subpatterns/strands and
checks the stored seqlet sequence against the reference FASTA before writing the
annotation table.

For local development:

```bash
python -m pip install -r requirements-dev.txt
python -m pytest
```

`hdf5plugin` is required for many production MoDISco HDF5 files because they may
use compressed HDF5 filters. The small synthetic tests do not require compressed
HDF5.

The seqlet-counting script additionally requires `modiscolite`, so run it in
the same ChromBPNet/MoDISco environment or Docker image used to generate the
MoDISco output.

### Back-annotation

```bash
python build_seqlet_annotation.py \
  --modisco GC_modisco_profile_v2.h5 \
  --bed GC_mm10.interpreted_regions.bed \
  --output GC_profile_seqlet_annotation \
  --window 1000 \
  --input-len 2114 \
  --genome mm10.fa
```

### Seqlet-count funnel

```bash
python count_leftover_seqlets.py \
  --scores GC_mm10.profile_scores.h5 \
  --modisco GC_modisco_profile_v2.h5 \
  --window 1000 \
  -n 100000
```

This reports candidate seqlets, positive/negative threshold survivors, an
estimated upper bound on clustering input after the `-n` cap, and seqlets saved
in final patterns. It does not rerun clustering.

### Seqlet visualization

```bash
# Atlas across all parent patterns. Keep the cap low on laptops.
python seqlet_viz.py atlas \
  --modisco GC_modisco_profile_v2.h5 \
  --per-pattern-cap 75 \
  --output umap_all_patterns.png

# Subpatterns inside one parent pattern.
python seqlet_viz.py subpattern \
  --modisco GC_modisco_profile_v2.h5 \
  --arm pos_patterns \
  --pattern pattern_0 \
  --max-seqlets 3000 \
  --output umap_pattern0_subpatterns.png

# Sequence-vs-attribution UMAP panels.
python seqlet_viz.py design1 \
  --modisco GC_modisco_profile_v2.h5 \
  --pattern pattern_0 \
  --output design1_pattern0.png

# Interpretable center-similarity scatter.
python seqlet_viz.py design2 \
  --modisco GC_modisco_profile_v2.h5 \
  --pattern pattern_0 \
  --output design2_pattern0.png
```

The UMAP visualizations use a pairwise Continuous Jaccard matrix, so memory and
time scale as O(N²). Reduce `--max-seqlets` or `--per-pattern-cap` before
running large plots locally.

### Upgraded MoDISco report PDF

The original MoDISco HTML report can reference images by relative paths, so
images may disappear if the HTML is moved away from its report image directory.
Generate a standalone PDF while the HTML is still next to its image files:

```bash
python upgrade_modisco_report.py \
  --report profile_report/motifs.html \
  --h5 GC_modisco_profile_v2.h5 \
  --meme JASPAR2024_CORE_non-redundant_pfms_meme.txt \
  --out-pdf motifs_upgraded.pdf
```

The PDF embeds the original report images that can be resolved from the HTML,
plus regenerated forward/reverse PFM logos from the MoDISco H5 file.

To also write self-contained HTML with base64-embedded images:

```bash
python upgrade_modisco_report.py \
  --report profile_report/motifs.html \
  --h5 GC_modisco_profile_v2.h5 \
  --meme JASPAR2024_CORE_non-redundant_pfms_meme.txt \
  --out-html motifs_upgraded.html \
  --embed-html-images
```

### GCP back-annotation test

On the GCP VM, run back-annotation with the exact files used for the MoDISco
run. The key check is that `--genome` verification passes before outputs are
trusted:

```bash
cd /home/jupyter/GC_chrombpnet/GC_contribs

python /path/to/Modisco_Annotation/build_seqlet_annotation.py \
  --modisco GC_modisco_profile_v2.h5 \
  --bed GC_mm10.interpreted_regions.bed \
  --output GC_profile_seqlet_annotation \
  --window 1000 \
  --input-len 2114 \
  --genome /home/jupyter/GC_chrombpnet/data/downloads/mm10.fa
```

Expected behavior:

- the script reports `MATCH=True` for sampled FASTA checks;
- it writes `GC_profile_seqlet_annotation.bed`;
- it writes `GC_profile_seqlet_annotation.csv`;
- the annotated seqlet count should match the final clustered seqlet total.

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
— **the patterns are built from only a fraction of the seqlets MoDISco actually
identifies.**

### What happens

MoDISco runs in stages, and seqlets are lost at each one:

1. **Seqlet identification** scans every position in the motif-discovery window
   and flags windows whose importance passes a (deliberately lenient) FDR
   threshold (`target_fdr=0.2` by default).
2. **Metacluster cap** — each metacluster is capped at `max_seqlets_per_metacluster`
   (the `-n` flag; library default **20,000**). Only the *strongest-attribution*
   seqlets enter clustering; the rest are discarded before clustering begins.
3. **Clustering** drops more (noise filtering, clusters too small to survive).

The cap is not discussed in the MoDISco README or technical note — it exists
only as a settable parameter. It is a reasonable engineering choice (clustering
millions of seqlets is infeasible) and it keeps the *strongest* seqlets, not
random ones. But a user reading the documentation would not know that most of
their signal was set aside.

### Concrete numbers from one run

A ChromBPNet profile-head model on mouse cerebellar granule cells
(158,710 peak regions), run with `-w 1000 -n 100000`:

| Stage | Seqlets | Note |
|---|---|---|
| Candidate seqlets identified | 478,887 | ~3.0 per 1 kb window |
| Survived pos/neg threshold | 478,887 | threshold dropped **0** here |
| **Assigned to final patterns** | **62,767** | **~13% of candidates** |

In this run the importance threshold removed nothing — every candidate was
confidently positive or negative. The entire gap from 478,887 down to 62,767 is
the **100,000-per-metacluster cap plus clustering**: there were 418,064 positive
candidates, but only the top 100,000 by strength entered clustering, and
clustering reduced those further.

> Note: the candidate count is sensitive to the `-w` (window) parameter. This
> run used `-w 1000`, so extraction scanned only the central 1 kb of each 2114 bp
> input. Counting on the full 2114 bp window instead gives ~3.6M candidates — an
> over-count that does not reflect what the real run saw. Always match `-w`.

## The complete seqlet funnel

MoDISco discards seqlets at every stage. Here is the full accounting for one
ChromBPNet profile-head run (mouse cerebellar granule cells, 158,710 peaks,
`-w 1000 -n 100000`):

| Stage | Seqlets | What happens |
|---|---|---|
| Candidates identified | 478,887 | windows passing the lenient FDR (~3 per 1 kb region) |
| Survived pos/neg threshold | 478,887 | split by sign; here the threshold dropped **0** |
| → positive metacluster | 418,064 | |
| → negative metacluster | 60,823 | |
| Entered clustering | ~160,823 | top **100,000** positive (cap) + all 60,823 negative |
| **Assigned to final patterns** | **62,767** | survived clustering into 53 patterns |

### Where seqlets are lost

Two separate bottlenecks, often confused:

**1. The metacluster cap (the largest loss).**
Each metacluster is capped at `-n` (here 100,000). The positive metacluster had
418,064 seqlets but only the **top 100,000 by attribution strength** entered
clustering — the other ~318,000 were discarded before clustering began. The
negative metacluster (60,823) was under the cap, so all of it entered.

**2. Clustering-stage quality filters (~98,000 dropped).**
Of the ~160,823 seqlets that entered clustering, only 62,767 survived. The
~98,000 difference is removed by three quality filters inside MoDISco:

- **noise filtering** — seqlets whose coarse and fine similarity disagree
  (look like noise, no consistent neighbors) are discarded;
- **small-cluster disbanding** — clusters below the minimum size are broken up,
  and their seqlets dropped unless similar enough to a surviving pattern;
- **low-information-content filtering** — whole patterns with too little
  sequence information are removed.

These three are *quality* filters — the seqlets entered clustering but never
formed or joined a clean motif. Their individual counts are **not saved** by
MoDISco, so only the combined ~98,000 is measurable; the per-cause breakdown
would require re-running clustering with added logging.

### Bottom line

Of 478,887 confidently-important candidate seqlets, only 62,767 (~13%) end up in
final patterns. The loss is dominated by the metacluster cap (a tractability
limit, not a biological one) and clustering-stage quality filters. **MoDISco
patterns therefore describe the strongest, cleanest, most common motifs — not
the full set of important sequence in the genome.**

### Why this matters

The ~3 candidate seqlets per region are a mix of real motif instances and
weaker-importance windows. Much of what is dropped past the cap is genuinely
lower-signal — **but not all of it.** Among the hundreds of thousands of
sub-cap-strength seqlets there are very likely *real but weaker* regulatory
elements — rarer motifs, weaker binding sites, cell-type-specific grammar with
modest attribution — discarded simply for not being in the top 100,000 by
strength.

In other words: **MoDISco patterns describe the strongest, most common motifs
well, but say little about the long tail of weaker signal.**

### Future direction

A robust method is needed to analyze the left-over seqlets — those that pass the
importance threshold but never enter clustering. Possible directions:

- cluster the discarded seqlets separately, in batches, rather than capping;
- match weak seqlets against the *strong* patterns already found, to rescue weak
  instances of known motifs;
- raise or remove the cap with a more scalable clustering backend.

Until then, any biological claim from a MoDISco run should be read as a statement
about the **dominant** motifs, explicitly not about the full set of important
sequence in the genome.

*(The script `count_leftover_seqlets.py` reproduces this funnel for any
modisco-lite run — pass the matching `--window`.)*
