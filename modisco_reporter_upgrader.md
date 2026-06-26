# Modisco_Annotation

Post-processing utilities for [TF-MoDISco](https://github.com/jmschrei/tfmodisco-lite)
motif discovery output. The main tool, `upgrade_modisco_report.py`, takes an existing
modisco HTML report and **appends annotation columns to it** without touching any of the
original content.

## What it adds

Starting from your standard `motifs.html`, the script keeps every original column
(CWM forward/reverse logos, `match0/1/2` IDs and their JASPAR logos, q-values) and
appends five new columns at the right edge:

| New column | Source | Notes |
|---|---|---|
| **PFM forward** | `.h5` `sequence` array | Information-content sequence logo (base-frequency, not contribution) |
| **PFM reverse** | `.h5` `sequence` array | Reverse-complement of the PFM |
| **matched TFs** | report match IDs + JASPAR `.meme` | Resolves `MA0671.2 → NFIX`; **reuses** the IDs already in the report, does not re-run TOMTOM |
| **#subpatterns** | `.h5` | Number of subpatterns in the cluster |
| **functional annotation** | `--notes` CSV (optional) | Blank by default; hook for downstream ontology / GO results |

The PFM logos complement the CWM logos already in the report: the CWM shows *what the
model weights*, the PFM shows *what bases are actually present*. Comparing the two is a
more grounded check than reading the CWM alone.

## Install

```bash
git clone https://github.com/wgzeng-maker/Modisco_Annotation.git
cd Modisco_Annotation
pip install -r requirements.txt
```

## Usage

```bash
python upgrade_modisco_report.py \
    --report  profile_report/motifs.html \
    --h5      modisco_results.h5 \
    --meme    JASPAR2024_CORE_non-redundant_pfms_meme.txt
```

Arguments:

| Flag | Required | Description |
|---|---|---|
| `--report` | yes | Existing modisco report HTML |
| `--h5` | yes | modisco results `.h5` (provides PFM arrays and subpatterns) |
| `--meme` | yes | JASPAR `.meme` file used for the original report (for ID → TF name) |
| `--out` | no | Output path (default: `<report stem>_upgraded.html`) |
| `--trim-frac` | no | Trim PFM logo to columns whose CWM signal ≥ this fraction of the max (default: `0.30`) |
| `--notes` | no | CSV with `pattern,annotation` columns to fill the functional-annotation column |

On startup the script prints a short `[verify]` block (parsed headers, match-ID
columns, an example ID → TF resolution). Glance at it to confirm the report and `.meme`
were read correctly before trusting the output.

## Adding functional annotations later

The `functional annotation` column is wired to an optional CSV so you don't edit code
when ontology results arrive:

```csv
pattern,annotation
pos_patterns.pattern_0,"enriched near neuronal differentiation genes"
pos_patterns.pattern_6,"synaptic signaling"
```

```bash
python upgrade_modisco_report.py --report ... --h5 ... --meme ... --notes notes.csv
```

## Note on image portability

The **new** PFM logos are base64-embedded, so they render in any browser. The
**original** report's logos are referenced by file path (e.g. `/work/.../*.png`), so the
original-logo columns only render in the environment that produced them (the paths must
resolve). The annotation columns are unaffected by this.
