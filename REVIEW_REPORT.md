# Review report: back-annotation, seqlet counting, and visualization

## Scope

The first pass focused on `build_seqlet_annotation.py`, because the genome
back-annotation table is the foundation for downstream co-occurrence, overlap,
and cross-dataset analyses. This follow-up also reviewed the notebooks that
produced the README plots/counts and wrapped the seqlet-counting and
visualization code into reusable command-line scripts.

## Findings addressed

1. `hdf5plugin` was imported unconditionally. This prevented uncompressed test
   fixtures from running on a lightweight laptop environment. The import is now
   optional, with a clearer error if a production compressed HDF5 file actually
   requires the plugin.
2. Subpattern iteration is preserved. This is the correct behavior for
   modisco-lite outputs where a parent pattern's seqlets are partitioned across
   subpatterns. To guard against silent omissions, the script now checks that
   the sum of subpattern `n_seqlets` equals the parent pattern `n_seqlets`.
3. Coordinate verification previously checked one seqlet. `--genome` now checks
   a sample across arms, patterns, subpatterns, and both strand orientations when
   available.
4. Basic tests now cover:
   - narrowPeak summit-column validation;
   - BED-style coordinate reconstruction;
   - subpattern count validation;
   - FASTA verification for forward and reverse-complement seqlets.
5. `seqlet_viz.py` now has a real CLI with `atlas`, `subpattern`, `design1`,
   and `design2` commands. Hard-coded notebook paths were removed.
6. Visualization defaults are more laptop-safe: the all-pattern atlas now uses a
   lower per-pattern cap, similarity matrices are float32, and the script
   refuses very large O(N²) pairwise matrices unless explicitly overridden.
7. `count_leftover_seqlets.py` now distinguishes measured quantities from the
   estimated metacluster-cap input. It reports candidate seqlets, pos/neg
   threshold survivors, optional estimated clustering entrants after `-n`, and
   final pattern-assigned seqlets.
8. The README and requirements now mention visualization dependencies and the
   special `modiscolite` runtime requirement for the counting script.
9. Heavy dependencies in the visualization/counting scripts are loaded lazily,
   so `--help` works even before the scientific Python environment is installed.

## Remaining review findings

These are intentionally left for later milestones:

- Package the scripts into an installable Python package with a top-level CLI.
- Make `upgrade_modisco_report` idempotent and rename the file without the
  upload suffix.
- Add a fixture that proves subpattern seqlet identity-level partitioning, not
  only count-level partitioning.
- Add integration tests for `seqlet_viz.py` once a small HDF5 fixture and the
  visualization dependencies are available in the test environment.

## Development notes

Install development dependencies with:

```bash
python -m pip install -r requirements-dev.txt
```

Run tests with:

```bash
python -m pytest
```

In this Codex session, syntax validation passed:

```bash
python3 -m py_compile \
  build_seqlet_annotation.py \
  count_leftover_seqlets.py \
  seqlet_viz.py \
  tests/test_build_seqlet_annotation.py \
  tests/test_count_leftover_seqlets.py
```

These lightweight checks also passed without the scientific dependencies:

```bash
python3 count_leftover_seqlets.py --help
python3 seqlet_viz.py --help
python3 -c 'from count_leftover_seqlets import estimate_clustering_entrants; assert estimate_clustering_entrants(418064,60823,100000)==160823'
```

Full pytest execution is pending until the development dependencies are
installed in an environment with `h5py`, `pandas`, `pyfaidx`, `umap-learn`,
`matplotlib`, and `pytest`.
