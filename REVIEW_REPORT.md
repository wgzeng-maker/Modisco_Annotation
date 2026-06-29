# Review report: first back-annotation milestone

## Scope

This pass focuses on `build_seqlet_annotation.py`, because the genome
back-annotation table is the foundation for downstream co-occurrence, overlap,
and cross-dataset analyses.

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

## Remaining review findings

These are intentionally left for later milestones:

- Package the scripts into a real CLI with subcommands.
- Make `seqlet_viz.py` laptop-safe by reducing default caps, documenting O(N²)
  memory, and avoiding unnecessary full-dataset reads.
- Make `upgrade_modisco_report` idempotent and rename the file without the
  upload suffix.
- Clarify `count_leftover_seqlets.py` documentation around the
  `max_seqlets_per_metacluster` cap; the current logic counts candidates and
  threshold survivors, not full clustering entrants.
- Add a fixture that proves subpattern seqlet identity-level partitioning, not
  only count-level partitioning.

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
python3 -m py_compile build_seqlet_annotation.py tests/test_build_seqlet_annotation.py
```

Full pytest execution is pending until the development dependencies are
installed in an environment with `h5py`, `pandas`, `pyfaidx`, and `pytest`.
