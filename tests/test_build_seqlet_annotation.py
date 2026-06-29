from pathlib import Path
import sys

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("pandas")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from build_seqlet_annotation import (  # noqa: E402
    build_annotation,
    load_peak_windows,
    verify_frame,
)


def _seq_to_onehot(seq):
    bases = "ACGT"
    arr = np.zeros((len(seq), 4), dtype="float32")
    for i, base in enumerate(seq):
        arr[i, bases.index(base)] = 1.0
    return arr


def _write_bed(path):
    path.write_text(
        "\n".join(
            [
                "chr1\t1000\t1500\tpeak_0\t100\t.\t5.0\t10.0\t8.0\t50",
                "chr2\t2000\t2500\tpeak_1\t100\t.\t5.0\t10.0\t8.0\t25",
            ]
        )
        + "\n"
    )


def _write_seqlets(group, example_idx, start, is_revcomp, width=10):
    example_idx = np.asarray(example_idx, dtype="int64")
    start = np.asarray(start, dtype="int64")
    is_revcomp = np.asarray(is_revcomp, dtype=bool)
    n = len(example_idx)
    sequence = np.zeros((n, width, 4), dtype="float32")
    sequence[:, :, 0] = 1.0
    group.create_dataset("sequence", data=sequence)
    group.create_dataset("contrib_scores", data=np.zeros_like(sequence))
    group.create_dataset("hypothetical_contribs", data=np.zeros_like(sequence))
    group.create_dataset("example_idx", data=example_idx)
    group.create_dataset("start", data=start)
    group.create_dataset("end", data=start + width)
    group.create_dataset("is_revcomp", data=is_revcomp)
    group.create_dataset("n_seqlets", data=np.array([n], dtype="int64"))


def _write_annotation_h5(path, mismatched_parent_count=False):
    with h5py.File(path, "w") as f:
        pos = f.create_group("pos_patterns")
        p0 = pos.create_group("pattern_0")
        parent_n = 4 if mismatched_parent_count else 3
        _write_seqlets(p0.create_group("seqlets"), [0] * parent_n, [0] * parent_n, [False] * parent_n)

        sp0 = p0.create_group("subpattern_0")
        _write_seqlets(sp0.create_group("seqlets"), [0, 1], [2, 5], [False, True])

        sp1 = p0.create_group("subpattern_1")
        _write_seqlets(sp1.create_group("seqlets"), [0], [15], [False])

        neg = f.create_group("neg_patterns")
        n0 = neg.create_group("pattern_0")
        _write_seqlets(n0.create_group("seqlets"), [1], [7], [False])


def test_build_annotation_maps_subpattern_seqlets_to_bed_coordinates(tmp_path):
    bed = tmp_path / "regions.bed"
    h5 = tmp_path / "modisco.h5"
    _write_bed(bed)
    _write_annotation_h5(h5)

    chrom, window_start = load_peak_windows(bed, input_len=100)
    ann = build_annotation(h5, chrom, window_start, crop_offset=40)

    assert len(ann) == 4
    assert ann["pattern"].tolist() == [
        "pos/pattern_0",
        "pos/pattern_0",
        "pos/pattern_0",
        "neg/pattern_0",
    ]
    assert ann["subpattern"].tolist() == [
        "subpattern_0",
        "subpattern_0",
        "subpattern_1",
        "",
    ]
    assert ann["chrom"].tolist() == ["chr1", "chr2", "chr1", "chr2"]
    assert ann["start"].tolist() == [1042, 2020, 1055, 2022]
    assert ann["end"].tolist() == [1052, 2030, 1065, 2032]
    assert ann["strand"].tolist() == ["+", "-", "+", "+"]


def test_build_annotation_fails_if_subpatterns_do_not_partition_parent_count(tmp_path):
    bed = tmp_path / "regions.bed"
    h5 = tmp_path / "modisco.h5"
    _write_bed(bed)
    _write_annotation_h5(h5, mismatched_parent_count=True)

    chrom, window_start = load_peak_windows(bed, input_len=100)
    with pytest.raises(ValueError, match="subpatterns contain 3 seqlets.*parent pattern reports 4"):
        build_annotation(h5, chrom, window_start, crop_offset=40)


def test_load_peak_windows_requires_narrowpeak_summit_column(tmp_path):
    bed = tmp_path / "bad.bed"
    bed.write_text("chr1\t0\t100\n")

    with pytest.raises(ValueError, match="10-column narrowPeak"):
        load_peak_windows(bed, input_len=100)


def test_verify_frame_checks_forward_and_reverse_complement_seqlets(tmp_path):
    pytest.importorskip("pyfaidx")
    bed = tmp_path / "regions.bed"
    h5 = tmp_path / "modisco.h5"
    fasta = tmp_path / "genome.fa"

    bed.write_text("chrTest\t0\t100\tpeak_0\t100\t.\t5.0\t10.0\t8.0\t10\n")

    genome = list("N" * 40)
    genome[7:11] = list("ACGT")
    genome[17:21] = list("AAGT")
    fasta.write_text(">chrTest\n" + "".join(genome) + "\n")

    with h5py.File(h5, "w") as f:
        p0 = f.create_group("pos_patterns").create_group("pattern_0")
        _write_seqlets(p0.create_group("seqlets"), [0, 0], [2, 12], [False, True], width=4)
        sp0 = p0.create_group("subpattern_0")
        sl = sp0.create_group("seqlets")
        sl.create_dataset("sequence", data=np.stack([_seq_to_onehot("ACGT"), _seq_to_onehot("ACTT")]))
        sl.create_dataset("contrib_scores", data=np.zeros((2, 4, 4), dtype="float32"))
        sl.create_dataset("hypothetical_contribs", data=np.zeros((2, 4, 4), dtype="float32"))
        sl.create_dataset("example_idx", data=np.array([0, 0], dtype="int64"))
        sl.create_dataset("start", data=np.array([2, 12], dtype="int64"))
        sl.create_dataset("end", data=np.array([6, 16], dtype="int64"))
        sl.create_dataset("is_revcomp", data=np.array([False, True]))
        sl.create_dataset("n_seqlets", data=np.array([2], dtype="int64"))

    chrom, window_start = load_peak_windows(bed, input_len=20)
    assert verify_frame(h5, fasta, chrom, window_start, crop_offset=5, max_checks=4) is True
