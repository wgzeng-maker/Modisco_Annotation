from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from count_leftover_seqlets import estimate_clustering_entrants  # noqa: E402


def test_estimate_clustering_entrants_applies_cap_per_metacluster():
    assert estimate_clustering_entrants(418_064, 60_823, 100_000) == 160_823


def test_estimate_clustering_entrants_is_none_without_cap():
    assert estimate_clustering_entrants(10, 5, None) is None
