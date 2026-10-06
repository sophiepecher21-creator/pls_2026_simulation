"""Check raw-data fingerprints and the safeguard against changed inputs."""

import pytest

from tissue_sim.provenance import (
    manifest,
    require_verified,
    verify,
    write_manifest,
)


def small_raw_directory(tmp_path):
    """Create one input file and save its original fingerprint."""
    root = tmp_path / "raw"
    root.mkdir()
    (root / "counts.csv").write_text("square_id,g1\n0,7\n")

    manifest_path = tmp_path / "raw_manifest.sha256"
    write_manifest(manifest(root), manifest_path)

    return root, manifest_path


def test_untouched_directory(tmp_path):
    root, manifest_path = small_raw_directory(tmp_path)

    assert verify(root, manifest_path) == {
        "changed": set(),
        "missing": set(),
        "added": set(),
    }
    require_verified(root, manifest_path)


def test_modified_file(tmp_path):
    root, manifest_path = small_raw_directory(tmp_path)
    (root / "counts.csv").write_text("square_id,g1\n0,999\n")

    assert verify(root, manifest_path)["changed"] == {"counts.csv"}


def test_missing_and_added_files(tmp_path):
    root, manifest_path = small_raw_directory(tmp_path)
    (root / "counts.csv").unlink()
    (root / "notes.txt").write_text("An unexpected file\n")

    assert verify(root, manifest_path) == {
        "changed": set(),
        "missing": {"counts.csv"},
        "added": {"notes.txt"},
    }


def test_changed_input_stops_processing(tmp_path):
    root, manifest_path = small_raw_directory(tmp_path)
    (root / "counts.csv").write_text("square_id,g1\n0,999\n")

    with pytest.raises(ValueError, match="Raw-data verification failed"):
        require_verified(root, manifest_path)
