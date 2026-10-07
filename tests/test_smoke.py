"""Check that the toy runs, reproduces its output, and rejects changed input."""

import json
import shutil

import pandas as pd
import pytest
import yaml

from tissue_sim.provenance import manifest, write_manifest
from tissue_sim.simulation import run


def prepare_run(tmp_path):
    """Create temporary inputs and a configuration for the toy."""
    raw_dir = tmp_path / "data" / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "example.csv").write_text("value\n7\n")

    manifest_path = raw_dir.parent / "raw_manifest.sha256"
    write_manifest(manifest(raw_dir), manifest_path)

    config = {
        "seed": 20260101,
        "n_cells": 20,
        "field_width_um": 320.0,
        "field_height_um": 320.0,
        "cell_types": ["type_a", "type_b"],
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    return config_path, raw_dir, tmp_path / "results"


def test_program_runs_end_to_end(tmp_path):
    config_path, raw_dir, out_dir = prepare_run(tmp_path)

    run(config_path, raw_dir, out_dir)

    assert (out_dir / "toy_cells.csv").is_file()
    assert (out_dir / "run_metadata.json").is_file()

    cells = pd.read_csv(out_dir / "toy_cells.csv")
    assert len(cells) == 20
    assert cells["cell_id"].is_unique
    assert cells["x_um"].between(0, 320, inclusive="left").all()
    assert cells["y_um"].between(0, 320, inclusive="left").all()

    metadata = json.loads((out_dir / "run_metadata.json").read_text())
    assert set(metadata) == {
        "git_commit",
        "git_clean",
        "config",
        "seed",
        "input_manifest_sha256",
        "created",
        "python",
    }
    assert metadata["seed"] == 20260101


def test_delete_and_rebuild_reproduces_cells(tmp_path):
    config_path, raw_dir, out_dir = prepare_run(tmp_path)

    run(config_path, raw_dir, out_dir)
    original = (out_dir / "toy_cells.csv").read_bytes()

    shutil.rmtree(out_dir)
    run(config_path, raw_dir, out_dir)

    assert (out_dir / "toy_cells.csv").read_bytes() == original


def test_changed_input_prevents_output(tmp_path):
    config_path, raw_dir, out_dir = prepare_run(tmp_path)
    (raw_dir / "example.csv").write_text("value\n999\n")

    with pytest.raises(ValueError, match="Raw-data verification failed"):
        run(config_path, raw_dir, out_dir)

    assert not out_dir.exists()
