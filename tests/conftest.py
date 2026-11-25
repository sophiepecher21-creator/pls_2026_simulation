"""Shared fixtures for small simulation runs."""

from dataclasses import asdict, replace

import pytest
import yaml

from tissue_sim.configuration import read_simulation_config
from tissue_sim.provenance import manifest, write_manifest


@pytest.fixture
def small_config():
    """Return a fast, low-density simulation configuration."""
    return replace(
        read_simulation_config(),
        n_cells=4,
        field_width_um=16.0,
        field_height_um=16.0,
        cell_radius_um=1.0,
        n_genes=6,
        n_marker_genes=2,
        mean_transcripts_per_type=(20.0, 40.0),
    )


@pytest.fixture
def run_inputs(tmp_path, small_config):
    """Write a tiny configuration and a checksummed raw input for the runner."""
    raw_dir = tmp_path / "data" / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "example.csv").write_text("value\n7\n", encoding="utf-8")
    write_manifest(manifest(raw_dir), raw_dir.parent / "raw_manifest.sha256")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(asdict(small_config)), encoding="utf-8")
    return config_path, raw_dir
