"""Tests for Week 10 configuration and scientific parameters."""

from dataclasses import asdict, replace
import json
import subprocess
import sys

import pytest
import yaml

from tissue_sim.configuration import (
    SimulationConfig,
    read_simulation_config,
)
from tissue_sim.provenance import manifest, write_manifest
from tissue_sim.simulation import simulate

def test_config_precedence(tmp_path):
    path = tmp_path / "config.yaml"
    values = asdict(SimulationConfig())
    values["seed"] = 123
    values["marker_boost"] = 5.0
    path.write_text(yaml.safe_dump(values))

    config = read_simulation_config(path)
    overridden = read_simulation_config(path, overrides={"seed": 42})

    assert config.seed == 123
    assert overridden.seed == 42
    assert overridden.marker_boost == 5.0

@pytest.mark.parametrize(
    "problem, key",
    [
        ("unknown", "marker_bost"),
        ("missing", "cell_types"),
        ("missing", "n_cells"),
        ("invalid", "distribution_shape"),
    ],
)
def test_invalid_config(tmp_path, problem, key):
    path = tmp_path / "invalid.yaml"
    values = asdict(SimulationConfig())

    if problem == "unknown":
        values[key] = 10
    elif problem == "missing":
        del values[key]
    else:
        values[key] = "unsupported"

    path.write_text(yaml.safe_dump(values))

    with pytest.raises(ValueError) as error:
        read_simulation_config(path)

    assert key in str(error.value)
    assert str(path) in str(error.value)
    if problem == "unknown":
        assert "did you mean 'marker_boost'?" in str(error.value)


def test_sequencing_depth():
    config = SimulationConfig()

    a = simulate(config)
    b = simulate(replace(config, sequencing_depth=2.0))

    assert a.counts.equals(b.counts)
    assert a.transcripts.equals(b.transcripts)

    total_a = a.reference_counts.iloc[:, 1:].to_numpy().sum()
    total_b = b.reference_counts.iloc[:, 1:].to_numpy().sum()

    assert total_b > total_a

def test_cli_reproducibility(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    write_manifest(
        manifest(raw), tmp_path / "raw_manifest.sha256"
    )

    config = tmp_path / "config.yaml"
    small_config = SimulationConfig(
        field_width_um=32.0,
        field_height_um=32.0,
        n_cells=4,
        cell_radius_um=3.0,
        n_genes=20,
        n_marker_genes=3,
        reference_n_per_type=5,
    )
    config.write_text(yaml.safe_dump(asdict(small_config)))

    def execute(name, seed):
        out = tmp_path / name
        subprocess.run(
            [
                sys.executable, "run.py",
                "--config", str(config),
                "--data-dir", str(raw),
                "--out-dir", str(out),
                "--seed", str(seed),
            ],
            check=True,
        )
        return out

    a = execute("run_a", 42)
    b = execute("run_b", 42)
    c = execute("run_c", 43)

    assert (a / "counts.csv").read_bytes() == (
        b / "counts.csv"
    ).read_bytes()

    assert (a / "counts.csv").read_bytes() != (
        c / "counts.csv"
    ).read_bytes()

    metadata = json.loads((c / "run_metadata.json").read_text())

    assert metadata["resolved_config"]["seed"] == 43
    assert "--seed 43" in metadata["command_line"]
