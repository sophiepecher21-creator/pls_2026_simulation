"""Tests for Week 10 configuration and scientific parameters."""

from dataclasses import asdict, replace
import json
import subprocess
import sys
from unittest.mock import Mock

import pytest
import yaml

from tissue_sim.configuration import read_simulation_config
from tissue_sim.place_cells import place_cells
from tissue_sim.provenance import manifest, write_manifest
from tissue_sim.simulation import simulate

def test_config_precedence(tmp_path):
    path = tmp_path / "config.yaml"
    defaults = read_simulation_config()
    values = asdict(defaults)
    values["seed"] = 123
    del values["n_cells"]
    path.write_text(yaml.safe_dump(values))

    config = read_simulation_config(path)
    overridden = read_simulation_config(path, overrides={"seed": 42})

    assert config.seed == 123
    assert overridden.seed == 42
    assert overridden.n_cells == defaults.n_cells

@pytest.mark.parametrize(
    "problem, key, reason",
    [
        ("unknown", "marker_bost", "Unknown configuration keys"),
        ("missing", "cell_types", "Missing required configuration keys"),
        ("invalid", "distribution_shape", "Unknown distribution_shape"),
    ],
)
def test_invalid_config(tmp_path, small_config, problem, key, reason):
    path = tmp_path / "invalid.yaml"
    values = asdict(small_config)

    if problem == "unknown":
        values[key] = 10
    elif problem == "missing":
        del values[key]
    else:
        values[key] = "unsupported"

    path.write_text(yaml.safe_dump(values))

    with pytest.raises(ValueError, match=reason) as error:
        read_simulation_config(path)

    assert key in str(error.value)
    assert str(path) in str(error.value)


def test_sequencing_depth(small_config):
    config = small_config

    a = simulate(config)
    b = simulate(replace(config, sequencing_depth=2.0))

    assert a.counts.equals(b.counts)
    assert a.transcripts.equals(b.transcripts)

    total_a = a.reference_counts.iloc[:, 1:].to_numpy().sum()
    total_b = b.reference_counts.iloc[:, 1:].to_numpy().sum()

    assert total_b > total_a

def test_placement_attempt_limit_comes_from_config(tmp_path, small_config):
    path = tmp_path / "config.yaml"
    values = asdict(small_config) | {"n_cells": 2, "max_attempts_per_cell": 3}
    path.write_text(yaml.safe_dump(values))
    config = read_simulation_config(path)

    # Always propose the same centre. Only the first candidate can be accepted.
    rng = Mock()
    rng.uniform.return_value = 4.0
    _, cells = place_cells([{"x_um": 4.0, "y_um": 4.0}], config, rng)

    assert len(cells) == 1
    assert rng.uniform.call_count == 2 * 2 * 3  # two coordinate draws per attempt


def test_cli_reproducibility(tmp_path, small_config):
    raw = tmp_path / "raw"
    raw.mkdir()
    write_manifest(
        manifest(raw), tmp_path / "raw_manifest.sha256"
    )

    config = tmp_path / "config.yaml"
    small_config = replace(small_config, reference_n_per_type=5)
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
