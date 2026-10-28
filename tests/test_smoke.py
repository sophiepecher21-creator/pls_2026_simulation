"""Check the connected pipeline, exact aggregation, and reproducibility."""

from dataclasses import replace
import json

import numpy as np
import pandas as pd
import pytest
import yaml

from tissue_sim.aggregation import aggregate_transcripts
from tissue_sim.configuration import SimulationConfig, read_simulation_config
from tissue_sim.generate_transcripts_with_coords import (
    generate_transcripts,
    get_square_id,
)
from tissue_sim.provenance import manifest, write_manifest
from tissue_sim.simulation import run, simulate


def small_config():
    """Return a fast, low-density simulation configuration."""
    return SimulationConfig(
        n_cells=4,
        field_width_um=16.0,
        field_height_um=16.0,
        cell_radius_um=1.0,
        n_genes=6,
        n_marker_genes=2,
    )


def test_pipeline_reproducible_and_conserves_counts():
    a = simulate(small_config())
    b = simulate(small_config())
    for name in (
        "coordinates",
        "cells",
        "transcripts",
        "counts",
        "composition",
        "domains",
    ):
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name))
    assert len(a.coordinates) == 64
    assert a.counts.iloc[:, 1:].to_numpy().sum() == len(a.transcripts)
    totals = (
        a.composition.groupby("square_id")
        .n_transcripts.sum()
        .reindex(a.coordinates.square_id, fill_value=0)
    )
    np.testing.assert_array_equal(totals, a.counts.iloc[:, 1:].sum(axis=1))
    assert not a.composition.duplicated(["square_id", "cell_type"]).any()


def test_known_answer_aggregation_includes_empty_squares():
    coordinates = pd.DataFrame({"square_id": ["s1", "s2"]})
    transcripts = pd.DataFrame(
        {
            "square_id": ["s1", "s1", "s1"],
            "cell_type": ["A", "A", "B"],
            "gene": ["g1", "g1", "g2"],
        }
    )
    counts, comp = aggregate_transcripts(coordinates, transcripts, ["g1", "g2"])
    assert counts.iloc[:, 1:].values.tolist() == [[2, 1], [0, 0]]
    assert comp.n_transcripts.tolist() == [2, 1]


def test_cells_without_square_centres_still_generate_and_boundaries_are_clipped():
    cells = pd.DataFrame(
        {
            "cell_id": [1],
            "cell_type": ["A"],
            "center_x_um": [0.01],
            "center_y_um": [0.01],
            "radius_um": [0.1],
        }
    )
    t = generate_transcripts(
        cells,
        {"A": np.array([1.0])},
        ["g"],
        {"A": 1000.0},
        2.0,
        4.0,
        4.0,
        np.random.default_rng(42),
    )
    assert 0 < len(t) < 1000
    assert t.x_um.ge(0).all() and t.y_um.ge(0).all()
    assert t.square_id.eq("sq_r000_c000").all()
    assert get_square_id(-0.1, 0, 2) == "sq_r000_c-01"


def test_zero_transcripts_retains_full_grid():
    result = simulate(replace(small_config(), mean_transcripts_per_type=(0.0, 0.0)))
    assert result.transcripts.empty and result.composition.empty
    assert result.counts.shape == (64, 7)
    assert result.counts.iloc[:, 1:].to_numpy().sum() == 0


def prepare_run(tmp_path):
    """Create temporary inputs so no supplied data is modified."""
    from dataclasses import asdict

    raw = tmp_path / "data/raw"
    raw.mkdir(parents=True)
    (raw / "example.csv").write_text("value\n7\n")
    write_manifest(manifest(raw), raw.parent / "raw_manifest.sha256")
    config = tmp_path / "config.yaml"
    values = asdict(small_config())
    # safe_dump produces ordinary lists for tuple settings.
    config.write_text(yaml.safe_dump(values))
    return config, raw, tmp_path / "results"


def test_run_exports_and_refuses_overwrite(tmp_path):
    config, raw, out = prepare_run(tmp_path)
    run(config, raw, out)
    for f in [
        "coordinates.csv",
        "counts.csv",
        "ground_truth/composition.csv",
        "ground_truth/domains.csv",
        "ground_truth/params.json",
        "diagnostics/cells.csv",
        "diagnostics/transcripts.csv",
        "run_metadata.json",
    ]:
        assert (out / f).is_file()
    params = json.loads((out / "ground_truth/params.json").read_text())
    assert params["seed"] == small_config().seed
    with pytest.raises(FileExistsError):
        run(config, raw, out)
    other = tmp_path / "repeat"
    run(config, raw, other)
    for path in out.rglob("*"):
        if path.is_file() and path.name != "run_metadata.json":
            assert path.read_bytes() == (other / path.relative_to(out)).read_bytes()


def test_changed_raw_data_prevents_output(tmp_path):
    config, raw, out = prepare_run(tmp_path)
    (raw / "example.csv").write_text("changed")
    with pytest.raises(ValueError, match="Raw-data verification failed"):
        run(config, raw, out)
    assert not out.exists()


def test_config_safe_loading(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text('seed: "42"\n')
    with pytest.raises(ValueError, match="seed"):
        read_simulation_config(path)
    path.write_text("!!python/object/apply:os.system [echo unsafe]")
    with pytest.raises(yaml.YAMLError):
        read_simulation_config(path)


def test_exact_8_to_16_aggregation():
    result = simulate(small_config())
    genes = result.counts.columns[1:]
    table = result.counts.merge(result.coordinates, on="square_id")

    def bins(k):
        return table.groupby([table.array_row // k, table.array_col // k])[
            list(genes)
        ].sum()

    fine = bins(4)
    coarse = bins(8)
    np.testing.assert_array_equal(fine.sum().to_numpy(), coarse.iloc[0].to_numpy())
