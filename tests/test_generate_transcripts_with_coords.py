"""Tests for tissue_sim.generate_transcripts_with_coords."""

import json
import math

import numpy as np
import pandas as pd
import pytest

from tissue_sim.generate_transcripts_with_coords import (
    generate_transcripts,
    get_square_id,
    load_expression_programs,
    place_transcripts_in_disk,
)


# =============================================================================
# Fixtures — small reusable inputs
# =============================================================================

@pytest.fixture
def two_genes_two_types():
    """Two cell types, two genes, simple 50/50 programs."""
    programs = {
        "type_a": np.array([0.8, 0.2]),
        "type_b": np.array([0.2, 0.8]),
    }
    gene_names = ["GENE_0001", "GENE_0002"]
    means = {"type_a": 20.0, "type_b": 30.0}
    return programs, gene_names, means


@pytest.fixture
def one_cell_dataframe():
    """A single cell of type_a placed at the center of a 100x100 µm field."""
    return pd.DataFrame([{
        "cell_id": 1,
        "cell_type": "type_a",
        "center_x_um": 50.0,
        "center_y_um": 50.0,
        "radius_um": 5.0,
    }])


@pytest.fixture
def two_cell_dataframe():
    """Two cells of different types."""
    return pd.DataFrame([
        {"cell_id": 1, "cell_type": "type_a", "center_x_um": 25.0, "center_y_um": 50.0, "radius_um": 5.0},
        {"cell_id": 2, "cell_type": "type_b", "center_x_um": 75.0, "center_y_um": 50.0, "radius_um": 5.0},
    ])


# =============================================================================
# Tests for get_square_id()
# =============================================================================

@pytest.mark.parametrize("x,y,size,expected", [
    (5.3,  1.0, 2.0, "sq_r000_c002"),   # x=5.3 -> col 2, y=1.0 -> row 0
    (0.0,  0.0, 2.0, "sq_r000_c000"),   # origin
    (1.99, 1.99, 2.0, "sq_r000_c000"),  # just inside first square
    (2.0,  2.0, 2.0, "sq_r001_c001"),   # exactly on boundary -> next square
    (10.0, 6.0, 2.0, "sq_r003_c005"),   # middle of grid
])
def test_get_square_id_known_values(x, y, size, expected):
    assert get_square_id(x, y, size) == expected


def test_get_square_id_floor_behaviour():
    """Floor division: 3.99 and 3.0 land in the same square."""
    assert get_square_id(3.0, 3.0, 2.0) == get_square_id(3.99, 3.99, 2.0)


# =============================================================================
# Tests for place_transcripts_in_disk()
# =============================================================================

def test_all_points_inside_disk():
    """Every generated point must lie within the disk."""
    rng = np.random.default_rng(0)
    cx, cy, radius = 10.0, 20.0, 5.0
    x, y = place_transcripts_in_disk(1000, cx, cy, radius, rng)
    dist = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    assert np.all(dist <= radius + 1e-10)


def test_points_centered_on_cell():
    """Mean of many points should be close to the cell center."""
    rng = np.random.default_rng(42)
    cx, cy, radius = 50.0, 30.0, 10.0
    x, y = place_transcripts_in_disk(10000, cx, cy, radius, rng)
    assert abs(x.mean() - cx) < 0.5
    assert abs(y.mean() - cy) < 0.5


def test_zero_transcripts_returns_empty():
    """Asking for 0 transcripts returns empty arrays."""
    rng = np.random.default_rng(0)
    x, y = place_transcripts_in_disk(0, 10.0, 10.0, 5.0, rng)
    assert len(x) == 0
    assert len(y) == 0


def test_uniform_radial_distribution():
    """
    For a uniform disk, the fraction of points within r < R/2
    should be close to (R/2)^2 / R^2 = 0.25 (area ratio).
    """
    rng = np.random.default_rng(7)
    cx, cy, radius = 0.0, 0.0, 10.0
    x, y = place_transcripts_in_disk(10000, cx, cy, radius, rng)
    dist = np.sqrt(x ** 2 + y ** 2)
    fraction_inner = (dist < radius / 2).mean()
    assert abs(fraction_inner - 0.25) < 0.02


def test_reproducibility():
    """Same seed gives same coordinates."""
    cx, cy, radius = 5.0, 5.0, 3.0
    x1, y1 = place_transcripts_in_disk(100, cx, cy, radius, np.random.default_rng(99))
    x2, y2 = place_transcripts_in_disk(100, cx, cy, radius, np.random.default_rng(99))
    np.testing.assert_array_equal(x1, x2)
    np.testing.assert_array_equal(y1, y2)


# =============================================================================
# Tests for load_expression_programs()
# =============================================================================

def test_load_expression_programs_basic(tmp_path):
    """Load a valid JSON and check programs, gene names, means."""
    data = {
        "type_a": [0.8, 0.2],
        "type_b": [0.3, 0.7],
        "gene_names": ["GENE_0001", "GENE_0002"],
        "mean_transcripts": {"type_a": 15, "type_b": 25},
    }
    path = tmp_path / "programs.json"
    path.write_text(json.dumps(data), encoding="utf-8")

    programs, gene_names, means = load_expression_programs(path)

    assert set(programs.keys()) == {"type_a", "type_b"}
    assert gene_names == ["GENE_0001", "GENE_0002"]
    assert means == {"type_a": 15, "type_b": 25}
    np.testing.assert_allclose(programs["type_a"], [0.8, 0.2])
    np.testing.assert_allclose(programs["type_b"], [0.3, 0.7])


def test_load_expression_programs_default_gene_names(tmp_path):
    """If gene_names is absent, defaults to gene_0, gene_1, ..."""
    data = {"type_a": [0.5, 0.3, 0.2]}
    path = tmp_path / "programs.json"
    path.write_text(json.dumps(data), encoding="utf-8")

    _, gene_names, _ = load_expression_programs(path)
    assert gene_names == ["gene_0", "gene_1", "gene_2"]


def test_load_expression_programs_normalises(tmp_path):
    """Programs that don't sum to exactly 1 are normalised."""
    data = {"type_a": [2.0, 2.0]}   # sums to 4, should normalise to [0.5, 0.5]
    path = tmp_path / "programs.json"
    path.write_text(json.dumps(data), encoding="utf-8")

    programs, _, _ = load_expression_programs(path)
    np.testing.assert_allclose(programs["type_a"], [0.5, 0.5])


def test_load_expression_programs_empty_means(tmp_path):
    """If mean_transcripts key is absent, returns empty dict."""
    data = {"type_a": [1.0]}
    path = tmp_path / "programs.json"
    path.write_text(json.dumps(data), encoding="utf-8")

    _, _, means = load_expression_programs(path)
    assert means == {}


# =============================================================================
# Tests for generate_transcripts()
# =============================================================================

def test_output_schema(one_cell_dataframe, two_genes_two_types):
    """Output DataFrame has the expected columns."""
    programs, gene_names, means = two_genes_two_types
    rng = np.random.default_rng(0)
    result = generate_transcripts(
        one_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    expected_columns = ["transcript_id", "square_id", "cell_id", "cell_type", "gene", "x_um", "y_um"]
    assert result.columns.tolist() == expected_columns


def test_output_is_empty_dataframe_for_empty_cells(two_genes_two_types):
    """Empty cell table returns empty DataFrame with correct columns."""
    programs, gene_names, means = two_genes_two_types
    empty_cells = pd.DataFrame(columns=["cell_id", "cell_type", "center_x_um", "center_y_um", "radius_um"])
    rng = np.random.default_rng(0)
    result = generate_transcripts(
        empty_cells, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    assert result.empty
    assert "transcript_id" in result.columns


def test_transcript_ids_are_unique(one_cell_dataframe, two_genes_two_types):
    """transcript_id must be unique across all transcripts."""
    programs, gene_names, means = two_genes_two_types
    rng = np.random.default_rng(1)
    result = generate_transcripts(
        one_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    assert result["transcript_id"].is_unique


def test_all_transcripts_inside_field(two_cell_dataframe, two_genes_two_types):
    """All transcript coordinates must be within the field boundaries."""
    programs, gene_names, means = two_genes_two_types
    rng = np.random.default_rng(2)
    result = generate_transcripts(
        two_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    assert (result["x_um"] >= 0).all()
    assert (result["x_um"] < 100.0).all()
    assert (result["y_um"] >= 0).all()
    assert (result["y_um"] < 100.0).all()


def test_all_transcripts_inside_their_cell_disk(two_cell_dataframe, two_genes_two_types):
    """Each transcript must lie within its cell's disk."""
    programs, gene_names, means = two_genes_two_types
    rng = np.random.default_rng(3)
    result = generate_transcripts(
        two_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    # merge with cell geometry to check distances
    cell_geom = two_cell_dataframe[["cell_id", "center_x_um", "center_y_um", "radius_um"]]
    merged = result.merge(cell_geom, on="cell_id")
    dist = np.sqrt(
        (merged["x_um"] - merged["center_x_um"]) ** 2 +
        (merged["y_um"] - merged["center_y_um"]) ** 2
    )
    assert (dist <= merged["radius_um"] + 1e-10).all()


def test_cell_type_matches_cell(two_cell_dataframe, two_genes_two_types):
    """Each transcript's cell_type must match its cell_id."""
    programs, gene_names, means = two_genes_two_types
    rng = np.random.default_rng(4)
    result = generate_transcripts(
        two_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    cell_types = two_cell_dataframe.set_index("cell_id")["cell_type"]
    for _, row in result.iterrows():
        assert row["cell_type"] == cell_types[row["cell_id"]]


def test_genes_come_from_program(two_cell_dataframe, two_genes_two_types):
    """Every gene in the output must be in gene_names."""
    programs, gene_names, means = two_genes_two_types
    rng = np.random.default_rng(5)
    result = generate_transcripts(
        two_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    assert result["gene"].isin(gene_names).all()


def test_square_id_matches_coordinates(one_cell_dataframe, two_genes_two_types):
    """square_id must correspond to the transcript's actual x/y position."""
    programs, gene_names, means = two_genes_two_types
    rng = np.random.default_rng(6)
    result = generate_transcripts(
        one_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0, rng
    )
    for _, row in result.iterrows():
        expected = get_square_id(row["x_um"], row["y_um"], 2.0)
        assert row["square_id"] == expected


def test_reproducibility_with_seed(one_cell_dataframe, two_genes_two_types):
    """Same seed produces identical transcript tables."""
    programs, gene_names, means = two_genes_two_types
    result1 = generate_transcripts(
        one_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0,
        np.random.default_rng(42)
    )
    result2 = generate_transcripts(
        one_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0,
        np.random.default_rng(42)
    )
    pd.testing.assert_frame_equal(result1, result2)


def test_different_seeds_give_different_results(one_cell_dataframe, two_genes_two_types):
    """Different seeds should produce different coordinates."""
    programs, gene_names, means = two_genes_two_types
    result1 = generate_transcripts(
        one_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0,
        np.random.default_rng(1)
    )
    result2 = generate_transcripts(
        one_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0,
        np.random.default_rng(2)
    )
    assert not result1["x_um"].equals(result2["x_um"])


def test_transcript_count_is_poisson_distributed():
    """
    Over many runs, the mean transcript count should be close
    to the Poisson mean set in mean_transcripts.
    """
    programs = {"type_a": np.array([1.0])}
    gene_names = ["GENE_0001"]
    lam = 50.0
    means = {"type_a": lam}
    cells = pd.DataFrame([{
        "cell_id": 1, "cell_type": "type_a",
        "center_x_um": 50.0, "center_y_um": 50.0, "radius_um": 5.0,
    }])
    counts = []
    for seed in range(200):
        result = generate_transcripts(
            cells, programs, gene_names, means, 2.0, 100.0, 100.0,
            np.random.default_rng(seed)
        )
        counts.append(len(result))
    observed_mean = np.mean(counts)
    # mean should be within 3 standard errors of lambda
    assert abs(observed_mean - lam) < 3 * math.sqrt(lam / 200)


def test_cell_at_field_boundary_transcripts_clipped():
    """
    A cell whose disk extends outside the field: only transcripts
    inside the field are kept.
    """
    programs = {"type_a": np.array([1.0])}
    gene_names = ["GENE_0001"]
    means = {"type_a": 100.0}
    # cell center is at x=1, right at the left edge — disk extends outside
    cells = pd.DataFrame([{
        "cell_id": 1, "cell_type": "type_a",
        "center_x_um": 1.0, "center_y_um": 50.0, "radius_um": 5.0,
    }])
    result = generate_transcripts(
        cells, programs, gene_names, means, 2.0, 100.0, 100.0,
        np.random.default_rng(0)
    )
    assert (result["x_um"] >= 0).all()


def test_input_dataframe_is_not_modified(two_cell_dataframe, two_genes_two_types):
    """generate_transcripts must not modify the input cells DataFrame."""
    programs, gene_names, means = two_genes_two_types
    original = two_cell_dataframe.copy(deep=True)
    generate_transcripts(
        two_cell_dataframe, programs, gene_names, means, 2.0, 100.0, 100.0,
        np.random.default_rng(0)
    )
    pd.testing.assert_frame_equal(two_cell_dataframe, original)