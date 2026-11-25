"""Check seed determinism and exact transcript aggregation relationships."""

from dataclasses import replace

import pandas as pd
import pytest

from tissue_sim.aggregation import aggregate_transcripts
from tissue_sim.configuration import SimulationConfig
from tissue_sim.simulation import simulate


@pytest.fixture
def tiny_config():
    """Return a small simulation that does not need external course data."""
    return SimulationConfig(
        n_cells=4,
        field_width_um=16.0,
        field_height_um=16.0,
        cell_radius_um=1.0,
        n_genes=6,
        n_marker_genes=2,
        reference_n_per_type=4,
    )


@pytest.mark.parametrize("seed", [0, 42, 20260101])
def test_same_seed_reproduces_all_simulation_tables(tiny_config, seed):
    """[reproducibility] A repeated seed reproduces every table and parameter."""
    config = replace(tiny_config, seed=seed)
    first = simulate(config)
    second = simulate(config)

    for name in (
        "coordinates",
        "cells",
        "transcripts",
        "counts",
        "composition",
        "domains",
        "reference_counts",
        "reference_labels",
    ):
        pd.testing.assert_frame_equal(
            getattr(first, name), getattr(second, name), check_exact=True
        )
    assert first.params == second.params


def test_seed_controls_random_sampling(tiny_config):
    """[reproducibility] A different seed changes sampled tissue geometry."""
    first = simulate(replace(tiny_config, seed=42))
    second = simulate(replace(tiny_config, seed=43))

    assert not first.cells[["center_x_um", "center_y_um"]].equals(
        second.cells[["center_x_um", "center_y_um"]]
    )


@pytest.fixture
def aggregation_inputs():
    """Use asymmetric counts, repeated contributions, and an empty square."""
    coordinates = pd.DataFrame({"square_id": ["s2", "empty", "s1"]})
    transcripts = pd.DataFrame(
        {
            "square_id": ["s1", "s1", "s1", "s2", "s2", "s2", "s2"],
            "cell_type": ["A", "A", "B", "A", "B", "B", "B"],
            "gene": ["g1", "g1", "g2", "g1", "g2", "g2", "g2"],
        }
    )
    return coordinates, transcripts, ["g2", "unused", "g1"]


def test_aggregation_preserves_exact_gene_and_type_counts(aggregation_inputs):
    """[known-answer] Every transcript contributes once to both output tables."""
    coordinates, transcripts, genes = aggregation_inputs
    counts, composition = aggregate_transcripts(coordinates, transcripts, genes)

    expected_counts = pd.DataFrame(
        {
            "square_id": ["s2", "empty", "s1"],
            "g2": [3, 0, 1],
            "unused": [0, 0, 0],
            "g1": [1, 0, 2],
        }
    )
    pd.testing.assert_frame_equal(counts, expected_counts)

    observed = composition.set_index(["square_id", "cell_type"])[
        "n_transcripts"
    ].to_dict()
    assert observed == {
        ("s1", "A"): 2,
        ("s1", "B"): 1,
        ("s2", "A"): 1,
        ("s2", "B"): 3,
    }
    assert not composition.duplicated(["square_id", "cell_type"]).any()
    assert counts[genes].to_numpy().sum() == len(transcripts)
    assert composition.n_transcripts.sum() == len(transcripts)


def test_duplicating_transcripts_doubles_both_aggregates(aggregation_inputs):
    """[metamorphic] Duplicating input doubles counts and preserves empty squares."""
    coordinates, transcripts, genes = aggregation_inputs
    counts, composition = aggregate_transcripts(coordinates, transcripts, genes)
    doubled_counts, doubled_composition = aggregate_transcripts(
        coordinates,
        pd.concat([transcripts, transcripts], ignore_index=True),
        genes,
    )

    expected_counts = counts.copy()
    expected_counts[genes] *= 2
    pd.testing.assert_frame_equal(doubled_counts, expected_counts)

    expected_composition = composition.copy()
    expected_composition["n_transcripts"] *= 2
    pd.testing.assert_frame_equal(doubled_composition, expected_composition)
def test_four_8um_bins_sum_to_each_16um_bin_gene_by_gene():
    """[property] Four adjacent 8 µm bins exactly form each 16 µm bin."""
    coordinates = pd.DataFrame(
        [
            {
                "square_id": f"s{row}_{col}",
                "array_row": row,
                "array_col": col,
            }
            for row in range(16)
            for col in range(24)
        ]
    )
    genes = ["row_gene", "col_gene"]
    transcripts = pd.DataFrame(
        [
            {
                "square_id": f"s{row}_{col}",
                "cell_type": "A",
                "gene": gene,
            }
            for row in range(16)
            for col in range(24)
            for gene, number in (
                ("row_gene", row + 1),
                ("col_gene", col + 1),
            )
            for _ in range(number)
        ]
    )
    counts, _ = aggregate_transcripts(coordinates, transcripts, genes)
    table = coordinates.merge(
        counts, on="square_id", validate="one_to_one"
    )

    def bin_counts(squares_per_side):
        return table.groupby(
            [
                table.array_row // squares_per_side,
                table.array_col // squares_per_side,
            ]
        )[genes].sum()

    # Each measurement square is 2 µm wide.
    fine = bin_counts(4)
    coarse = bin_counts(8)

    assert fine.index.tolist() == [
        (row, col) for row in range(4) for col in range(6)
    ]
    assert coarse.index.tolist() == [
        (row, col) for row in range(2) for col in range(3)
    ]

    # For bin row i: four rows contain 4i+1 through 4i+4
    # transcripts each, repeated across four columns.
    for (row, col), values in fine.iterrows():
        assert values.tolist() == [64 * row + 40, 64 * col + 40]

    for (row, col), values in coarse.iterrows():
        children = fine.loc[
            [
                (2 * row, 2 * col),
                (2 * row, 2 * col + 1),
                (2 * row + 1, 2 * col),
                (2 * row + 1, 2 * col + 1),
            ],
            genes,
        ]
        pd.testing.assert_series_equal(
            children.sum(), values, check_names=False
        )
        assert values.tolist() == [
            512 * row + 288,
            512 * col + 288,
        ]
