import numpy as np
import pytest
from prometheus_client.utils import INF

from tissue_sim.expression_programs import *


def test_baseline_probabilities_sum_to_one():
    rng = np.random.default_rng(42)

    probabilities = make_baseline_distribution(
        n_genes=100,
        baseline_spread=1.0,
        random_generator=rng,
    )

    assert len(probabilities) == 100
    assert np.all(probabilities >= 0)
    assert np.all(np.isfinite(probabilities))
    assert np.isclose(probabilities.sum(), 1.0)


def test_baseline_is_reproducible():
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)

    p1 = make_baseline_distribution(100, 1.0, rng1)
    p2 = make_baseline_distribution(100, 1.0, rng2)

    np.testing.assert_array_equal(p1, p2)


@pytest.mark.parametrize("n_genes", [0, -2, 1.2, float("inf")])
def test_invalid_number_of_genes_is_rejected(n_genes):
    rng = np.random.default_rng(42)

    with pytest.raises(ValueError):
        make_baseline_distribution(n_genes, 1.0, rng)


def test_negative_spread_is_rejected():
    rng = np.random.default_rng(42)

    with pytest.raises(ValueError):
        make_baseline_distribution(100, -1.0, rng)


def test_marker_genes_are_boosted():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])

    program = boost_marker_genes(
        baseline=baseline,
        marker_genes_indices=[0, 2],
        marker_boost=5.0,
    )

    assert program[0] > baseline[0]
    assert program[2] > baseline[2]


def test_non_marker_genes_keep_their_relative_proportions():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])

    program = boost_marker_genes(baseline, [0, 2], 5.0)

    # The ratio between non-marker genes should remain unchanged.
    assert np.isclose(program[1] / program[3], baseline[1] / baseline[3])
    assert program[1] < baseline[1]
    assert program[3] < baseline[3]


def test_expression_program_sums_to_one():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])

    program = boost_marker_genes(baseline, [0, 2], 5.0)

    assert np.isclose(program.sum(), 1.0)


# boost_marker_genes() does not accidentally change the original baseline distribution
def test_original_baseline_is_not_modified():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])
    original = baseline.copy()

    boost_marker_genes(baseline, [0, 2], 5.0)

    np.testing.assert_array_equal(baseline, original)


def test_invalid_marker_index_is_rejected():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])

    with pytest.raises(ValueError):
        boost_marker_genes(baseline, [4], 5.0)


def test_marker_one_stays_the_same():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])
    program = boost_marker_genes(baseline, [0, 2], 1.0)
    np.testing.assert_array_equal(program, baseline)


def test_marker_boost_below_one_is_rejected():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])

    with pytest.raises(ValueError):
        boost_marker_genes(baseline, [0], 0.5)


def test_marker_genes_are_generated_for_each_cell_type():
    rng = np.random.default_rng(42)

    markers = generate_marker_genes(
        cell_types=["type_a", "type_b"],
        n_genes=10,
        n_marker_genes=3,
        random_generator=rng,
    )

    assert set(markers) == {"type_a", "type_b"}
    print("Marker genes:")
    for cell_type, marker_indices in markers.items():
        print(f"{cell_type}: {sorted(marker_indices)}")

    for marker_indices in markers.values():
        assert len(marker_indices) == 3
        assert len(np.unique(marker_indices)) == 3
        assert np.all(marker_indices >= 0)
        assert np.all(marker_indices < 10)
    assert len(np.intersect1d(markers["type_a"], markers["type_b"])) == 0


def test_marker_gene_selection_is_reproducible():
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)

    markers1 = generate_marker_genes(
        ["type_a", "type_b"], 10, 3, rng1
    )
    markers2 = generate_marker_genes(
        ["type_a", "type_b"], 10, 3, rng2
    )

    for cell_type in markers1:
        np.testing.assert_array_equal(
            markers1[cell_type],
            markers2[cell_type],
        )


@pytest.mark.parametrize("n_marker_genes", [-1, 11, 1.5])
def test_invalid_number_of_marker_genes_is_rejected(n_marker_genes):
    rng = np.random.default_rng(42)

    with pytest.raises(ValueError):
        generate_marker_genes(
            ["type_a"],
            n_genes=10,
            n_marker_genes=n_marker_genes,
            random_generator=rng,
        )

def test_not_enough_genes_for_non_overlapping_markers():
    rng = np.random.default_rng(42)

    with pytest.raises(ValueError, match="Not enough genes"):
        generate_marker_genes(
            cell_types=["type_a", "type_b"],
            n_genes=10,
            n_marker_genes=6,
            random_generator=rng,
        )


def test_expression_programs_are_generated_for_each_cell_type():
    cell_types = ["type_a", "type_b"]
    baseline = np.array([0.1, 0.2, 0.3, 0.4])
    marker_genes = {
        "type_a": [0],
        "type_b": [3],
    }

    programs = generate_cell_type_programs(
        cell_types=cell_types,
        baseline=baseline,
        marker_genes=marker_genes,
        marker_boost=5.0,
    )

    assert set(programs.keys()) == set(cell_types)

    for program in programs.values():
        assert len(program) == len(baseline)
        assert np.all(np.isfinite(program))
        assert np.all(program >= 0)
        assert np.isclose(program.sum(), 1.0)


def test_marker_genes_have_higher_probabilities():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])
    marker_genes = {
        "type_a": [0],
        "type_b": [3],
    }

    programs = generate_cell_type_programs(
        cell_types=["type_a", "type_b"],
        baseline=baseline,
        marker_genes=marker_genes,
        marker_boost=5.0,
    )

    assert programs["type_a"][0] > baseline[0]
    assert programs["type_b"][3] > baseline[3]


def test_generating_programs_does_not_modify_baseline():
    baseline = np.array([0.1, 0.2, 0.3, 0.4])
    original = baseline.copy()

    generate_cell_type_programs(
        cell_types=["type_a", "type_b"],
        baseline=baseline,
        marker_genes={
            "type_a": [0],
            "type_b": [3],
        },
        marker_boost=5.0,
    )

    np.testing.assert_array_equal(baseline, original)


def test_missing_cell_type_markers_are_rejected():
    with pytest.raises(ValueError, match="exactly one entry"):
        generate_cell_type_programs(
            cell_types=["type_a", "type_b"],
            baseline=np.array([0.1, 0.2, 0.3, 0.4]),
            marker_genes={"type_a": [0]},
            marker_boost=5.0,
        )


def test_duplicate_cell_types_are_rejected():
    with pytest.raises(ValueError, match="unique names"):
        generate_cell_type_programs(
            cell_types=["type_a", "type_a"],
            baseline=np.array([0.1, 0.2, 0.3, 0.4]),
            marker_genes={"type_a": [0]},
            marker_boost=5.0,
        )