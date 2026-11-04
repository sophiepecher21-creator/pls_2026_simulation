"""Test domain labels, cell extraction, and cell-type assignment.

Check returned DataFrames, reproducibility of random sampling, and that input
DataFrames remain unchanged. Detailed placement validation belongs upstream.
"""

import numpy as np
import pandas as pd
import pytest

from tissue_sim.domains import (
    assign_cell_types,
    build_domains_dataframe,
    cells_from_placement,
    get_domain_id,
)


def placement():
    """Two squares share a cell; the third square is empty."""
    return pd.DataFrame({
        "square_id": ["s1", "s2", "s3"],
        "x_um": [98, 102, 150],
        "y_um": [1, 1, 1],
        "cell_id": ["1", "1", "0"],
        "center_x_um": [101.0, 101.0, np.nan],
        "center_y_um": [1.0, 1.0, np.nan],
        "radius_um": [5.0, 5.0, np.nan],
    })


def test_domain_boundaries():
    """Boundaries belong to the right-hand strip; y does not affect labels."""
    labels = [
        get_domain_id(x, 1, 2, 200)
        for x in [0, 99, 100, 199]
    ]

    assert labels == [0, 0, 1, 1]
    assert get_domain_id(100, 50, 2, 200) == 1


def test_single_domain():
    """With one domain, every valid position receives zero."""
    labels = [
        get_domain_id(x, 1, 1, 200)
        for x in [0, 100, 199]
    ]

    assert labels == [0, 0, 0]


def test_square_domains():
    """Square labels preserve row order, index and input data."""
    frame = placement()
    frame.index = [7, 7, 2]
    before = frame.copy(deep=True)

    result = build_domains_dataframe(frame, 2, 200)

    assert result.columns.tolist() == ["square_id", "domain"]
    assert result["square_id"].tolist() == ["s1", "s2", "s3"]
    assert result["domain"].tolist() == ["0", "1", "1"]
    assert result.index.tolist() == [7, 7, 2]
    pd.testing.assert_frame_equal(frame, before)


def test_cells_from_placement():
    """Repeated records collapse and empty squares are excluded."""
    frame = placement()
    before = frame.copy(deep=True)

    cells = cells_from_placement(frame)

    assert cells["cell_id"].tolist() == ["1"]
    assert cells["center_x_um"].tolist() == [101.0]
    assert cells["radius_um"].tolist() == [5.0]
    pd.testing.assert_frame_equal(frame, before)


def test_conflicting_cell_records():
    """One cell ID cannot describe two different centres."""
    frame = placement()
    frame.loc[1, "center_x_um"] = 110

    with pytest.raises(ValueError, match="Conflicting"):
        cells_from_placement(frame)


def test_cell_types():
    """Cell types use cell centres and preserve existing data."""
    cells = pd.DataFrame({
        "cell_id": ["c1", "c2", "c3"],
        "center_x_um": [98, 100, 150],
        "center_y_um": [1, 1, 1],
        "radius_um": [5, 5, 5],
    }, index=[7, 7, 2])
    before = cells.copy(deep=True)

    result = assign_cell_types(
        cells,
        2,
        200,
        ["A", "B"],
        np.eye(2),
        np.random.default_rng(42),
    )

    assert result["domain"].tolist() == [0, 1, 1]
    assert result["cell_type"].tolist() == ["A", "B", "B"]
    pd.testing.assert_frame_equal(result[cells.columns], before)
    pd.testing.assert_frame_equal(cells, before)


def test_sampling_reproducibility():
    """The same seed produces the same sampled types."""
    cells = pd.DataFrame({
        "center_x_um": [1] * 20,
        "center_y_um": [1] * 20,
    })

    first = assign_cell_types(
        cells, 1, 200, ["A", "B"], [[0.5, 0.5]],
        np.random.default_rng(42),
    )
    second = assign_cell_types(
        cells, 1, 200, ["A", "B"], [[0.5, 0.5]],
        np.random.default_rng(42),
    )

    pd.testing.assert_frame_equal(first, second)


def test_empty_tables():
    """Empty inputs return empty tables with the required columns."""
    cells = pd.DataFrame(columns=["center_x_um", "center_y_um"])
    squares = pd.DataFrame(columns=["square_id", "x_um", "y_um"])

    annotated = assign_cell_types(
        cells, 1, 200, ["A"], [[1]],
        np.random.default_rng(42),
    )
    labels = build_domains_dataframe(squares, 1, 200)

    assert annotated.empty
    assert annotated.columns.tolist() == [
        "center_x_um", "center_y_um", "domain", "cell_type"
    ]
    assert labels.empty
    assert labels.columns.tolist() == ["square_id", "domain"]


def test_sampling_frequencies():
    """Observed frequencies roughly follow each domain's probabilities."""
    n = 4000
    cells = pd.DataFrame({
        "center_x_um": np.repeat([50, 150], n),
        "center_y_um": 1,
    })

    result = assign_cell_types(
        cells,
        2,
        200,
        ["A", "B"],
        [[0.8, 0.2], [0.2, 0.8]],
        np.random.default_rng(42),
    )

    for domain, expected in [(0, 0.8), (1, 0.2)]:
        types = result.loc[result["domain"] == domain, "cell_type"]
        observed = (types == "A").mean()

        assert abs(observed - expected) < 0.04


def test_shared_rng_advances():
    """Sampling advances the supplied generator."""
    cells = pd.DataFrame({
        "center_x_um": [1] * 20,
        "center_y_um": [1] * 20,
    })
    rng = np.random.default_rng(42)
    reference = np.random.default_rng(42)

    assign_cell_types(
        cells, 1, 200, ["A", "B"], [[0.5, 0.5]], rng,
    )

    assert rng.random() != reference.random()


def test_dataframe_pipeline():
    """Placement, extraction and annotation work together in memory."""
    frame = placement()
    before = frame.copy(deep=True)

    cells = cells_from_placement(frame)
    annotated_cells = assign_cell_types(
        cells,
        2,
        200,
        ["A", "B"],
        np.eye(2),
        np.random.default_rng(42),
    )
    square_domains = build_domains_dataframe(frame, 2, 200)

    assert isinstance(annotated_cells, pd.DataFrame)
    assert isinstance(square_domains, pd.DataFrame)
    assert annotated_cells["cell_id"].tolist() == ["1"]
    assert annotated_cells["domain"].tolist() == [1]
    assert annotated_cells["cell_type"].tolist() == ["B"]
    assert square_domains["domain"].tolist() == ["0", "1", "1"]
    pd.testing.assert_frame_equal(frame, before)