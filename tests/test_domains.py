"""Codex-restored tests for the spatial-domain component only."""

import copy

import numpy as np
import pandas as pd
import pytest

from tissue_sim.domains import (
    assign_cell_types,
    build_domains_dataframe,
    domains_from_csv,
    get_domain_id,
)


def grid(rows, cols):
    """Independent mock 2 micrometer measurement grid."""
    row, col = np.indices((rows, cols))
    return pd.DataFrame(
        {
            "square_id": [f"sq_r{r:03d}_c{c:03d}" for r, c in zip(row.flat, col.flat, strict=True)],
            "array_row": row.ravel(),
            "array_col": col.ravel(),
            "x_um": 2 * col.ravel() + 1,
            "y_um": 2 * row.ravel() + 1,
        }
    )


def annotate(frame, probabilities=None, rng=None):
    """Short fixture helper using two domains and two cell types."""
    if probabilities is None:
        probabilities = np.eye(2)
    if rng is None:
        rng = np.random.default_rng(42)
    return assign_cell_types(frame, 2, 200, ["A", "B"], probabilities, rng)


@pytest.mark.parametrize(
    "x,expected", [(0, 0), (50, 0), (99, 0), (100, 1), (150, 1), (200, 2), (299, 2)]
)
def test_known_boundaries(x, expected):
    assert get_domain_id(x, 50, 3, 300) == expected


@pytest.mark.parametrize("width,K", [(1.0, 10), (320.0, 3), (37.0, 7)])
def test_floating_point_boundary_neighbors(width, K):
    for i in range(1, K):
        boundary = width * i / K
        assert get_domain_id(np.nextafter(boundary, 0), 0, K, width) == i - 1
        assert get_domain_id(boundary, 0, K, width) == i
        assert get_domain_id(np.nextafter(boundary, width), 0, K, width) == i
    assert get_domain_id(np.nextafter(width, 0), 0, K, width) == K - 1


@pytest.mark.parametrize("K", [0, -1, 1.5, True, np.bool_(True), "3", None])
def test_invalid_K(K):
    with pytest.raises((TypeError, ValueError), match="K"):
        get_domain_id(0, 0, K, 10)


@pytest.mark.parametrize("width", [0, -10, np.nan, np.inf, -np.inf, "10", None, True])
def test_invalid_width(width):
    with pytest.raises((TypeError, ValueError), match="field_width_um"):
        get_domain_id(0, 0, 1, width)


@pytest.mark.parametrize(
    "x,y",
    [
        (-1, 0),
        (10, 0),
        (11, 0),
        (np.nan, 0),
        (np.inf, 0),
        (0, np.nan),
        (0, np.inf),
        (0, -1),
        (None, 0),
        ("2", 0),
        (True, 0),
    ],
)
def test_invalid_coordinates(x, y):
    with pytest.raises((TypeError, ValueError)):
        get_domain_id(x, y, 2, 10)


def test_single_domain_and_y_independence():
    for x in [0, 4, np.nextafter(10.0, 0)]:
        assert get_domain_id(x, 0, np.int64(1), 10) == 0
        assert get_domain_id(x, 10000, 1, 10) == 0
    assert get_domain_id(7, 0, 2, 10) == get_domain_id(7, 10000, 2, 10) == 1


def test_cell_crossing_boundary_and_known_types():
    cells = pd.DataFrame(
        {
            "cell_id": ["left", "boundary", "right"],
            "center_x_um": [98, 100, 102],
            "center_y_um": [10] * 3,
            "radius_um": [5] * 3,
        }
    )
    result = annotate(cells)
    assert result["domain"].tolist() == [0, 1, 1]
    assert result["cell_type"].tolist() == ["A", "B", "B"]
    assert get_domain_id(103, 10, 2, 200) == 1


@pytest.mark.parametrize(
    "rows,cols,K", [(160, 160, 3), (20, 20, 3), (3, 5, 1), (4, 3, 10)]
)
def test_complete_coverage_and_spatial_contiguity(rows, cols, K):
    coordinates = grid(rows, cols)
    before = coordinates.copy(deep=True)
    labels = build_domains_dataframe(coordinates, K, 2 * cols)
    pd.testing.assert_frame_equal(coordinates, before)
    assert labels.columns.tolist() == ["square_id", "domain"]
    assert labels["square_id"].tolist() == coordinates["square_id"].tolist()
    assert len(labels) == rows * cols
    assert labels["square_id"].is_unique
    assert not labels.isna().any().any()
    assert labels["domain"].map(lambda d: isinstance(d, str)).all()
    spatial = labels["domain"].astype(int).to_numpy().reshape(rows, cols)
    assert ((spatial >= 0) & (spatial < K)).all()
    # Check spatial arrangement: all rows agree and each occupied strip has
    # consecutive columns, giving one connected rectangle in the 2-D grid.
    assert np.all(spatial == spatial[0])
    assert np.all(np.diff(spatial, axis=1) >= 0)
    for domain in range(K):
        columns = np.flatnonzero(spatial[0] == domain)
        assert np.all(np.diff(columns) == 1)
    if (rows, cols, K) == (160, 160, 3):
        assert labels["domain"].value_counts().sort_index().tolist() == [
            8480,
            8640,
            8480,
        ]
    if (rows, cols, K) == (20, 20, 3):
        assert labels["domain"].value_counts().sort_index().tolist() == [140, 120, 140]
    if K > cols:
        assert len(np.unique(spatial)) < K  # Domains without square centers are valid.


@pytest.mark.parametrize(
    "probabilities",
    [
        [[1, 0]],
        [0.5, 0.5],
        [[[0.5, 0.5]], [[0.5, 0.5]]],
        [[0.5, 0.5], [0.5]],
        [[-0.1, 1.1], [0.5, 0.5]],
        [[np.nan, 0], [0.5, 0.5]],
        [[np.inf, 0], [0.5, 0.5]],
        [[0.4, 0.4], [0.5, 0.5]],
        [[0, 0], [0.5, 0.5]],
        [[0.5, 0.5], [0.5, 0.500001]],
        [["bad", 0], [0.5, 0.5]],
    ],
)
def test_invalid_probabilities_without_rng_consumption(probabilities):
    cells = pd.DataFrame({"center_x_um": [1], "center_y_um": [1]})
    rng = np.random.default_rng(42)
    state = copy.deepcopy(rng.bit_generator.state)
    with pytest.raises(ValueError, match="domain_probabilities"):
        annotate(cells, probabilities, rng)
    assert rng.bit_generator.state == state


def test_probability_roundoff_and_immutability():
    probabilities = np.array([[0.5, 0.500000001], [0.5, 0.5]])
    original = probabilities.copy()
    result = annotate(pd.DataFrame({"center_x_um": [1], "center_y_um": [1]}), probabilities)
    assert result["cell_type"].iloc[0] in ["A", "B"]
    np.testing.assert_array_equal(probabilities, original)


@pytest.mark.parametrize("names", [[], ["A", "A"], [""], [" "], [1], "AB", None])
def test_invalid_cell_type_names(names):
    with pytest.raises((TypeError, ValueError), match="cell_types"):
        assign_cell_types(
            pd.DataFrame({"center_x_um": [1], "center_y_um": [1]}),
            1,
            10,
            names,
            [[1]],
            np.random.default_rng(1),
        )


def test_cell_square_consistency_immutability_and_order():
    cells = pd.DataFrame(
        {
            "cell_id": ["c2", "c1", "c3"],
            "center_x_um": [150, 98, 100],
            "center_y_um": [1, 1, 1],
            "radius_um": [4, 5, 3],
            "extra": ["keep", "also keep", "unchanged"],
        },
        index=[7, 7, 2],
    )
    before = cells.copy(deep=True)
    result = annotate(cells)
    pd.testing.assert_frame_equal(cells, before)
    pd.testing.assert_frame_equal(result[cells.columns], cells)
    squares = cells.rename(columns={
        "cell_id": "square_id",
        "center_x_um": "x_um",
        "center_y_um": "y_um",
    })    
    square_before = squares.copy(deep=True)
    labels = build_domains_dataframe(squares, 2, 200)
    pd.testing.assert_frame_equal(squares, square_before)
    assert labels.index.tolist() == [7, 7, 2]
    assert labels["square_id"].tolist() == ["c2", "c1", "c3"]
    assert (
        labels["domain"].astype(int).tolist() == result["domain"].tolist() == [1, 0, 1]
    )


@pytest.mark.parametrize(
    "ids",
    [["sq", "sq"], ["sq", None], ["sq", np.nan], ["sq", ""], ["sq", " "], ["sq", 1]],
)
def test_invalid_square_ids(ids):
    with pytest.raises(ValueError, match="square_id"):
        build_domains_dataframe(
            pd.DataFrame({"square_id": ids, "x_um": [1, 3], "y_um": [1, 1]}), 2, 10
        )


@pytest.mark.parametrize("bad", [np.nan, np.inf, None, "bad", -1, 200])
def test_invalid_table_coordinates_without_rng_consumption(bad):
    squares = pd.DataFrame({
        "square_id": ["s1", "s2"],
        "x_um": [1, bad],
        "y_um": [1, 1],
    })
    cells = squares.rename(columns={
        "x_um": "center_x_um",
        "y_um": "center_y_um",
    })

    rng = np.random.default_rng(123)
    state = copy.deepcopy(rng.bit_generator.state)

    with pytest.raises(ValueError, match="row position 1"):
        annotate(cells, rng=rng)
    assert rng.bit_generator.state == state

    with pytest.raises(ValueError, match="row position 1"):
        build_domains_dataframe(squares, 2, 200)


def test_empty_tables_and_geometry_validation():
    cells = pd.DataFrame(columns=["center_x_um", "center_y_um", "cell_id"])
    rng = np.random.default_rng(1)
    state = copy.deepcopy(rng.bit_generator.state)
    result = assign_cell_types(cells, 1, 10, ["A"], [[1]], rng)
    assert result.empty
    assert pd.api.types.is_integer_dtype(result["domain"])
    assert rng.bit_generator.state == state
    squares = cells.rename(columns={
        "cell_id": "square_id",
        "center_x_um": "x_um",
        "center_y_um": "y_um",
    })
    labels = build_domains_dataframe(squares, 1, 10)
    assert labels.empty
    assert labels.columns.tolist() == ["square_id", "domain"]
    for K, width in [(0, 10), (1, 0), (True, 10), (1, np.nan)]:
        with pytest.raises((TypeError, ValueError)):
            assign_cell_types(cells, K, width, ["A"], [[1]], rng)
        with pytest.raises((TypeError, ValueError)):
            build_domains_dataframe(squares, K, width)


def test_required_columns_annotation_collisions_and_rng_type():
    cells = pd.DataFrame({
        "center_x_um": [1],
        "center_y_um": [1],
    })
    squares = pd.DataFrame({"x_um": [1], "y_um": [1]})

    for column in ["center_x_um", "center_y_um"]:
        with pytest.raises(
            ValueError, match="center_x_um and center_y_um"
        ):
            annotate(cells.drop(columns=column))

    for column in ["x_um", "y_um"]:
        with pytest.raises(ValueError):
            build_domains_dataframe(squares.drop(columns=column), 1, 10)

    with pytest.raises(ValueError, match="square_id"):
        build_domains_dataframe(squares, 1, 10)

    for column in ["domain", "cell_type"]:
        with pytest.raises(ValueError, match="already contains"):
            annotate(cells.assign(**{column: "old"}))

    with pytest.raises(TypeError, match="rng"):
        assign_cell_types(cells, 1, 10, ["A"], [[1]], None)

    duplicated = pd.concat([cells, cells[["center_x_um"]]], axis=1)
    with pytest.raises(ValueError, match="column names"):
        annotate(duplicated)


@pytest.mark.parametrize("width,K,n_types", [(13, 1, 1), (80, 4, 2), (15, 7, 5)])
def test_general_geometry_and_type_counts(width, K, n_types):
    names = [f"type_{i}" for i in range(n_types)]
    cells = pd.DataFrame(
        {"center_x_um": width * (np.arange(K) + 0.5) / K, "center_y_um": np.zeros(K)}
    )
    probabilities = np.zeros((K, n_types))
    probabilities[np.arange(K), np.arange(K) % n_types] = 1
    result = assign_cell_types(
        cells, K, width, names, probabilities, np.random.default_rng(5)
    )
    assert result["domain"].tolist() == list(range(K))
    assert result["cell_type"].tolist() == [names[i % n_types] for i in range(K)]


def test_statistical_frequencies_for_every_domain():
    n = 6000
    probabilities = np.array([[0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]])
    cells = pd.DataFrame({"center_x_um": np.repeat([50, 150, 250], n), "center_y_um": 1})
    result = assign_cell_types(
        cells, 3, 300, ["A", "B", "C"], probabilities, np.random.default_rng(42)
    )
    for domain, expected in enumerate(probabilities):
        observed = (
            result.loc[result["domain"] == domain, "cell_type"]
            .value_counts()
            .reindex(["A", "B", "C"], fill_value=0)
            .to_numpy()
            / n
        )
        # Six binomial standard errors: conservative across all nine comparisons.
        tolerance = 6 * np.sqrt(expected * (1 - expected) / n)
        assert np.all(np.abs(observed - expected) < tolerance)


def test_reproducibility_and_shared_rng_advancement():
    cells = pd.DataFrame({"center_x_um": np.ones(100), "center_y_um": np.ones(100)})
    rng = np.random.default_rng(42)
    probabilities = [[0.5, 0.5], [0.5, 0.5]]
    first = annotate(cells, probabilities, rng)
    repeat = annotate(cells, probabilities, np.random.default_rng(42))
    pd.testing.assert_frame_equal(first, repeat)
    second = annotate(cells, probabilities, rng)
    assert not first["cell_type"].equals(second["cell_type"])


def test_csv_contract_round_trip(tmp_path):
    labels = build_domains_dataframe(grid(3, 5), 3, 10)
    path = tmp_path / "domains.csv"
    labels.to_csv(path, index=False)
    assert path.read_text().splitlines()[0] == "square_id,domain"
    loaded = pd.read_csv(path, dtype={"square_id": "string", "domain": "string"})
    assert loaded["square_id"].tolist() == labels["square_id"].tolist()
    assert loaded["domain"].tolist() == labels["domain"].tolist()


@pytest.mark.parametrize("names", [{"A", "B"}, frozenset({"A", "B"}), {"A": 1, "B": 2}])
def test_reject_unordered_or_mapping_cell_types(names):
    cells = pd.DataFrame({"center_x_um": [1], "center_y_um": [1]})
    rng = np.random.default_rng(42)
    state = copy.deepcopy(rng.bit_generator.state)
    with pytest.raises(TypeError, match="cell_types"):
        assign_cell_types(cells, 1, 10, names, [[1, 0]], rng)
    assert rng.bit_generator.state == state


@pytest.mark.parametrize(
    "probabilities",
    [
        np.array([[1 + 2j, 0j]]),
        [[1 + 2j, 0j]],
        np.array([[1 + 2j, 0j]], dtype=object),
    ],
)
def test_reject_complex_probabilities(probabilities):
    cells = pd.DataFrame({"center_x_um": [1], "center_y_um": [1]})
    rng = np.random.default_rng(42)
    state = copy.deepcopy(rng.bit_generator.state)
    with pytest.raises(ValueError, match="domain_probabilities"):
        assign_cell_types(cells, 1, 10, ["A", "B"], probabilities, rng)
    assert rng.bit_generator.state == state


@pytest.mark.parametrize("names", [["A", "B"], ("A", "B"), np.array(["A", "B"])])
def test_ordered_cell_types_keep_probability_column_mapping(names):
    cells = pd.DataFrame({"center_x_um": [50, 150], "center_y_um": [1, 1]})
    result = assign_cell_types(
        cells, 2, 200, names, np.eye(2), np.random.default_rng(42)
    )
    assert result["cell_type"].tolist() == ["A", "B"]


def test_coordinates_csv_is_updated_in_place(tmp_path):
    path = tmp_path / "coordinates.csv"
    original = pd.DataFrame({
        "square_id": ["s1", "s2", "s3"],
        "array_row": [0, 0, 0],
        "array_col": [0, 1, 2],
        "x_um": [98, 102, 150],
        "y_um": [1, 1, 1],
        "cell_id": ["1", "1", "0"],
        "center_x_um": [101.0, 101.0, np.nan],
        "center_y_um": [1.0, 1.0, np.nan],
        "radius_um": [5.0, 5.0, np.nan],
    })
    original.to_csv(path, index=False)

    domains_from_csv(
        path, 2, 200, ["A", "B"], np.eye(2),
        np.random.default_rng(42),
    )

    updated = pd.read_csv(
        path,
        dtype={
            "square_id": "string",
            "cell_id": "string",
            "domain": "string",
            "cell_type": "string",
        },
    )

    # Original data and row order are preserved.
    pd.testing.assert_frame_equal(
        updated[original.columns],
        original.astype({"square_id": "string", "cell_id": "string"}),
    )
    assert updated.columns.tolist() == [
        *original.columns, "domain", "cell_type"
    ]

    # Square s1 is left of the boundary; its cell centre is right.
    assert updated["domain"].tolist() == ["0", "1", "1"]

    # Both squares belonging to cell 1 receive its single cell type.
    assert updated["cell_type"].iloc[:2].tolist() == ["B", "B"]
    assert pd.isna(updated["cell_type"].iloc[2])

    # No additional output CSV was created.
    assert list(tmp_path.iterdir()) == [path]
