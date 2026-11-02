"""Check conformance with generated fixtures independent of the course pack."""

from dataclasses import replace
import json

import numpy as np
import pandas as pd
import pytest

from tissue_sim.configuration import SimulationConfig
from tissue_sim.reference import generate_reference
from tissue_sim.simulation import simulate, write_results
from tissue_sim.validate import validate_dataset


@pytest.fixture
def dataset(tmp_path):
    config = SimulationConfig(
        n_cells=4,
        field_width_um=16.0,
        field_height_um=16.0,
        cell_radius_um=1.0,
        n_genes=6,
        n_marker_genes=2,
    )
    output = tmp_path / "dataset"
    write_results(simulate(config), output)
    return output


def test_generated_dataset_conforms(dataset, capsys):
    assert validate_dataset(dataset) is None
    assert capsys.readouterr() == ("", "")


def test_empty_composition_is_zero_not_missing(tmp_path):
    config = replace(SimulationConfig(), mean_transcripts_per_type=(0.0, 0.0))
    output = tmp_path / "zeros"
    write_results(simulate(config), output)
    validate_dataset(output)
    assert pd.read_csv(output / "ground_truth/composition.csv").empty


@pytest.mark.parametrize(
    "file,change,reason",
    [
        (
            "metadata.json",
            lambda m: m.update(contract_version="2.0"),
            "contract_version",
        ),
        ("metadata.json", lambda m: m.update(n_squares=3), "n_squares"),
        ("metadata.json", lambda m: m.update(square_size_um=4), "square_size_um"),
        ("ground_truth/params.json", lambda m: m.update(seed=99), "seed/cell_types"),
        (
            "ground_truth/params.json",
            lambda m: m.update(expression_programs=[[0]]),
            "expression_programs",
        ),
        (
            "ground_truth/params.json",
            lambda m: m.update(marker_genes={}),
            "marker_genes",
        ),
    ],
)
def test_invalid_json_values(dataset, file, change, reason):
    path = dataset / file
    value = json.loads(path.read_text())
    change(value)
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match=reason):
        validate_dataset(dataset)


@pytest.mark.parametrize(
    "case,reason",
    [
        ("duplicate_position", "without gaps or duplicates"),
        ("bad_coordinate", "x_um.*micrometres"),
        ("row_order", "square_id mismatch"),
        ("fractional_count", "non-negative integer"),
        ("missing_type", "reference must cover"),
        ("gene_order", "gene columns must match"),
        ("duplicate_pair", "duplicate.*pair"),
        ("unknown_square", "must occur in coordinates"),
        ("missing_domain", "square_id mismatch"),
    ],
)
def test_invalid_tables(dataset, case, reason):
    if case in {"duplicate_position", "bad_coordinate"}:
        path = dataset / "coordinates.csv"
        table = pd.read_csv(path)
        if case == "duplicate_position":
            table.loc[1, ["array_row", "array_col"]] = table.loc[
                0, ["array_row", "array_col"]
            ].values
        else:
            table.loc[0, "x_um"] = 1.5
    elif case in {"row_order", "fractional_count"}:
        path = dataset / "counts.csv"
        table = pd.read_csv(path, dtype=str)
        if case == "row_order":
            table = table.iloc[::-1]
        else:
            table.iloc[0, 1] = "1.5"
    elif case == "missing_type":
        path = dataset / "reference/reference_labels.csv"
        table = pd.read_csv(path)
        table["cell_type"] = "type_a"
    elif case == "gene_order":
        path = dataset / "reference/reference_counts.csv"
        table = pd.read_csv(path)
        table = table[[table.columns[0], *table.columns[1:][::-1]]]
    elif case == "missing_domain":
        path = dataset / "ground_truth/domains.csv"
        table = pd.read_csv(path).iloc[1:]
    else:
        path = dataset / "ground_truth/composition.csv"
        table = pd.read_csv(path)
        if case == "duplicate_pair":
            table = pd.concat([table, table.iloc[:1]])
        else:
            table.loc[0, "square_id"] = "unknown_square"
    table.to_csv(path, index=False)
    with pytest.raises(ValueError, match=reason):
        validate_dataset(dataset)


def test_missing_and_malformed_json(dataset):
    path = dataset / "metadata.json"
    path.write_text("{bad")
    with pytest.raises(ValueError, match="metadata.json.*readable JSON"):
        validate_dataset(dataset)
    path.unlink()
    with pytest.raises(ValueError, match="metadata.json.*readable JSON"):
        validate_dataset(dataset)


def test_reference_known_answer_and_fresh_sampling():
    programs = {"A": np.array([1.0, 0.0]), "B": np.array([0.0, 1.0])}
    args = (programs, ["g1", "g2"], {"A": 1000, "B": 2000}, 80)
    a, labels = generate_reference(*args, np.random.default_rng(42))
    b, _ = generate_reference(*args, np.random.default_rng(42))
    pd.testing.assert_frame_equal(a, b)
    assert labels.cell_type.value_counts().to_dict() == {"A": 80, "B": 80}
    assert a.loc[:79, "g2"].eq(0).all()
    assert a.loc[80:, "g1"].eq(0).all()
    assert a.loc[:79, "g1"].nunique() > 1
    assert abs(a.loc[80:, "g2"].mean() - 2000) < 30


def test_notebook_layer_checks(tmp_path):
    """The three notebook exercises work on its deliberately partial 3x3 fixture."""
    import csv

    from tissue_sim.validate import (
        check_composition_sums,
        check_counts,
        check_ids_match,
        read_csv,
    )

    def write(name, rows):
        with (tmp_path / name).open("w", newline="") as stream:
            csv.writer(stream).writerows(rows)

    (tmp_path / "ground_truth").mkdir()
    ids = [f"sq_r{r}_c{c}" for r in range(3) for c in range(3)]
    counts = [["square_id", "GENE_1", "GENE_2", "GENE_3"]]
    coords = [["square_id", "array_row", "array_col", "x_um", "y_um"]]
    comp = [["square_id", "cell_type", "n_transcripts"]]
    for i, square in enumerate(ids):
        r, c = divmod(i, 3)
        counts.append([square, i, 0, 1])
        coords.append([square, r, c, 2 * c + 1, 2 * r + 1])
        total = i + 1
        a, b = (
            (total, 0)
            if c == 0
            else ((0, total) if c == 2 else (total - total // 2, total // 2))
        )
        comp.extend([[square, "A", a], [square, "B", b]])
    write("counts.csv", counts)
    write("coordinates.csv", coords)
    write("ground_truth/composition.csv", comp)
    assert read_csv(tmp_path / "counts.csv")[0]["GENE_1"] == "0"
    assert check_counts(tmp_path) is None
    assert check_ids_match(tmp_path) is None
    assert check_composition_sums(tmp_path) is None

    counts[4][2] = -6
    write("counts.csv", counts)
    with pytest.raises(ValueError, match=r"counts.csv.*sq_r1_c0.*GENE_2='-6'"):
        check_counts(tmp_path)
    counts[4][2] = 0
    write("counts.csv", counts)

    coords.append(["sq_r9_c9", 9, 9, 19, 19])
    write("coordinates.csv", coords)
    with pytest.raises(ValueError, match="sq_r9_c9"):
        check_ids_match(tmp_path)
    write("coordinates.csv", coords[:-1])

    comp[-2][2], comp[-1][2] = 9, 3
    write("ground_truth/composition.csv", comp)
    with pytest.raises(ValueError, match=r"sq_r2_c2.*sum=12.*total=9"):
        check_composition_sums(tmp_path)
