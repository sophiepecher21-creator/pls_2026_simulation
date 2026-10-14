import csv

import pytest

from tissue_sim.grid_write_coordinates import write_coordinates
from tissue_sim.make_grid_csv import make_grid_csv, read_config


def test_read_config(tmp_path):
    """The reader returns the values stored in the YAML file."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "field_width_um: 4\n"
        "field_height_um: 2\n"
        "square_size_um: 2\n",
        encoding="utf-8",
    )

    assert read_config(config_path) == {
        "field_width_um": 4,
        "field_height_um": 2,
        "square_size_um": 2,
    }


def test_write_coordinates_known_answer(tmp_path):
    """The writer produces the expected CSV columns and values."""
    grid = {
        "sq_r000_c000": {
            "array_row": 0,
            "array_col": 0,
            "x_um": 1.0,
            "y_um": 1.0,
        },
        "sq_r000_c001": {
            "array_row": 0,
            "array_col": 1,
            "x_um": 3.0,
            "y_um": 1.0,
        },
    }
    output_path = tmp_path / "coordinates.csv"

    write_coordinates(grid, output_path)

    with output_path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))

    assert rows == [
        ["square_id", "array_row", "array_col", "x_um", "y_um"],
        ["sq_r000_c000", "0", "0", "1.0", "1.0"],
        ["sq_r000_c001", "0", "1", "3.0", "1.0"],
    ]


def test_write_coordinates_preserves_existing_file(tmp_path):
    """An existing output is protected even if the new grid is different."""
    output_path = tmp_path / "coordinates.csv"
    original = b"existing results\n"
    output_path.write_bytes(original)

    with pytest.raises(FileExistsError):
        write_coordinates({}, output_path)

    assert output_path.read_bytes() == original


def test_make_grid_csv_known_answer(tmp_path):
    """The full workflow creates a folder and writes the expected grid."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "field_width_um: 4\n"
        "field_height_um: 2\n"
        "square_size_um: 2\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "new_results" / "coordinates.csv"

    returned_path = make_grid_csv(config_path, output_path)

    assert returned_path == output_path

    with output_path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))

    assert rows == [
        ["square_id", "array_row", "array_col", "x_um", "y_um"],
        ["sq_r000_c000", "0", "0", "1.0", "1.0"],
        ["sq_r000_c001", "0", "1", "3.0", "1.0"],
    ]
