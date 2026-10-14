import csv
import pytest

from tissue_sim.make_grid_csv import (GridConfig,main,make_grid,read_config,write_coordinates)


def test_make_grid_known_answer():
	''' A 4x4 µm grid matches four squares from the reference csv '''
	expected = {
       	 "sq_r000_c000": {"array_row": 0,"array_col": 0,"x_um": 1.0,"y_um": 1.0},
       	 "sq_r000_c001": {"array_row": 0,"array_col": 1,"x_um": 3.0,"y_um": 1.0},
       	 "sq_r001_c000": {"array_row": 1,"array_col": 0,"x_um": 1.0,"y_um": 3.0},
       	 "sq_r001_c001": {"array_row": 1,"array_col": 1,"x_um": 3.0,"y_um": 3.0}
    	}
	got = make_grid(field_width_um=4,field_height_um=4,square_size_um=2)
	assert got == expected

# Include a second parameterization test to validate the inputs - no negative numbers, etc. allowed

@pytest.mark.parametrize(
    "width, height, square_size",
    [
        (0, 4, 2),
        (4, -1, 2),
        (4, 4, 0),
        (5, 4, 2),
        (float("nan"), 4, 2),
        (float("inf"), 4, 2),
    ],
)
def test_make_grid_rejects_invalid_dimensions(width, height, square_size):
    with pytest.raises(ValueError):
        make_grid(width, height, square_size)

def test_read_config(tmp_path):
    """The reader returns validated settings from the YAML file."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "field_width_um: 4\n"
        "field_height_um: 2\n"
        "square_size_um: 2\n",
        encoding="utf-8",
    )

    assert read_config(config_path) == GridConfig(
        field_width_um=4,
        field_height_um=2,
        square_size_um=2,
    )


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


def test_main_known_answer(tmp_path, monkeypatch, capsys):
    """The terminal workflow writes the expected grid."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "field_width_um: 4\n"
        "field_height_um: 2\n"
        "square_size_um: 2\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "new_results" / "coordinates.csv"

# Monkeypatch temporarily supplies the terminal arguments that main() reads. Pytest will restore them afterwards
    monkeypatch.setattr(
        "sys.argv",
        [
            "make_grid_csv",
            "--config",
            str(config_path),
            "--output",
            str(output_path),
        ],
    )

    main()

    with output_path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))

    assert rows == [
        ["square_id", "array_row", "array_col", "x_um", "y_um"],
        ["sq_r000_c000", "0", "0", "1.0", "1.0"],
        ["sq_r000_c001", "0", "1", "3.0", "1.0"],
    ]
# Capsys captures printed message so we can checm it. All of the files are created in tmp_path
    assert capsys.readouterr().out == f"Wrote {output_path}\n"
