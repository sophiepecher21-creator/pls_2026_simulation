import argparse
from pathlib import Path
import csv
import yaml
import math
from numbers import Real
from dataclasses import dataclass

@dataclass(frozen=True)
class GridConfig:
    """Validated settings used to construct the grid."""
    field_width_um: float
    field_height_um: float
    square_size_um: float

def read_config(path):
    """Read and validate the grid settings from a YAML configuration."""
    with Path(path).open(encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    if not isinstance(config, dict):
        raise ValueError("The configuration must contain a YAML mapping.")

    required = ("field_width_um","field_height_um","square_size_um",)

    missing = [name for name in required if name not in config]
    if missing:
        raise ValueError(
            f"Missing configuration settings: {', '.join(missing)}"
        )

    for name in required:
        value = config[name]

        if (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(
                f"{name} must be a positive, finite number."
            )

    return GridConfig(field_width_um=config["field_width_um"],field_height_um=config["field_height_um"],square_size_um=config["square_size_um"])

def make_grid(field_width_um, field_height_um, square_size_um):
    """Return measurement-square IDs, indices, and centre coordinates.

    Field dimensions must be positive, finite multiples of square_size_um.
    Row and column indices start at zero.
    """
    dimensions = (field_width_um, field_height_um, square_size_um)

    # Reject invalid dimensions before calculating the grid.
    for value in dimensions:
        if (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(
                "Field dimensions and square size must be positive, finite numbers."
            )

    col_count = field_width_um / square_size_um
    row_count = field_height_um / square_size_um

    if col_count != round(col_count) or row_count != round(row_count):
    	raise ValueError("Field dimensions must be multiples of square_size_um.")

    n_cols = round(col_count)
    n_rows = round(row_count)
    grid = {}

    for row in range(n_rows):
        for col in range(n_cols):
            square_id = f"sq_r{row:03d}_c{col:03d}"

            grid[square_id] = {
                "array_row": row,
                "array_col": col,
                "x_um": (col + 0.5) * square_size_um,
                "y_um": (row + 0.5) * square_size_um,
            }

    return grid

def write_coordinates(grid, path):
    """Write grid coordinates to a new CSV, refusing to overwrite a file."""
    with open(path, "x", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["square_id", "array_row", "array_col", "x_um", "y_um"])

        for square_id, square in grid.items():
            writer.writerow([square_id,square["array_row"],square["array_col"],square["x_um"],square["y_um"]])

def main():
    """Read terminal arguments and generate the coordinates file."""
    parser = argparse.ArgumentParser(
        description="Generate a measurement-grid coordinates CSV."
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    config = read_config(args.config)

    grid = make_grid(
        field_width_um=config.field_width_um,
        field_height_um=config.field_height_um,
        square_size_um=config.square_size_um,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_coordinates(grid, args.output)

    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
