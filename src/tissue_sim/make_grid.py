import math
from numbers import Real


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
