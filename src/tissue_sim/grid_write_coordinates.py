import csv


def write_coordinates(grid, path):
    """Write grid coordinates to a new CSV, refusing to overwrite a file."""
    with open(path, "x", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["square_id", "array_row", "array_col", "x_um", "y_um"])

        for square_id, square in grid.items():
            writer.writerow([square_id,square["array_row"],square["array_col"],square["x_um"],square["y_um"]])
