import numpy as np
from scipy.spatial import cKDTree

import reader_csv

coordinates = reader_csv.read("TestCSV_path")


def place_cells(coordinates, distribution_shape="uniform", Cell_density=0.01,
                cell_radius=None, seed=None):
    """
    Place round, non-overlapping cells on a set of coordinates.

    coordinates        : list/array of rows from the CSV. Row 3 (index 2) = x,
                         row 4 (index 3) = y.
    distribution_shape : how cell centers are distributed over the area:
                         "uniform"   - evenly random over the bounding box
                         "gaussian"  - concentrated around the middle
                         "clustered" - several dense groups
    Cell_density       : cells per unit area (units of the coordinates squared).
                         Number of cells N = Cell_density * area of bounding box.
    cell_radius        : radius of each round blob. Default is derived from the
                         density so that cells fit without touching.
    seed               : optional random seed for reproducible placement.

    Returns the coordinates with a 5th row appended: for each coordinate the
    cell number (1..N) it belongs to, or 0 if it lies in no cell.
    """
    rng = np.random.default_rng(seed)

    rows = [list(r) for r in coordinates]
    x = np.asarray(rows[2], dtype=float)
    y = np.asarray(rows[3], dtype=float)
    pts = np.column_stack([x, y])

    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    width, height = xmax - xmin, ymax - ymin
    if width == 0 or height == 0:
        raise ValueError("Coordinates have no extent in x or y.")

    n_cells = max(1, int(round(Cell_density * width * height)))
    if cell_radius is None:
        cell_radius = 0.25 / np.sqrt(Cell_density)

    # --- sample candidate cell centers according to the chosen shape --------
    cx, cy = xmin + width / 2, ymin + height / 2
    if distribution_shape == "clustered":
        k = max(2, int(np.sqrt(n_cells) / 2))
        cluster_centers = np.column_stack([
            rng.uniform(xmin, xmax, k), rng.uniform(ymin, ymax, k)])

    def sample_center():
        if distribution_shape == "uniform":
            return np.array([rng.uniform(xmin, xmax), rng.uniform(ymin, ymax)])
        if distribution_shape == "gaussian":
            return np.array([rng.normal(cx, width / 6), rng.normal(cy, height / 6)])
        if distribution_shape == "clustered":
            c = cluster_centers[rng.integers(len(cluster_centers))]
            return c + rng.normal(0, [width / 15, height / 15])
        raise ValueError(f"Unknown distribution_shape: {distribution_shape!r}")

    # --- rejection sampling so that cells never overlap ---------------------
    centers = []
    attempts, max_attempts = 0, n_cells * 200
    while len(centers) < n_cells and attempts < max_attempts:
        attempts += 1
        c = sample_center()
        if not (xmin <= c[0] <= xmax and ymin <= c[1] <= ymax):
            continue
        if centers and np.min(np.linalg.norm(np.array(centers) - c, axis=1)) < 2 * cell_radius:
            continue
        centers.append(c)

    if len(centers) < n_cells:
        print(f"Warning: only placed {len(centers)} of {n_cells} cells "
              f"(area too crowded for radius {cell_radius:.3g}).")

    # --- assign every coordinate to the cell whose blob contains it ---------
    cell_ids = np.zeros(len(pts), dtype=int)
    if centers:
        dist, nearest = cKDTree(np.array(centers)).query(pts)
        inside = dist <= cell_radius
        cell_ids[inside] = nearest[inside] + 1  # number cells 1..N

    # --- append as the 5th row ---------------------------------------------
    if len(rows) > 4:
        rows[4] = cell_ids.tolist()
    else:
        rows.append(cell_ids.tolist())
    return rows


if __name__ == "__main__":
    result = place_cells(coordinates, "uniform", Cell_density=0.01, seed=42)
    print("cells placed:", max(result[4]), "| points inside cells:",
          sum(1 for v in result[4] if v))