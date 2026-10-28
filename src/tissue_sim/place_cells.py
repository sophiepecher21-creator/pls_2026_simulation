import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

DEFAULT_CONFIG = "config/default.yaml"

# columns this script adds to (or overwrites in) the CSV
CELL_COLUMNS = ["cell_id", "center_x_um", "center_y_um", "radius_um"]


# ---------------------------------------------------------------------------
# config + csv helpers
# ---------------------------------------------------------------------------


def read_coordinates(path):
    """
    Read the grid CSV (square_id, array_row, array_col, x_um, y_um) and return
    a list of dicts, one per measurement square.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"CSV not found: {path}")
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        missing = {"x_um", "y_um"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path} is missing column(s): {', '.join(sorted(missing))}"
            )
        return list(reader)


def _write_csv(path, rows):
    """Write rows to `path`, replacing the file (via a temp file, so a crash
    cannot leave a half-written CSV behind)."""
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# cell placement
# ---------------------------------------------------------------------------


def place_cells(coordinates, config, rng):
    """
    Place round, non-overlapping cells on the field and label the grid squares.
    Parameters are supplied by a SimulationConfig object:

        field_width_um     : width of the whole field in um (x from 0 to width)
        field_height_um    : height of the whole field in um (y from 0 to height)
        cell_radius_um     : radius of each round cell in um
        n_cells            : number of cells to place
        distribution_shape : "uniform", "gaussian" or "clustered"
        rng                : shared NumPy random generator supplied by the caller

    coordinates : list of dicts (one per measurement square) with the columns
                  "x_um" and "y_um", with numeric values. config is a SimulationConfig; rng is the shared Generator.

    Returns (records, cells):
      records : copy of the rows (all original columns kept, e.g. square_id,
                array_row, array_col, x_um, y_um) plus cell_id (1..N, or 0 if
                the square center lies in no cell), center_x_um, center_y_um
                and radius_um (empty if cell_id is 0)
      cells   : list of dicts with cell_id, center_x_um, center_y_um,
                radius_um
    """
    if not coordinates:
        raise ValueError("No coordinates in input.")

    field_width_um = config.field_width_um
    field_height_um = config.field_height_um
    cell_radius = config.cell_radius_um
    n_cells = config.n_cells
    distribution_shape = config.distribution_shape
    if distribution_shape not in ("uniform", "gaussian", "clustered"):
        raise ValueError(f"Unknown distribution_shape: {distribution_shape!r}")

    records = [dict(r) for r in coordinates]
    pts = np.array([[float(r["x_um"]), float(r["y_um"])] for r in records])

    xmin, xmax = 0.0, field_width_um
    ymin, ymax = 0.0, field_height_um
    width, height = xmax - xmin, ymax - ymin
    if width <= 0 or height <= 0:
        raise ValueError("field_width_um and field_height_um must be positive.")
    if cell_radius <= 0:
        raise ValueError("cell_radius_um must be positive.")

    cx, cy = xmin + width / 2, ymin + height / 2
    if distribution_shape == "clustered":
        k = max(2, int(np.sqrt(n_cells) / 2))
        cluster_centers = np.column_stack(
            [rng.uniform(xmin, xmax, k), rng.uniform(ymin, ymax, k)]
        )

    # --- rejection sampling: draw a candidate center, keep it only if it lies
    #     inside the field and does not overlap an already placed cell ---------
    centers = []
    attempts, max_attempts = 0, n_cells * 200
    while len(centers) < n_cells and attempts < max_attempts:
        attempts += 1

        if distribution_shape == "uniform":
            c = np.array([rng.uniform(xmin, xmax), rng.uniform(ymin, ymax)])
        elif distribution_shape == "gaussian":
            c = np.array([rng.normal(cx, width / 6), rng.normal(cy, height / 6)])
        else:  # "clustered"
            c = cluster_centers[rng.integers(len(cluster_centers))]
            c = c + rng.normal(0, [width / 15, height / 15])

        if not (xmin <= c[0] <= xmax and ymin <= c[1] <= ymax):
            continue
        if (
            centers
            and np.min(np.linalg.norm(np.array(centers) - c, axis=1)) < 2 * cell_radius
        ):
            continue
        centers.append(c)

    # --- assign every square center to the cell whose blob contains it ------
    cell_ids = np.zeros(len(pts), dtype=int)
    if centers:
        dist, nearest = cKDTree(np.array(centers)).query(pts)
        inside = dist <= cell_radius
        cell_ids[inside] = nearest[inside] + 1  # number cells 1..N

    # --- write the cell info into every row (empty for squares outside) -----
    for record, cell_id in zip(records, cell_ids, strict=True):
        cell_id = int(cell_id)
        record["cell_id"] = cell_id
        if cell_id:
            center = centers[cell_id - 1]
            record["center_x_um"] = float(center[0])
            record["center_y_um"] = float(center[1])
            record["radius_um"] = cell_radius
        else:
            record["center_x_um"] = ""
            record["center_y_um"] = ""
            record["radius_um"] = ""

    cells = [
        {
            "cell_id": i,
            "center_x_um": float(c[0]),
            "center_y_um": float(c[1]),
            "radius_um": cell_radius,
        }
        for i, c in enumerate(centers, start=1)
    ]
    return records, cells


# ---------------------------------------------------------------------------
# command line interface
# ---------------------------------------------------------------------------


def main():
    """Run the standalone cell-placement command on a grid CSV."""
    p = argparse.ArgumentParser(
        description="Place round cells on the measurement grid. The input CSV "
        "is updated in place (columns: " + ", ".join(CELL_COLUMNS) + ")."
    )
    p.add_argument(
        "--input", required=True, help="grid CSV (coordinates.csv), edited in place"
    )
    p.add_argument(
        "--config",
        default=DEFAULT_CONFIG,
        help=f"YAML config file (default: {DEFAULT_CONFIG})",
    )
    args = p.parse_args()

    from tissue_sim.configuration import read_simulation_config

    config = read_simulation_config(args.config)
    coordinates = read_coordinates(args.input)
    records, cells = place_cells(
        coordinates, config, np.random.default_rng(config.seed)
    )
    _write_csv(args.input, records)

    requested = config.n_cells
    if len(cells) < requested:
        print(
            f"Warning: only placed {len(cells)} of {requested} cells "
            f"(field too crowded for the cell radius)."
        )
    print(f"cells placed: {len(cells)}")
    print(
        f"squares inside a cell: {sum(1 for r in records if r['cell_id'])} of {len(records)}"
    )
    print(f"updated {args.input}")


if __name__ == "__main__":
    main()
