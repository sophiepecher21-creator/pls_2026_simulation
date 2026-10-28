import argparse
import csv
from pathlib import Path

import numpy as np
import yaml
from scipy.spatial import cKDTree

DEFAULT_CONFIG = "config/default.yaml"

# columns this script adds to (or overwrites in) the CSV
CELL_COLUMNS = ["cell_id", "center_x_um", "center_y_um", "radius_um"]


# ---------------------------------------------------------------------------
# config + csv helpers
# ---------------------------------------------------------------------------

def load_config(config_path=DEFAULT_CONFIG):
    """Load the YAML config file into a dict."""
    config_path = Path(config_path)
    if not config_path.is_file():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with config_path.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def _require(config, key):
    value = config.get(key)
    if value is None:
        raise KeyError(f"'{key}' is missing in the config file.")
    return value


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
            raise ValueError(f"{path} is missing column(s): {', '.join(sorted(missing))}")
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

def _number_of_cells(config):
    """Number of cells, read from 'n_cells' in the config file."""
    n_cells = int(_require(config, "n_cells"))
    if n_cells < 1:
        raise ValueError("The number of cells must be at least 1.")
    return n_cells


def place_cells(coordinates, config_path=DEFAULT_CONFIG):
    """
    Place round, non-overlapping cells on the field and label the grid squares.
    All parameters are read from the config file:

        field_width_um     : width of the whole field in um (x from 0 to width)
        field_height_um    : height of the whole field in um (y from 0 to height)
        cell_radius_um     : radius of each round cell in um
        n_cells            : number of cells to place
        distribution_shape : "uniform", "gaussian" or "clustered"
        seed               : integer or null (null = different every run)

    coordinates : list of dicts (one per measurement square) with the columns
                  "x_um" and "y_um", e.g. from read_coordinates().

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

    config = load_config(config_path)
    field_width_um = float(_require(config, "field_width_um"))
    field_height_um = float(_require(config, "field_height_um"))
    cell_radius = float(_require(config, "cell_radius_um"))
    n_cells = _number_of_cells(config)
    distribution_shape = str(_require(config, "distribution_shape"))
    if distribution_shape not in ("uniform", "gaussian", "clustered"):
        raise ValueError(f"Unknown distribution_shape: {distribution_shape!r}")
    seed = config.get("seed")

    rng = np.random.default_rng(seed)

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
        cluster_centers = np.column_stack([
            rng.uniform(xmin, xmax, k), rng.uniform(ymin, ymax, k)])

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
        if centers and np.min(np.linalg.norm(np.array(centers) - c, axis=1)) < 2 * cell_radius:
            continue
        centers.append(c)

    # --- assign every square center to the cell whose blob contains it ------
    cell_ids = np.zeros(len(pts), dtype=int)
    if centers:
        dist, nearest = cKDTree(np.array(centers)).query(pts)
        inside = dist <= cell_radius
        cell_ids[inside] = nearest[inside] + 1  # number cells 1..N

    # --- write the cell info into every row (empty for squares outside) -----
    for record, cell_id in zip(records, cell_ids):
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
    p = argparse.ArgumentParser(
        description="Place round cells on the measurement grid. The input CSV "
                    "is updated in place (columns: " + ", ".join(CELL_COLUMNS) + ").")
    p.add_argument("--input", required=True, help="grid CSV (coordinates.csv), edited in place")
    p.add_argument("--config", default=DEFAULT_CONFIG,
                   help=f"YAML config file (default: {DEFAULT_CONFIG})")
    args = p.parse_args()

    coordinates = read_coordinates(args.input)
    records, cells = place_cells(coordinates, args.config)
    _write_csv(args.input, records)

    requested = _number_of_cells(load_config(args.config))
    if len(cells) < requested:
        print(f"Warning: only placed {len(cells)} of {requested} cells "
              f"(field too crowded for the cell radius).")
    print(f"cells placed: {len(cells)}")
    print(f"squares inside a cell: {sum(1 for r in records if r['cell_id'])} of {len(records)}")
    print(f"updated {args.input}")


if __name__ == "__main__":
    main()