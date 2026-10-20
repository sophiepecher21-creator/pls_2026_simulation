import argparse
import csv
from pathlib import Path

import numpy as np
import yaml
from scipy.spatial import cKDTree

DEFAULT_CONFIG = "config/default.yaml"


# ---------------------------------------------------------------------------
# config + csv helpers
# ---------------------------------------------------------------------------

def _find_key(data, key):
    """Find a key anywhere in a nested dict (config values may be nested)."""
    if isinstance(data, dict):
        if key in data:
            return data[key]
        for value in data.values():
            found = _find_key(value, key)
            if found is not None:
                return found
    return None


def load_config(config_path=DEFAULT_CONFIG):
    """Load the YAML config file into a dict."""
    config_path = Path(config_path)
    if not config_path.is_file():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with config_path.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def _require(config, key):
    value = _find_key(config, key)
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
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# cell placement
# ---------------------------------------------------------------------------

def _number_of_cells(config):
    """Number of cells: 'n_cells' if set, otherwise density * field area."""
    n_cells = _find_key(config, "n_cells")
    if n_cells is None:
        density = _find_key(config, "cell_density_per_um2")
        if density is None:
            raise KeyError("Give either 'n_cells' or 'cell_density_per_um2' in the config file.")
        area = float(_require(config, "field_width_um")) * float(_require(config, "field_height_um"))
        n_cells = round(float(density) * area)
    n_cells = int(n_cells)
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
        n_cells            : concrete number of cells (optional; if set, it
                             overrides the density)
        cell_density_per_um2 : cells per um^2; N = density * width * height
                             (used only when n_cells is not set)
        distribution_shape : "uniform", "gaussian" or "clustered"
        seed               : integer or null (null = different every run)

    coordinates : list of dicts (one per measurement square) with the columns
                  "x_um" and "y_um", e.g. from read_coordinates().

    Returns (coordinates, cells):
      coordinates : copy of the records with a new column "cell_id"
                    (1..N, or 0 if the square center lies in no cell)
      cells       : list of dicts with cell_id, center_x_um, center_y_um,
                    radius_um
    """
    config = load_config(config_path)
    field_width_um = float(_require(config, "field_width_um"))
    field_height_um = float(_require(config, "field_height_um"))
    cell_radius = float(_require(config, "cell_radius_um"))
    n_cells = _number_of_cells(config)
    distribution_shape = str(_require(config, "distribution_shape"))
    if distribution_shape not in ("uniform", "gaussian", "clustered"):
        raise ValueError(f"Unknown distribution_shape: {distribution_shape!r}")
    seed = _find_key(config, "seed")

    rng = np.random.default_rng(seed)

    records = [dict(r) for r in coordinates]
    pts = np.array([[float(r["x_um"]), float(r["y_um"])] for r in records])

    xmin, xmax = 0.0, field_width_um
    ymin, ymax = 0.0, field_height_um
    width, height = xmax - xmin, ymax - ymin
    if width <= 0 or height <= 0:
        raise ValueError("field_width_um and field_height_um must be positive.")

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

    for record, cell_id in zip(records, cell_ids):
        record["cell_id"] = int(cell_id)

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
    p = argparse.ArgumentParser(description="Place round cells on the measurement grid.")
    p.add_argument("--input", required=True, help="grid CSV (coordinates.csv)")
    p.add_argument("--config", default=DEFAULT_CONFIG,
                   help=f"YAML config file (default: {DEFAULT_CONFIG})")
    p.add_argument("--output", help="CSV with cell_id column "
                                    "(default: <input folder>/coordinates_with_cells.csv)")
    p.add_argument("--cells-output", help="CSV with the cell table "
                                          "(default: <input folder>/cells.csv)")
    args = p.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path.with_name("coordinates_with_cells.csv")
    cells_path = Path(args.cells_output) if args.cells_output else in_path.with_name("cells.csv")

    coordinates = read_coordinates(in_path)
    records, cells = place_cells(coordinates, args.config)

    _write_csv(out_path, records)
    if cells:
        _write_csv(cells_path, cells)

    requested = _number_of_cells(load_config(args.config))
    if len(cells) < requested:
        print(f"Warning: only placed {len(cells)} of {requested} cells "
              f"(field too crowded for the cell radius).")
    print(f"cells placed: {len(cells)}")
    print(f"squares inside a cell: {sum(1 for r in records if r['cell_id'])} of {len(records)}")
    print(f"wrote {out_path}")
    print(f"wrote {cells_path}")


if __name__ == "__main__":
    main()