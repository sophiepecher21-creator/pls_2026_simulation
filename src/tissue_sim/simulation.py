"""Run the Week 3 toy: two cell types in two vertical domains."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from tissue_sim.provenance import require_verified, stamp


def run(config_path, raw_dir, out_dir):
    """Verify supplied data, place toy cells, and record provenance."""
    config_path = Path(config_path)
    raw_dir = Path(raw_dir)
    out_dir = Path(out_dir)

    # The manifest is next to raw/, outside the read-only directory.
    manifest_path = raw_dir.parent / "raw_manifest.sha256"
    require_verified(raw_dir, manifest_path)

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    seed = config["seed"]
    n_cells = config["n_cells"]
    width = config["field_width_um"]
    height = config["field_height_um"]
    cell_types = config["cell_types"]

    if not isinstance(n_cells, int) or n_cells <= 0:
        raise ValueError("n_cells must be a positive integer.")
    if width <= 0 or height <= 0:
        raise ValueError("Field dimensions must be positive.")
    if len(cell_types) != 2:
        raise ValueError("The Week 3 toy requires exactly two cell types.")

    rng = np.random.default_rng(seed)
    x = rng.uniform(0, width, n_cells)
    y = rng.uniform(0, height, n_cells)

    # Left half is domain 0; right half is domain 1.
    domains = (x >= width / 2).astype(int)
    types = np.asarray(cell_types)[domains]

    cells = pd.DataFrame(
        {
            "cell_id": [f"cell_{i:04d}" for i in range(n_cells)],
            "x_um": x,
            "y_um": y,
            "domain": domains,
            "cell_type": types,
        }
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    cells.to_csv(out_dir / "toy_cells.csv", index=False)

    metadata = stamp(config_path, seed, manifest_path)
    (out_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    return out_dir / "toy_cells.csv"


def main():
    """Read terminal arguments and run the toy."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    output = run(args.config, args.data_dir, args.out_dir)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
