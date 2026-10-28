"""Generate positioned transcripts from complete typed cells and expression programs.

The pipeline calls generate_transcripts() directly with in-memory objects.
The optional standalone CLI reads diagnostics/cells.csv, coordinates.csv and
an expression-program JSON file, then writes transcripts and composition.
"""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


# =============================================================================
# Step 1 - load input files
# =============================================================================


def load_expression_programs(path):
    """
    Read expression_programs.json.
    Returns:
        programs        — dict: cell_type -> numpy array of gene probabilities (p_t)
        gene_names      — list of gene name strings
        mean_transcripts — dict: cell_type -> Poisson mean (may be empty)
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    # these two keys are metadata, not expression programs — pop them out first
    gene_names = data.pop("gene_names", None)
    mean_transcripts = data.pop("mean_transcripts", {})

    # everything remaining is a cell type -> probability list
    programs = {}
    for cell_type, probs in data.items():
        arr = np.array(probs, dtype=float)
        # normalise defensively so probabilities always sum exactly to 1
        programs[cell_type] = arr / arr.sum()

    # if no gene names were provided, just call them gene_0, gene_1, ...
    n_genes = len(next(iter(programs.values())))
    if gene_names is None:
        gene_names = [f"gene_{i}" for i in range(n_genes)]

    return programs, gene_names, mean_transcripts


# =============================================================================
# Step 2 - place transcripts inside a cell disk
# =============================================================================


def place_transcripts_in_disk(n, cx, cy, radius, rng):
    """
    Generate n random positions uniformly inside a disk.

    Uses the polar coordinate method from the checklist:
        r     = radius * sqrt(U),   U ~ Uniform(0, 1)
        theta ~ Uniform(0, 2*pi)
        x_t   = r * cos(theta) + cx
        y_t   = r * sin(theta) + cy

    Why sqrt(U) and not just U?
        Without it, points would cluster near the centre.
        The sqrt corrects for the fact that outer rings have more area.
    """
    # random radius: sqrt gives uniform distribution over disk area
    r = radius * np.sqrt(
        rng.uniform(0.0, 1.0, n)
    )  # random distance from center, 0 to R
    # random angle: full circle
    theta = rng.uniform(0.0, 2.0 * math.pi, n)

    # convert polar -> cartesian, offset by cell center
    x_t = r * np.cos(theta) + cx  # cx is the cell center
    y_t = r * np.sin(theta) + cy  # cy is the cell center

    return x_t, y_t


# =============================================================================
# Step 3 - find which square a transcript belongs to
# =============================================================================


def get_square_id(x, y, square_size):
    """
    Given a transcript's (x, y) position in µm, return the square_id
    of the 2-µm bin it falls into.

    Uses floor division — same convention as make_grid_csv.py.
    Example: x=5.3, square_size=2 -> col=2 -> sq_r???_c002
    """
    col = math.floor(x / square_size)
    row = math.floor(y / square_size)
    return f"sq_r{row:03d}_c{col:03d}"


# =============================================================================
# Step 4 - main loop: generate all transcripts
# =============================================================================


def generate_transcripts(
    cells,
    programs,
    gene_names,
    mean_transcripts,
    square_size,
    field_width_um,
    field_height_um,
    rng,
):
    """Return in-field transcripts from a complete typed cell DataFrame.

    Cells contain cell_id, cell_type, center_x_um, center_y_um, radius_um.
    Programs map type names to gene probabilities; mean_transcripts maps types
    to Poisson means. Sampling uses only the supplied shared Generator.
    Disks may cross the field boundary; only transcripts inside the field survive.
    Input tables are unchanged. The returned table retains its schema when empty.
    """
    import pandas as pd

    columns = [
        "transcript_id",
        "square_id",
        "cell_id",
        "cell_type",
        "gene",
        "x_um",
        "y_um",
    ]
    transcripts = []
    for cell in cells.sort_values("cell_id").itertuples(index=False):
        p_t = programs[cell.cell_type]
        n = int(rng.poisson(mean_transcripts[cell.cell_type]))
        if not n:
            continue
        genes = rng.choice(len(gene_names), size=n, p=p_t)
        x, y = place_transcripts_in_disk(
            n, cell.center_x_um, cell.center_y_um, cell.radius_um, rng
        )
        inside = (x - cell.center_x_um) ** 2 + (
            y - cell.center_y_um
        ) ** 2 <= cell.radius_um**2
        inside &= (x >= 0) & (x < field_width_um) & (y >= 0) & (y < field_height_um)
        for tx, ty, gi in zip(x[inside], y[inside], genes[inside], strict=True):
            transcripts.append(
                {
                    "transcript_id": len(transcripts),
                    "square_id": get_square_id(float(tx), float(ty), square_size),
                    "cell_id": cell.cell_id,
                    "cell_type": cell.cell_type,
                    "gene": gene_names[gi],
                    "x_um": float(tx),
                    "y_um": float(ty),
                }
            )
    return pd.DataFrame(transcripts, columns=columns)


# =============================================================================
# Step 5 - write output files
# =============================================================================


def write_csv(path, rows, fieldnames):
    """
    Write a list of dicts to a CSV file.
    Uses a temp file so a crash can never leave a half-written file behind.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(path).with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    # only rename to final name once writing is complete
    tmp.replace(path)


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================


def main():
    """Run transcript generation from explicitly supplied standalone input files."""
    parser = argparse.ArgumentParser(
        description="Step 5 — generate transcripts in space for each cell."
    )
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--programs", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("output"))
    parser.add_argument(
        "--mean-transcripts",
        type=float,
        default=20.0,
        help="Fallback Poisson mean if not in the JSON (default: 20)",
    )
    parser.add_argument(
        "--square-size",
        type=float,
        default=2.0,
        help="Size of one measurement square in µm (default: 2.0)",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility"
    )
    args = parser.parse_args()

    import pandas as pd

    coordinates = pd.read_csv(args.coordinates)
    cells = pd.read_csv(args.cells)
    required = {"cell_id", "cell_type", "center_x_um", "center_y_um", "radius_um"}
    if not required.issubset(cells.columns):
        raise ValueError(
            "--cells must be the complete typed cells table from simulation.py"
        )
    programs, gene_names, means = load_expression_programs(args.programs)
    means = {name: means.get(name, args.mean_transcripts) for name in programs}
    transcripts = generate_transcripts(
        cells,
        programs,
        gene_names,
        means,
        args.square_size,
        float(coordinates.x_um.max() + args.square_size / 2),
        float(coordinates.y_um.max() + args.square_size / 2),
        np.random.default_rng(args.seed),
    )
    from tissue_sim.aggregation import aggregate_transcripts

    _, composition = aggregate_transcripts(coordinates, transcripts, gene_names)
    write_csv(
        args.out_dir / "transcripts.csv",
        transcripts.to_dict("records"),
        list(transcripts.columns),
    )
    write_csv(
        args.out_dir / "composition.csv",
        composition.to_dict("records"),
        list(composition.columns),
    )


if __name__ == "__main__":
    main()
