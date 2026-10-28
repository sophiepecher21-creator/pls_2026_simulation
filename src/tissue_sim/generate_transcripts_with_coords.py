""" Step 5: Transcripts in Space

What this script does:
    For every cell in the tissue:
        1. Draw a random number of transcripts  (Poisson distribution)
        2. For each transcript, draw a gene     (from the expression program p_t)
        3. Give it a random position inside the cell disk
        4. Exclude any transcript outside the disk 
        5. Record which 2-µm square the transcript landed in

Inputs:
    coordinates.csv         — grid squares + cell geometry (center, radius)
    toy_cells.csv           — which cell_id has which cell_type
    expression_programs.json — probability distribution over genes, per cell type

Outputs:
    transcripts.csv         — one row per transcript (position, gene, cell)
    composition.csv         — one row per square, with transcript count

Usage:
    python generate_transcripts.py \
        --coordinates  coordinates.csv \
        --cells        toy_cells.csv \
        --programs     expression_programs.json \
        --out-dir      output/ \
        --seed         42
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

def load_coordinates(path):
    """
    Read coordinates.csv.
    Returns a dict: square_id -> {cell_id, center_x_um, center_y_um, radius_um}
    We only need the cell geometry here (center + radius), not the square position,
    because transcript coordinates are generated from the cell, not the square.
    """
    coords = {}
    with open(path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            cell_id = int(row["cell_id"])
            coords[row["square_id"]] = {
                "cell_id":     cell_id,
                # cell geometry — empty string if background square (cell_id == 0)
                "center_x_um": float(row["center_x_um"]) if row["center_x_um"] else None,
                "center_y_um": float(row["center_y_um"]) if row["center_y_um"] else None,
                "radius_um":   float(row["radius_um"])   if row["radius_um"]   else None,
            }
    return coords


def load_cell_types(path):
    """
    Read toy_cells.csv (produced by simulation.py).
    Returns a dict: cell_id (int) -> cell_type (str), e.g. {1: "A", 2: "B"}
    """
    cell_types = {}
    with open(path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            # cell_id may be an integer or a string like "cell_0001"
            # strip any non-numeric characters to always get a plain integer
            raw = row["cell_id"]
            cid = int("".join(ch for ch in raw if ch.isdigit()))
            cell_types[cid] = row["cell_type"]
    return cell_types


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
    gene_names       = data.pop("gene_names", None)
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
    r     = radius * np.sqrt(rng.uniform(0.0, 1.0, n)) # random distance from center, 0 to R
    # random angle: full circle
    theta = rng.uniform(0.0, 2.0 * math.pi, n)

    # convert polar -> cartesian, offset by cell center
    x_t = r * np.cos(theta) + cx        # cx is the cell center
    y_t = r * np.sin(theta) + cy        # cy is the cell center

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
    col = int(x / square_size)
    row = int(y / square_size)
    return f"sq_r{row:03d}_c{col:03d}"


# =============================================================================
# Step 4 - main loop: generate all transcripts
# =============================================================================

def generate_transcripts(coords, cell_types, programs, gene_names,
                         mean_transcripts, default_mean, square_size, rng):
    """
    Loop over every cell and generate its transcripts.
    Returns:
        transcripts  — list of dicts, one per transcript
        composition  — list of dicts, one per (square, cell_type) with counts
    """

    # --- collect unique cell geometries from coordinates.csv ----------------
    # multiple squares can belong to the same cell, but the geometry is the same
    # so we just need one entry per cell_id
    cell_geometry = {}
    for sq_info in coords.values():
        cid = sq_info["cell_id"]
        if cid == 0:
            continue  # skip background squares (not inside any cell)
        if cid not in cell_geometry:
            cell_geometry[cid] = {
                "center_x_um": sq_info["center_x_um"],
                "center_y_um": sq_info["center_y_um"],
                "radius_um":   sq_info["radius_um"],
            }

    transcripts = []   # will hold one dict per transcript
    sq_counts   = {}   # (square_id, cell_id, cell_type) -> count, for composition.csv
    transcript_id = 0  # global counter across all transcripts

    # --- loop over every cell ------------------------------------------------
    for cid, geom in sorted(cell_geometry.items()):

        # get the cell type for this cell (e.g. "A" or "B")
        cell_type = cell_types.get(cid)
        if cell_type is None:
            print(f"Warning: cell_id={cid} has no cell_type, skipping.")
            continue

        # get the expression program (probability distribution over genes)
        p_t = programs.get(cell_type)
        if p_t is None:
            print(f"Warning: cell_type='{cell_type}' has no expression program, skipping.")
            continue

        # unpack cell geometry
        cx     = geom["center_x_um"]
        cy     = geom["center_y_um"]
        radius = geom["radius_um"]

        # ----- B1: draw total number of transcripts for this cell ------------
        # Poisson distribution: lam is the mean number of transcripts per cell
        # each cell type can have a different mean (set in the JSON)
        lam     = mean_transcripts.get(cell_type, default_mean)
        n_total = int(rng.poisson(lam))

        if n_total == 0:
            continue  # this cell happened to get zero transcripts this run

        # ----- B2: draw a gene for each transcript ---------------------------
        # rng.choice picks n_total gene indices at random, weighted by p_t
        gene_indices = rng.choice(len(gene_names), size=n_total, p=p_t)

        # ----- B3: place each transcript inside the cell disk ----------------
        x_t, y_t = place_transcripts_in_disk(n_total, cx, cy, radius, rng)

        # ----- B4: exclude any transcript outside the disk (safety check) ----
        # with the sqrt method all points are inside by construction,
        # but we check anyway as required by the spec
        dist_squared = (x_t - cx)**2 + (y_t - cy)**2
        inside       = dist_squared <= radius**2
        x_t          = x_t[inside]
        y_t          = y_t[inside]
        gene_indices = gene_indices[inside]

        # ----- B5: record each transcript ------------------------------------
        for x, y, gi in zip(x_t, y_t, gene_indices):

            # find which 2-µm square this transcript landed in
            sq_id = get_square_id(float(x), float(y), square_size)
            gene  = gene_names[gi]

            # add to transcript list
            transcripts.append({
                "transcript_id": transcript_id,
                "square_id":     sq_id,
                "cell_id":       cid,
                "cell_type":     cell_type,
                "gene":          gene,
                "x_um":          round(float(x), 4),
                "y_um":          round(float(y), 4),
            })

            # increment count for this square (used to build composition.csv)
            key = (sq_id, cid, cell_type)
            sq_counts[key] = sq_counts.get(key, 0) + 1

            transcript_id += 1

    # --- build composition rows from the square counts -----------------------
    composition = [
        {
            "square_id":     sq_id,
            "cell_id":       cid,
            "cell_type":     ctype,
            "n_transcripts": count,
        }
        for (sq_id, cid, ctype), count in sorted(sq_counts.items())
    ]

    return transcripts, composition


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
    parser = argparse.ArgumentParser(
        description="Step 5 — generate transcripts in space for each cell."
    )
    parser.add_argument("--coordinates",      type=Path, required=True)
    parser.add_argument("--cells",            type=Path, required=True)
    parser.add_argument("--programs",         type=Path, required=True)
    parser.add_argument("--out-dir",          type=Path, default=Path("output"))
    parser.add_argument("--mean-transcripts", type=float, default=20.0,
                        help="Fallback Poisson mean if not in the JSON (default: 20)")
    parser.add_argument("--square-size",      type=float, default=2.0,
                        help="Size of one measurement square in µm (default: 2.0)")
    parser.add_argument("--seed",             type=int,   default=None,
                        help="Random seed for reproducibility")
    args = parser.parse_args()

    # load all three input files
    print(f"Loading coordinates    : {args.coordinates}")
    coords = load_coordinates(args.coordinates)

    print(f"Loading cell types     : {args.cells}")
    cell_types = load_cell_types(args.cells)

    print(f"Loading programs       : {args.programs}")
    programs, gene_names, mean_transcripts = load_expression_programs(args.programs)

    # create the random number generator (seed makes results reproducible)
    rng = np.random.default_rng(args.seed)

    # run the main generation
    print("\nGenerating transcripts …")
    transcripts, composition = generate_transcripts(
        coords           = coords,
        cell_types       = cell_types,
        programs         = programs,
        gene_names       = gene_names,
        mean_transcripts = mean_transcripts,
        default_mean     = args.mean_transcripts,
        square_size      = args.square_size,
        rng              = rng,
    )

    # write outputs
    write_csv(
        args.out_dir / "transcripts.csv",
        transcripts,
        ["transcript_id", "square_id", "cell_id", "cell_type", "gene", "x_um", "y_um"],
    )
    write_csv(
        args.out_dir / "composition.csv",
        composition,
        ["square_id", "cell_id", "cell_type", "n_transcripts"],
    )

    print(f"\nDone.")
    print(f"  Total transcripts : {len(transcripts)}")
    print(f"  Squares covered   : {len(composition)}")
    print(f"  transcripts.csv  -> {args.out_dir / 'transcripts.csv'}")
    print(f"  composition.csv  -> {args.out_dir / 'composition.csv'}")


if __name__ == "__main__":
    main()
