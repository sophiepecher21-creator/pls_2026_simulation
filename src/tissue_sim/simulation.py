"""Connect the simulation stages in memory and export their final tables."""

import argparse
import sys
import datetime
from importlib.metadata import version
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import hashlib
import platform
import sys

import numpy as np
import pandas as pd

from tissue_sim.aggregation import aggregate_transcripts
from tissue_sim.configuration import DEFAULT_CONFIG_PATH, read_simulation_config
from tissue_sim.domains import assign_cell_types, build_domains_dataframe
from tissue_sim.expression_programs import (
    generate_cell_type_programs,
    generate_marker_genes,
    make_baseline_distribution,
)
from tissue_sim.generate_transcripts_with_coords import generate_transcripts
from tissue_sim.make_grid_csv import make_grid
from tissue_sim.place_cells import place_cells
from tissue_sim.provenance import (git_clean,git_commit,read_manifest,require_verified,stamp,)
from tissue_sim.reference import generate_reference
from tissue_sim.validate import validate_dataset



@dataclass
class SimulationResult:
    """In-memory tables and model parameters produced by one simulation."""

    coordinates: pd.DataFrame
    cells: pd.DataFrame
    transcripts: pd.DataFrame
    counts: pd.DataFrame
    composition: pd.DataFrame
    domains: pd.DataFrame
    params: dict
    reference_counts: pd.DataFrame
    reference_labels: pd.DataFrame


def simulate(config):
    """Return simulated tissue and independent reference tables without file I/O.

    config contains model settings in micrometres and captured transcripts per
    cell. One seeded generator drives all stages. Returned composition records
    true transcript counts per square/type, not cell counts or inferred fractions.
    """
    rng = np.random.default_rng(config.seed)
    grid = make_grid(
        config.field_width_um, config.field_height_um, config.square_size_um
    )
    coordinates = (
        pd.DataFrame.from_dict(grid, orient="index")
        .rename_axis("square_id")
        .reset_index()
    )
    _, placed = place_cells(coordinates.to_dict("records"), config, rng)
    if len(placed) != config.n_cells:
        raise ValueError(
            f"Only placed {len(placed)} of {config.n_cells} cells; reduce density or radius"
        )
    cells = assign_cell_types(
        pd.DataFrame(placed),
        config.n_domains,
        config.field_width_um,
        config.cell_types,
        config.domain_probabilities,
        rng,
    )
    domains = build_domains_dataframe(
        coordinates, config.n_domains, config.field_width_um
    )
    baseline = make_baseline_distribution(config.n_genes, config.baseline_spread, rng)
    markers = generate_marker_genes(
        config.cell_types, config.n_genes, config.n_marker_genes, rng
    )
    programs = generate_cell_type_programs(
        config.cell_types, baseline, markers, config.marker_boost
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        program_correlation = np.atleast_2d(
            np.corrcoef([programs[t] for t in config.cell_types])
        )

    program_correlation = np.where(
        np.isfinite(program_correlation),
        program_correlation,
        None,
    )

    genes = [f"GENE_{i + 1:04d}" for i in range(config.n_genes)]
    means = dict(zip(config.cell_types, config.mean_transcripts_per_type, strict=True))
    transcripts = generate_transcripts(
        cells,
        programs,
        genes,
        means,
        config.square_size_um,
        config.field_width_um,
        config.field_height_um,
        rng,
    )
    counts, composition = aggregate_transcripts(coordinates, transcripts, genes)
    reference_counts, reference_labels = generate_reference(
        programs, genes, means, config.reference_n_per_type, rng,
        sequencing_depth=config.sequencing_depth,
    )
    params = asdict(config)
    params.update(
        expression_programs=[programs[t].tolist() for t in config.cell_types],
        marker_genes={t: markers[t].tolist() for t in config.cell_types},
        n_markers_per_type=config.n_marker_genes,
        marker_gene_boost=config.marker_boost,
        baseline=baseline.tolist(),
        expression_program_correlation=program_correlation.tolist(),
    )
    return SimulationResult(
        coordinates,
        cells,
        transcripts,
        counts,
        composition,
        domains,
        params,
        reference_counts,
        reference_labels,
    )


def write_results(result, out_dir):
    """Export final tables, refusing to overwrite an existing output directory.

    Diagnostic cell/transcript tables are separate from aggregate ground truth.
    Coordinates are in micrometres; counts and composition contain captured
    transcripts. Export includes independent reference cells and contract metadata.
    Validate the written directory before returning; raises ValueError on failure.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    truth = out_dir / "ground_truth"
    truth.mkdir()
    reference = out_dir / "reference"
    reference.mkdir()
    diagnostics = out_dir / "diagnostics"
    diagnostics.mkdir()
    for table, path in (
        (result.coordinates, out_dir / "coordinates.csv"),
        (result.counts, out_dir / "counts.csv"),
        (result.composition, truth / "composition.csv"),
        (result.domains, truth / "domains.csv"),
        (result.reference_counts, reference / "reference_counts.csv"),
        (result.reference_labels, reference / "reference_labels.csv"),
        (result.cells, diagnostics / "cells.csv"),
        (result.transcripts, diagnostics / "transcripts.csv"),
    ):
        table.to_csv(path, index=False)
    (truth / "params.json").write_text(
        json.dumps(result.params, indent=2) + "\n", encoding="utf-8"
    )

    params = result.params
    clean = git_clean()
    commit = git_commit()
    metadata = {
        "contract_version": "1.0",
        "square_size_um": params["square_size_um"],
        "grid": {
            "n_rows": int(params["field_height_um"] / params["square_size_um"]),
            "n_cols": int(params["field_width_um"] / params["square_size_um"]),
        },
        "field_um": {
            "width": params["field_width_um"],
            "height": params["field_height_um"],
        },
        "n_squares": len(result.coordinates),
        "n_genes": params["n_genes"],
        "cell_types": list(params["cell_types"]),
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "generator": {
            "name": "tissue-sim",
            "version": version("tissue-sim"),
            "git_commit": (commit if clean else f"{commit}-dirty")
            if commit
            else "unavailable",
            "git_clean": clean,
        },
        "random_seed": params["seed"],
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    validate_dataset(out_dir)


def run(config_path, raw_dir, out_dir, overrides=None):
    """Verify supplied raw data, simulate in memory, export tables and provenance."""
    raw_dir = Path(raw_dir)

    manifest_path = raw_dir.parent / "raw_manifest.sha256"
    require_verified(raw_dir, manifest_path)
    config = read_simulation_config(config_path, overrides=overrides)

    if Path(out_dir).exists():
        raise FileExistsError(f"Choose a new output directory: {out_dir}")
    result = simulate(config)
    metadata = stamp(config_path, config.seed, manifest_path)
    metadata.update(
        command=sys.argv,
        platform=platform.platform(),
        packages={
            name: version(name) for name in ("tissue-sim", "numpy", "pandas")
        },
        config_sha256=hashlib.sha256(Path(config_path).read_bytes()).hexdigest(),
        config_resolved=asdict(config),
        raw_data_dir=str(raw_dir),
        raw_data_manifest=read_manifest(manifest_path),
        output_status="validated: Data Contract 1.0",
    )
    write_results(result, out_dir)
    (Path(out_dir) / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    return result


def main():
    """Read terminal arguments and run the connected simulation pipeline."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="CLI overrides take precedence over YAML values. "
               "Omitted overrides retain the corresponding YAML values.",
    )

    parser.add_argument(
        "--config", type=Path,
        default=Path("config/default.yaml"),
        help="Path to the YAML configuration file. All parameters must be explicitly defined."
    )

    parser.add_argument(
        "--data-dir", type=Path,
        default=Path("data/raw"),
        help="Directory containing the read-only course data."
    )

    parser.add_argument(
        "--out-dir", type=Path, required=True,
        help="New directory for simulation results."
    )

    parser.add_argument(
        "--seed", type=int,
        help="Override the random seed."
    )
    parser.add_argument(
        "--marker-boost", type=float,
        help="Override marker gene expression fold-change."
    )
    parser.add_argument(
        "--n-marker-genes", type=int,
        help="Override the number of marker genes per cell type."
    )
    parser.add_argument(
        "--sequencing-depth", type=float,
        help="Override reference sequencing depth multiplier."
    )
    parser.add_argument(
        "--reference-n-per-type", type=int,
        help="Override the number of reference cells per type."
    )

    args = parser.parse_args()

    overrides = {
        name: getattr(args, name)
        for name in (
            "seed",
            "marker_boost",
            "n_marker_genes",
            "sequencing_depth",
            "reference_n_per_type",
        )
        if getattr(args, name) is not None
    }

    result = run(
        args.config, args.data_dir, args.out_dir,
        overrides=overrides,
    )

    print(
        f"Wrote {len(result.coordinates)} squares and {len(result.transcripts)} transcripts to {args.out_dir}"
    )


if __name__ == "__main__":
    main()
