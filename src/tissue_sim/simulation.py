"""Connect the simulation stages in memory and export their final tables."""

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tissue_sim.aggregation import aggregate_transcripts
from tissue_sim.configuration import read_simulation_config
from tissue_sim.domains import assign_cell_types, build_domains_dataframe
from tissue_sim.expression_programs import (
    generate_cell_type_programs,
    generate_marker_genes,
    make_baseline_distribution,
)
from tissue_sim.generate_transcripts_with_coords import generate_transcripts
from tissue_sim.make_grid_csv import make_grid
from tissue_sim.place_cells import place_cells
from tissue_sim.provenance import require_verified, stamp


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


def simulate(config):
    """Run scientific stages with one shared seeded Generator and no file I/O."""
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
    params = asdict(config)
    params.update(
        expression_programs=[programs[t].tolist() for t in config.cell_types],
        marker_genes={t: markers[t].tolist() for t in config.cell_types},
        n_markers_per_type=config.n_marker_genes,
        marker_gene_boost=config.marker_boost,
        baseline=baseline.tolist(),
    )
    return SimulationResult(
        coordinates, cells, transcripts, counts, composition, domains, params
    )


def write_results(result, out_dir):
    """Export final tables, refusing to overwrite an existing output directory.

    Diagnostic cell/transcript tables are separate from aggregate ground truth.
    This milestone does not yet produce reference/ or a full contract dataset.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    truth = out_dir / "ground_truth"
    truth.mkdir()
    diagnostics = out_dir / "diagnostics"
    diagnostics.mkdir()
    for table, path in (
        (result.coordinates, out_dir / "coordinates.csv"),
        (result.counts, out_dir / "counts.csv"),
        (result.composition, truth / "composition.csv"),
        (result.domains, truth / "domains.csv"),
        (result.cells, diagnostics / "cells.csv"),
        (result.transcripts, diagnostics / "transcripts.csv"),
    ):
        table.to_csv(path, index=False)
    (truth / "params.json").write_text(
        json.dumps(result.params, indent=2) + "\n", encoding="utf-8"
    )


def run(config_path, raw_dir, out_dir):
    """Verify supplied raw data, simulate in memory, export tables and provenance."""
    raw_dir = Path(raw_dir)
    manifest_path = raw_dir.parent / "raw_manifest.sha256"
    require_verified(raw_dir, manifest_path)
    config = read_simulation_config(config_path)
    if Path(out_dir).exists():
        raise FileExistsError(f"Choose a new output directory: {out_dir}")
    result = simulate(config)
    metadata = stamp(config_path, config.seed, manifest_path)
    metadata["output_status"] = (
        "partial: reference experiment and contract validator pending"
    )
    write_results(result, out_dir)
    (Path(out_dir) / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main():
    """Read terminal arguments and run the connected simulation pipeline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.config, args.data_dir, args.out_dir)
    print(
        f"Wrote {len(result.coordinates)} squares and {len(result.transcripts)} transcripts to {args.out_dir}"
    )


if __name__ == "__main__":
    main()
