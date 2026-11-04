"""Simulate a separate single-cell experiment from the model's type programs."""

import pandas as pd


def generate_reference(programs, genes, means, n_per_type, rng):
    """Return independently sampled reference count and label DataFrames.

    programs maps types to true gene probabilities, means to captured transcripts
    per whole cell. Draw a fresh Poisson total and multinomial gene counts per
    reference cell, using rng; never copy sampled tissue cells or transcripts.
    Rows share cell_id order; gene order follows genes. n_per_type controls the
    number of profiled cells per type (80 by default in SimulationConfig).
    """
    rows, labels = [], []
    for cell_type, probabilities in programs.items():
        for _ in range(n_per_type):
            cell_id = f"ref_{len(rows) + 1:04d}"
            counts = rng.multinomial(rng.poisson(means[cell_type]), probabilities)
            rows.append([cell_id, *counts])
            labels.append([cell_id, cell_type])
    return (
        pd.DataFrame(rows, columns=["cell_id", *genes]),
        pd.DataFrame(labels, columns=["cell_id", "cell_type"]),
    )
