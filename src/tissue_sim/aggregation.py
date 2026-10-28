"""Aggregate in-memory transcripts to instrument-resolution tables."""

import numpy as np
import pandas as pd


def aggregate_transcripts(coordinates, transcripts, gene_names):
    """Return dense square-by-gene counts and sparse square/type composition.

    Preserve coordinate and gene order, including all empty squares and genes.
    Composition has one positive count per square/type pair and no cell IDs.
    """
    ids = coordinates.square_id.tolist()
    if not transcripts.square_id.isin(ids).all():
        raise ValueError("Transcript references an unknown square")
    if not transcripts.gene.isin(gene_names).all():
        raise ValueError("Transcript references an unknown gene")
    matrix = np.zeros((len(ids), len(gene_names)), dtype=np.int64)
    rows = transcripts.square_id.map({v: i for i, v in enumerate(ids)}).to_numpy(
        dtype=int
    )
    cols = transcripts.gene.map({v: i for i, v in enumerate(gene_names)}).to_numpy(
        dtype=int
    )
    np.add.at(matrix, (rows, cols), 1)
    counts = pd.DataFrame(matrix, columns=gene_names)
    counts.insert(0, "square_id", ids)
    composition = (
        transcripts.groupby(["square_id", "cell_type"], sort=False)
        .size()
        .rename("n_transcripts")
        .reset_index()
    )
    order = {v: i for i, v in enumerate(ids)}
    composition = composition.sort_values(
        "square_id", key=lambda column: column.map(order), kind="stable"
    ).reset_index(drop=True)
    return counts, composition
