import numpy as np
from pathlib import Path
import yaml

def make_baseline_distribution(n_genes, baseline_spread, random_generator):
    """
    Generate a long-tailed baseline probability distribution over genes using a log-normal distribution.

    Args:
        n_genes (int): Number of genes in the distribution.
        baseline_spread (float): Standard deviation (sigma) of the baseline distribution used to generate log-normal distribution.
        random_generator (numpy.random.Generator): Random number generator used to generate log-normal distribution.

    Returns:
        (numpy.ndarray): One-dimensional baseline probability distribution over genes.
    """

    # Validate n_genes
    if (
            isinstance(n_genes, bool)
            or not isinstance(n_genes, (int, np.integer))
            or n_genes <= 0
    ):
        raise ValueError("n_genes must be a positive integer.")

    # Validate baseline_spread
    if (
        isinstance(baseline_spread, bool)
        or not isinstance(baseline_spread, (int, float))
        or not np.isfinite(baseline_spread)
        or baseline_spread < 0
    ):
        raise ValueError("baseline_spread must be a finite, non-negative number.")

    # Validate random_generator
    if not isinstance(random_generator, np.random.Generator):
        raise ValueError("random_generator must be a NumPy Generator.")

    weights = random_generator.lognormal(
        mean=0.0,
        sigma=baseline_spread,
        size=n_genes,
    )

    probabilities = weights / weights.sum()
    return probabilities


def boost_marker_genes(baseline, marker_genes_indices, marker_boost):
    """
    Increase probability of marker genes and renormalize the probabilities.

    Args:
        baseline (numpy.ndarray): One-dimensional baseline probability distribution over genes.
        marker_genes_indices (list): Indices of marker genes to boost.
        marker_boost (int or float): Multiplicative factor applied to marker genes, must be at least 1.

    Returns:
        (numpy.ndarray): Normalized expression probability distribution with boosted marker
        genes.

    """

    # Validate marker indices before converting them to integers
    raw_indices = np.asarray(marker_genes_indices)

    if raw_indices.ndim != 1:
        raise ValueError("Marker indices must be a 1D sequence of integers.")

    if not np.issubdtype(raw_indices.dtype, np.integer):
        raise ValueError("Marker indices must be integers.")

    marker_indices = raw_indices.astype(int)

    # Validate baseline
    baseline = np.asarray(baseline, dtype=float)

    if baseline.ndim != 1 or len(baseline) == 0:
        raise ValueError("Baseline distribution must be a non-empty 1D array.")

    if (
            not np.all(np.isfinite(baseline))
            or np.any(baseline < 0)
            or not np.isclose(baseline.sum(), 1.0)
    ):
        raise ValueError("Baseline must contain valid probabilities summing to 1.")

    if (
            np.any(marker_indices < 0)
            or np.any(marker_indices >= len(baseline))
            or len(np.unique(marker_indices)) != len(marker_indices)
    ):
        raise ValueError("Marker indices must be unique, valid gene indices.")

    # Validate marker_boost
    if (
            isinstance(marker_boost, bool)
            or not isinstance(marker_boost, (int, float))
            or not np.isfinite(marker_boost)
            or marker_boost < 1
    ):
        raise ValueError("marker_boost must be finite and at least 1.")

    boosted = baseline.copy()
    boosted[marker_indices] *= marker_boost

    return boosted / boosted.sum()


def generate_marker_genes(cell_types, n_genes, n_marker_genes, random_generator):
    """
    Select non-overlapping marker gene indices for each cell type.

    Args:
        cell_types (sequence of str): Names of the cell types.
        n_genes (int): Total number of genes.
        n_marker_genes (int): Number of marker genes to select for each cell type.
        random_generator (numpy.random.Generator): Random number generator used for marker gene selection.

    Returns:
         (dict): Dictionary mapping each cell type to a NumPy array of marker
        gene indices.
    """

    # Validate cell_types
    if len(cell_types) == 0:
        raise ValueError("cell_types must contain at least one cell type.")

    if len(set(cell_types)) != len(cell_types):
        raise ValueError("cell_types must contain unique names.")

    # Validate n_genes
    if (
        isinstance(n_genes, bool)
        or not isinstance(n_genes, (int, np.integer))
        or n_genes <= 0
    ):
        raise ValueError("n_genes must be a positive integer.")

    # Validate n_marker_genes
    if (
        isinstance(n_marker_genes, bool)
        or not isinstance(n_marker_genes, (int, np.integer))
        or not 0 <= n_marker_genes <= n_genes
    ):
        raise ValueError(
            "n_marker_genes must be an integer between 0 and n_genes."
        )

    # Validate random_generator
    if not isinstance(random_generator, np.random.Generator):
        raise ValueError("random_generator must be a NumPy Generator.")


    if len(cell_types) * n_marker_genes > n_genes:
        raise ValueError(
            "Not enough genes for non-overlapping marker sets."
        )

    available_genes = list(range(n_genes))
    marker_genes = {}

    for cell_type in cell_types:
        selected = random_generator.choice(
            available_genes,
            size=n_marker_genes,
            replace=False,
        )
        marker_genes[cell_type] = selected

        # Selected genes do not overlap between cell types
        selected_set = set(selected)
        available_genes = [gene for gene in available_genes if gene not in selected_set]

    return marker_genes


def generate_cell_type_programs(cell_types, baseline, marker_genes, marker_boost):
    """
    Generate a normalized expression program for each cell type.

    Args:
        cell_types (sequence of str): Names of the cell types.
        baseline (numpy.ndarray): One-dimensional baseline probability distribution over genes.
        marker_genes (dict): Dictionary mapping each cell type to its marker gene indices.
        marker_boost (int or float): Multiplicative factor applied to marker genes.

    Returns:
        (dict): Dictionary mapping each cell type to a NumPy array containing
        its normalized expression probabilities.
    """

    # Validate cell_types
    if len(cell_types) == 0:
        raise ValueError("cell_types must contain at least one cell type.")

    if len(set(cell_types)) != len(cell_types):
        raise ValueError("cell_types must contain unique names.")

    # Validate baseline_distribution
    baseline = np.asarray(baseline, dtype=float)

    if baseline.ndim != 1 or baseline.size == 0:
        raise ValueError(
            "baseline_distribution must be a non-empty 1D array."
        )

    # Validate marker_genes
    if not isinstance(marker_genes, dict):
        raise ValueError("marker_genes must be a dictionary.")

    if set(marker_genes) != set(cell_types):
        raise ValueError(
            "marker_genes must contain exactly one entry for each cell type."
        )

    expression_programs = {}

    for cell_type in cell_types:
        expression_programs[cell_type] = boost_marker_genes(
            baseline=baseline,
            marker_genes_indices=marker_genes[cell_type],
            marker_boost=marker_boost,
        )

    return expression_programs



def main():
    """Build expression programs for the default configuration as a standalone example."""
    # Load configuration
    config_path = (
            Path(__file__).resolve().parents[2]
            / "config"
            / "default.yaml"
    )
    with config_path.open() as file:
        config = yaml.safe_load(file)

    # Get parameters from config
    cell_types = config["cell_types"]
    n_genes = config["n_genes"]
    n_marker_genes = config["n_marker_genes"]
    baseline_spread = config["baseline_spread"]
    marker_boost = config["marker_boost"]

    # Reproducible random number generator
    rng = np.random.default_rng(config["seed"])

    # 1. Generate the shared baseline distribution
    baseline = make_baseline_distribution(
        n_genes=n_genes,
        baseline_spread=baseline_spread,
        random_generator=rng,
    )

    # 2. Select non-overlapping marker genes
    marker_genes = generate_marker_genes(
        cell_types=cell_types,
        n_genes=n_genes,
        n_marker_genes=n_marker_genes,
        random_generator=rng,
    )

    # 3. Generate an expression program for each cell type
    expression_programs = generate_cell_type_programs(
        cell_types=cell_types,
        baseline=baseline,
        marker_genes=marker_genes,
        marker_boost=marker_boost,
    )

    # Convert outputs to ordinary Python lists
    expression_programs = {
        cell_type: probabilities.tolist()
        for cell_type, probabilities in expression_programs.items()
    }
    #print(expression_programs)


if __name__ == "__main__":
    main()