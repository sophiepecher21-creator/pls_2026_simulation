
def make_baseline_distribution(number_of_genes, baseline_parameter, random_generator):
    """
    Make the baseline distribution for one cell

    Args:
        number_of_genes (int): Number of genes
        baseline_parameter (float): Baseline parameter
        random_generator (RandomGenerator): Random generator

    Returns:
        (list): a probability for each gene, summing to 1
    """


    # test: correct length, finite non-negative values, probabilities sum to 1
    return []

def make_expression_program(baseline_distribution, marker_genes_indices, marker_boost):
    """
    Make the expression program for one cell

    Args:
        baseline_distribution (list): baseline distribution
        marker_genes_indices (list): marker genes indices
        marker_boost (int): marker gene boost

    Returns:
        (list): a normalized expression program with boosted marker probabilities.
    """

    # test: markers are boosted, non-marker genes retain their baseline proportions, probabilities sum to 1
    return []

def make_expression_programs(cell_types, expression_parameters, random_generator):
    """
    Make the expression program for each cell type

    Args:
        cell_types (list): cell types
        expression_parameters (list): expression parameters
        random_generator (RandomGenerator): Random generator

    Returns:
        (dict): a dictionary that maps each cell type to its expression program
    """

    # test: all cell types have programs, marker sets are valid, and the same seed reproduces the same results
    return {}