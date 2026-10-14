def draw_transcript_count(mean_transcripts, random_generator):
    """ Determine how many transcripts are captured from one cell """
    pass

def sample_transcript_genes(number_of_transcripts, gene_program, random_generator):
    """Assign a gene to each transcript according to the cell type's expression program"""
    pass

def sample_positions_in_cell(number_of_transcripts, center_x, center_y, radius, random_generator):
    """Generate a random spatial position for each transcript inside a cell's circular area"""
    pass

def simulate_cell_transcripts(cell, expression_programs, mean_transcripts_by_type, random_generator):
    """Generate all transcripts belonging to one cell."""
    pass

def simulate_tissue_transcripts(cells, expression_programs, mean_transcripts_by_type, seed):
    """Generate transcripts for all cells in the tissue using one reproducible seed."""
    pass