# PLS_2026_Simulation — specifications

## What it does

<What does the program do >
  
## Inputs

<What it reads. Refer to the Data Contract rather than inventing a format.>

## Outputs

<What it writes. Again, per the Data Contract.>

## Acceptance criteria

Concrete, checkable statements of "how we will know it is right". Tag each with the kind of
check that enforces it (smoke, known-answer, property, metamorphic, characterization,
schema/validation, reproducibility).

1. Does the simulator produce a coordinates .csv file which agrees with the technological specs. 2 um squares, 160x160 with each square id being unique and correlating correctly with the grid layout. Also, x & y need to have correct coords (x0 = 1 um, y0 = 1 um, and x160 = 320 um, y160 = 320 um)
   - Known answer test: Binning output to 8 µm and 16 µm, four 8 µm bins sum exactly to their 16 µm bin, gene by gene
   - Property test: Every (row, col) appears exactly once; µm coordinates match the indices
2. Metamorphic check: If we change the formatting of gene IDs does everything else stay the same? Could someone using a different reference genome use this tool?
   - Known answer test: For each of a cell type's marker genes, that gene's proportion is higher in that type's program than in any other type's regardless of reference genome
   - Property test: No gene is a marker for two different types
3. Property Check/Reproducibility Check: The simulation input should be the same every time, except the counts themselves can change. The counts should only change such that the sum of counts for a barcode are within the specs of what 10X visium claims the readout should be.
   - Property check: The total counts over all squares equal the number of simulated transcripts inside the field
   - Property check: Ground truth matches the counts: each square's per-type transcript counts sum to its total in `counts.csv`
   - Reproducibility test: Two runs, same seed ⇒ identical files
   - Property test: every `p_t`, and the baseline, sums to 1 over genes
4. Spatial domains are functioning as intended: Partition the field into `K` number of vertical strips of equal width,
   numbered 0 to `K-1` from left to right. Strips are contiguous by construction.
   - Property test: domains are contiguous, each domain forms a single connected region.
5. The model learns how to predict cell types correctly when supplemented with training data
   - Reference is independent, property test: no cell used in the tissue appears in `reference/`
   - Metamorphic test: Signatures estimated from `reference/` approach the true `p_t` as reference size grows
   - Metamorphic test: Lowering `marker_gene_boost` raises the correlation between two types' programs, and a boost of 1 makes every program equal to the baseline
   - Separation -> recoverability: A reference analyzer recovers domains at ARI > 0.8 on easy settings (it will get worse as the settings get harder) 

## First known answer

Input: When we input information to create a square of cells at 2-D positions, with spatial domains, cell types that are dominant in certain domains, and some genes/marker genes given the spatial constraints of a Visium grid...

Expected output: The output should, at minimum, be able to bin output to 8 µm and 16 µm, four 8 µm bins sum exactly to their 16 µm bin, gene by gene. Every (row, col) appears exactly once; µm coordinates match the indices. 
