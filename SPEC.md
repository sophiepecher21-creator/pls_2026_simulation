# <project name> — specification

## What it does

`tissue-sim` simulates a 2D tissue of cells that are organized in spatial domains and express genes stochastically depending on their cell type. The program then mimics a Visium HD experiment by assigning the generated transcripts to a grid of 2 µm squares. In addition to the simulation, the program provides ground truth about the underlying tissue and independent reference data, so that spatial transcriptomics analysis methods can be tested and evaluated.

## Inputs

<What it reads. Refer to the Data Contract rather than inventing a format.>

## Outputs

<What it writes. Again, per the Data Contract.>

## Acceptance criteria

Concrete, checkable statements of "how we will know it is right". Tag each with the kind of
check that enforces it (smoke, known-answer, property, metamorphic, characterization,
schema/validation, reproducibility).

1. `[check type]` <criterion>
2. `[check type]` <criterion>
3. `[check type]` <criterion>

## First known answer

Input: <the smallest input whose correct output you can state by hand>
Expected output: <what a correct program must produce for it>
