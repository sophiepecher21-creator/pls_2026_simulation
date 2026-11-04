# PLS_2026_Simulation — specifications

## What it does

`tissue-sim` simulates a 2D tissue of cells that are organized in spatial domains and express genes stochastically depending on their cell type. The program then mimics a Visium HD experiment by assigning the generated transcripts to a grid of 2 µm squares. In addition to the simulation, the program provides ground truth about the underlying tissue and independent reference data, so that spatial transcriptomics analysis methods can be tested and evaluated.

## Inputs

seed : Predefined randomness seeding 
Field size µm: 0-N gid size
N: Number of Cells
Celltypes: 0-N celltypes 
Genes: 0-N Genes defined for each cell type
Marker genes: Defined marker gene for a cell type, can be multiple
Marker gene bost: Makes expression of marker genes stronger/weaker

The Week 3 toy reads config/default.yaml for its seed, cell count, field
dimensions, and two cell-type names. Before generating output, it verifies
all files under data/raw/ against data/raw_manifest.sha256. The supplied
reference dataset is retained as an external resource; the toy does not
read its ground truth.

## Outputs

Directory layout produced by the simulator:

```text
<dataset_name>/ 
├── metadata.json    # Experiment specifications
├── coordinates.csv   # Position of each square on the grid
├── counts.csv    # Number of transcripts detected for each gene in each spatial square (squares x genes)
├── reference/   # Reference for scRNA-seq
│  ├── reference_counts.csv    # Counts for the simulated reference cells
│  └── reference_labels.csv   # True cell type for each reference cell
└── ground_truth/   # "Answer key" - analyzes must NOT read this
   ├── params.json   # Parameters used to generate the simulation
   ├── composition.csv   # Transcript counts per square and cell type
   └── domains.csv   # True spatial domain label for each square
```

## Acceptance criteria

The table maps the current Weeks 6–7 acceptance criteria to existing checks.
Coverage gaps are noted explicitly.

| Claim | Kind of check | Test file and function |
| --- | --- | --- |
| 1. At least two function docstrings cover grid geometry and transcript-to-square assignment. | Known-answer / doctest | `src/tissue_sim/make_grid_csv.py::make_grid`; `src/tissue_sim/generate_transcripts_with_coords.py::get_square_id`. Run with `pytest --doctest-modules`. |
| 2. README file paths, local links and code-block file arguments exist. | Schema / validation | In `tests/test_documentation.py`: `test_readme_paths_exist`, `test_path_checker_known_answer`, `test_checker_catches_command_and_link_paths`. |
| 3. Dataset validation enforces Contract 1.0, valid tables, parameter consistency and transcript conservation; absent composition pairs mean zero. | Schema / validation | In `tests/test_validate.py`: `test_generated_dataset_conforms`, `test_empty_composition_is_zero_not_missing`, `test_invalid_json_values`, `test_invalid_tables`, `test_missing_and_malformed_json`, `test_notebook_layer_checks`. These exercise conformance and selected malformed inputs. |
| 4. Simulation exports all required files, validates output and independently samples configurable reference cells. | Integration / known-answer | In `tests/test_smoke.py`: `test_run_exports_and_refuses_overwrite`. In `tests/test_validate.py`: `test_generated_dataset_conforms`, `test_reference_known_answer_and_fresh_sampling`. Direct checks of automatic validation before success, independence from tissue and configurable reference sample size remain to be confirmed. |
| 5. Six separate tests reject dataset_A–E for their specific faults and silently accept dataset_OK. External tests skip if the pack is absent; generated-fixture checks run independently. | Regression | In `tests/test_validate_pack.py`: `test_dataset_A`, `test_dataset_B`, `test_dataset_C`, `test_dataset_D`, `test_dataset_E`, `test_dataset_OK`. In `tests/test_validate.py`: `test_generated_dataset_conforms`, `test_invalid_json_values`, `test_invalid_tables`, `test_missing_and_malformed_json`. |
| 6. Fixed configuration and seed reproduce numerical exports, including reference data; metadata records dirty Git provenance. | Reproducibility / provenance | In `tests/test_smoke.py`: `test_pipeline_reproducible_and_conserves_counts`, `test_run_exports_and_refuses_overwrite`. A direct check of dirty Git metadata remains to be confirmed. |
| 7. Supplied inputs retain their hashes, both packs are covered by the manifest and verification precedes output. | Provenance | In `tests/test_provenance.py`: `test_untouched_directory`, `test_modified_file`, `test_missing_and_added_files`, `test_changed_input_stops_processing`. In `tests/test_smoke.py`: `test_changed_raw_data_prevents_output`. These use temporary inputs; actual supplied-pack hashes and manifest coverage require a separate check. |


## Week 3 acceptance criteria

1. [smoke test] The current toy runs end to end and writes toy_cells.csv
   and run_metadata.json to the requested output directory.
2. [reproducibility test] A teammate can create the environment from
   environment.yml, install the package with pip install -e ., and run
   pytest successfully without manually installing additional dependencies.
   Tests use temporary inputs rather than the external reference dataset.
3. [schema / validation check] The program refuses to run if a raw file
   is changed, missing, or added relative to the recorded manifest.
4. [reproducibility test] Deleting results and rerunning with the same
   configuration and seed reproduces toy_cells.csv byte for byte.
   The provenance timestamp records each run and may differ.

## One file we cannot regenerate

The course-supplied
data/raw/reference_dataset/dataset/ground_truth/composition.csv
is given data. Our project cannot recreate that exact supplied file from
its own code and configuration. We preserve it read-only and record its
checksum. Ground truth generated by our own simulator is a separate,
generated artifact. 

## First known answer

Input: When we input information to create a square of cells at 2-D positions, with spatial domains, cell types that are dominant in certain domains, and some genes/marker genes given the spatial constraints of a Visium grid...

Expected output: The output should, at minimum, be able to bin output to 8 µm and 16 µm, four 8 µm bins sum exactly to their 16 µm bin, gene by gene. Every (row, col) appears exactly once; µm coordinates match the indices. 


## Week 6 connected pipeline

The Week 3 toy is superseded by the connected pipeline. The original criteria
above describe that historical milestone. Current acceptance criteria:

1. A single simulation command connects grid, placement, domains, expression,
   transcript generation and aggregation through in-memory objects.
2. Repeated runs with the same configuration and seed reproduce numerical tables.
3. Gene counts and square/type composition both conserve in-field transcripts;
   empty squares remain in the grid and gene matrix.
4. A known-answer aggregation combines contributions across cells, and exact
   8-to-16 µm aggregation holds gene by gene.
5. Unsafe YAML constructors and malformed settings are rejected before sampling.
6. Raw-data verification still precedes output; existing result directories are protected.

Independent reference sampling and full contract metadata/validation remain
pending; the current output is explicitly a partial dataset. See README for its
layout. `diagnostics/` retains cell/transcript details outside aggregate ground truth.

## Weeks 6 and 7 acceptance criteria (supersede partial-output status above)

1. [known-answer test] `pytest --doctest-modules` runs at least two function
   docstrings, covering grid geometry and transcript-to-square assignment.
2. [schema / validation check] README inline file paths, local Markdown links
   and file arguments in code blocks point to files that exist.
3. [schema / validation check] `validate_dataset(path)` checks Contract 1.0
   metadata, required files/columns, full unique grid coverage, coordinate units,
   IDs and ordering, non-negative integer gene counts, gene alignment, labeled
   reference coverage, ground truth, parameter consistency and exact per-square
   transcript conservation. Missing sparse composition pairs mean zero.
4. [integration test] The simulator writes every required contract file and calls
   validation on that directory before returning success. Independent reference
   cells are sampled from true programs with fresh random draws, never copied
   from the tissue. Reference sample size is configurable (default 80 per type).
5. [regression test] Six separate tests reject dataset_A through dataset_E for
   their specific faults and accept dataset_OK silently. External fixtures remain
   read-only and untracked; these six tests skip explicitly if the pack is absent.
   Generated-fixture tests always exercise conformance and malformed inputs.
6. [reproducibility test] Fixed seed/configuration reproduces numerical exports,
   including the independent reference experiment. UTC timestamps may differ.
   Metadata explicitly marks dirty Git provenance rather than implying clean code.
7. [provenance check] The original reference files retain their recorded hashes.
   The manifest additionally covers the new broken dataset pack. No supplied
   input is altered; simulation verifies both packs before writing output.
