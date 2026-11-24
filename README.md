# Programming for Life Sciences 2026: tissue-sim

**Group:** Sophie, Sven, Ronja, Martina · **Track:** simulation

tissue-sim generates a two-dimensional tissue, measures its transcripts on a
2 µm grid, and produces independent single-cell reference data and aggregate
ground truth for testing spatial transcriptomics methods. Stages pass Python
objects in memory. Every exported dataset is checked with `validate_dataset(path)`
before the run succeeds. File formats are defined in `docs/data_contract.pdf`.

## Install

### Get the project files
Before installing, download or clone this repository onto your computer.

- **Download ZIP:** On the repository page, select **Code → Download ZIP**, then extract the downloaded archive.
- **Git:** Clone the repository using its repository URL.

Open a terminal and navigate into the extracted or cloned project folder:

```bash
cd path/to/pls_2026_simulation
```

Replace `path/to/pls_2026_simulation` with the location of the folder on your computer. If you downloaded a ZIP, the extracted folder may have a suffix such as `-main`.

Check that this folder contains `environment.yml` and `pyproject.toml`.

### Install the environment and project

Run each command separately. If a command fails, resolve the error before continuing.

```bash
conda env create -f environment.yml
conda activate pls2026
python -m pip install -e .
```

After activation, your terminal prompt should show `(pls2026)`. Run the remaining commands in this README from the project folder unless instructed otherwise.

The environment recipe is `environment.yml`; model settings are in
`config/default.yaml`. Python 3.11 or newer is required.

## Prepare course data

Obtain the course reference dataset and broken dataset pack from ADAM. Copy the
reference_dataset directory and broken_datasets directory into data/raw, preserving
their contents. You need to create the raw folder. These external inputs are read-only and excluded from Git.
Do not include Finder metadata or notebook checkpoints. The committed
`data/raw_manifest.sha256` fingerprints both packs; the simulation refuses changed,
missing or additional raw files. Tests of the six course datasets skip explicitly
if the pack has not been installed; generated-fixture validator tests always run.

Your data path should now look something like './data/raw/' and contain both the reference dataset and the broken datasets.

## Run

Display all available command-line options:

```bash
python run.py --help
```

From the repository root, choose an unused output directory:

```bash
python run.py --out-dir results/run_001
```
Replace `run_001` with any directory name of your choice. The output directory is created automatically and must not already exist.

By default, the program uses `config/default.yaml` and `data/raw`.
These paths can be changed using `--config` and `--data-dir`.

Selected scientific parameters can be overridden directly:

```bash
python run.py --out-dir results/run_002 --seed 42 --marker-boost 5.0
```

Configuration precedence (lowest to highest priority):
Python defaults → YAML configuration → CLI overrides.

The configuration must explicitly define `cell_types`, `domain_probabilities` and `mean_transcripts_per_type`. Unknown keys and invalid parameter values are rejected.

The original module-based command remains supported.

The simulation writes contract metadata, square coordinates, a square-by-gene transcript count matrix,
independently sampled reference counts and labels, and ground truth containing
parameters, sparse square/type transcript counts and square domain labels.
Cell and transcript diagnostic tables are kept separately from aggregate truth.
See `docs/data_contract.pdf` for the exact directory layout and schemas.

The supplied raw files are verified for provenance; they are not used to generate
the tissue. Reference cells use the same true expression programs as the tissue
but are sampled afresh, with 80 cells per type by default. Their number is controlled
by `reference_n_per_type`. The `sequencing_depth` multiplier affects reference transcript recovery only.
Expression-program correlations are stored in ground_truth/params.json,
with undefined correlations represented as null.

The resolved configuration and command-line invocation are recorded
in run_metadata.json. Fixed configuration, seed, software environment and stage
order reproduce the numerical tables; metadata timestamps vary. Metadata marks a
modified working tree explicitly, so a dirty run is not attributed to a clean commit.

Validate an existing output independently (success is silent):

```bash
python -m tissue_sim.validate results/run_001
```

Replace `results/run_001` with the output directory you want to validate.

For Python use, import `validate_dataset` from `tissue_sim.validate`. It raises
`ValueError` naming the file, offending value and failed requirement.

## Check

```bash
pytest --doctest-modules
ruff check .
```

The first command runs ordinary tests and executable docstring examples, including
the README path check and one test per supplied broken dataset. The path check
covers inline file paths, local Markdown links and file arguments in code blocks;
it cannot establish that a new reader can follow the instructions.

Acceptance criteria are in `SPEC.md`; the AI and verification records are in
`ai_log/`.
