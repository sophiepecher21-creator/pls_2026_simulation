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

```bash
conda env create -f environment.yml
conda activate pls2026
pip install -e .
```

The environment recipe is `environment.yml`; model settings are in
`config/default.yaml`. Python 3.11 or newer is required.

## Prepare course data

Obtain the course reference dataset and broken dataset pack from ADAM. Copy the
reference_dataset directory and broken_datasets directory into data/raw, preserving
their contents. These external inputs are read-only and excluded from Git.
Do not include Finder metadata or notebook checkpoints. The committed
`data/raw_manifest.sha256` fingerprints both packs; the simulation refuses changed,
missing or additional raw files. Tests of the six course datasets skip explicitly
if the pack has not been installed; generated-fixture validator tests always run.

## Run

From the repository root, choose an unused output directory:

```bash
python -m tissue_sim.simulation \
  --config config/default.yaml \
  --data-dir data/raw \
  --out-dir results/run_001
```

The command creates its output directory and refuses to overwrite it. It writes
contract metadata, square coordinates, a square-by-gene transcript count matrix,
independently sampled reference counts and labels, and ground truth containing
parameters, sparse square/type transcript counts and square domain labels.
Cell and transcript diagnostic tables are kept separately from aggregate truth.
See `docs/data_contract.pdf` for the exact directory layout and schemas.

The supplied raw files are verified for provenance; they are not used to generate
the tissue. Reference cells use the same true expression programs as the tissue
but are sampled afresh, with 80 cells per type by default. Their number is controlled
by `reference_n_per_type`. Fixed configuration, seed, software environment and stage
order reproduce the numerical tables; metadata timestamps vary. Metadata marks a
modified working tree explicitly, so a dirty run is not attributed to a clean commit.

Validate an existing output independently (success is silent):

```bash
python -m tissue_sim.validate results/run_001
```

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

Acceptance criteria are in `SPEC.md`; the AI and verification record is in
`ai_log/sophiepecher.md`.
