# Programming for Life Sciences 2026, Simulation Path

**Group:** Sophie, Sven, Ronja, Martina
**Track:** simulation (`tissue-sim`)

> Replace everything in angle brackets. This README is graded, and the criterion is simple:
> can someone who has never met you clone this repository and reproduce a result?
>
> **One repository per group.** One of you creates it from the course template with git's *Use this
> template*; the other two are added under *Settings → Collaborators* and clone it. Everyone making their own copy is the tangle we would spend a session undoing.

## What this does

tissue-sim simulates a two-dimensional tissue for testing spatial
transcriptomics analysis methods. The Week 6 pipeline builds the measurement
grid, places non-overlapping cells, assigns domains and types, builds expression
programs, generates positioned transcripts, and aggregates them to final tables.
Stages exchange Python objects; no intermediate CSV is read by the pipeline.

This is a partial dataset: the independent reference experiment, contract
metadata, and full dataset validator are planned for the next milestone.

## Install

```bash
conda env create -f environment.yml
conda activate pls2026
pip install -e .
```

The `-e` points the install at the source files, so an edit takes effect without reinstalling.
Someone who only wants to *run* this project rather than work on it can leave it out; see the
course document *Running a project from GitHub*.

## Run

```bash
python -m tissue_sim.simulation \
  --config config/default.yaml \
  --data-dir data/raw \
  --out-dir results/connected_run_001
```

## Repository layout

```
src/        the package you write
tests/      the checks that prove it works
config/     parameters — no magic numbers in the source
SPEC.md     what this program must do, and how you will know it is correct
ai_log/     your AI + verification log — one file per person, named for you
```

## Data

The reference dataset is **not** in this repository. Refer to it by path; never commit data.

The course-provided reference_dataset directory must be obtained separately
and copied to data/raw/reference_dataset/. It is not committed to Git.
Exclude .DS_Store and .ipynb_checkpoints when copying it, then make data/raw/
read-only with chmod -R a-w data/raw.

The committed data/raw_manifest.sha256 records the supplied files'
fingerprints. Every run verifies data/raw/ against this manifest before
creating output. The Week 3 toy verifies these supplied files but does not
use their contents to generate cells.

Supplied ground truth remains under
data/raw/reference_dataset/dataset/ground_truth/. It is given data that
this project cannot regenerate. Ground truth produced by our own simulator
will instead be generated output.

Choose a new output directory for each run; existing directories are protected.
The connected run writes:

```text
results/connected_run_001/
├── coordinates.csv           # complete square grid, five contract columns
├── counts.csv                # all squares × all genes, including zeros
├── ground_truth/
│   ├── composition.csv       # sparse square/type transcript counts
│   ├── domains.csv           # square-centre domain labels
│   └── params.json           # settings, baseline, programs, marker indices
├── diagnostics/
│   ├── cells.csv             # all placed cells with geometry, type and domain
│   └── transcripts.csv       # individual in-field transcripts
└── run_metadata.json         # provenance and partial-dataset status
```

`config/default.yaml` is loaded safely once into `SimulationConfig`. Use YAML
numbers and lists, not quoted numbers or Python expressions. Domain probability
rows correspond to domains; columns and transcript means follow `cell_types` order.
Transcript means describe captured RNA per type and can differ between types.
The new defaults use two mixed domains and means of 20 and 40 transcripts.

One shared NumPy generator is seeded once per run. Fixed configuration, seed,
software environment and stage order reproduce the numerical tables exactly.
This establishes a new baseline compared with the former independent scripts,
which restarted their generators. The provenance timestamp may differ.
Transcripts outside the field are discarded; all placed cells contribute even
if their disk contains no measurement-square centre. Composition combines
contributions from different cells of the same type. Crowded placements that
cannot fit every requested cell raise an error instead of silently shrinking
the tissue.

For Python use, call `read_simulation_config(path)` followed by `simulate(config)`
from `tissue_sim.configuration` and `tissue_sim.simulation`, respectively.
`simulate` performs no file I/O and returns a `SimulationResult` containing
DataFrames and parameters. `run(config_path, raw_dir, out_dir)` verifies raw data,
calls `simulate`, and writes outputs. The raw-data check is retained for project
provenance; these supplied files do not drive the simulation.

Standalone grid, placement and transcript commands remain available as file I/O
wrappers. `place_cells(coordinates, config, rng)` now takes a configuration object
and Generator. `generate_transcripts` now takes the complete typed cell DataFrame;
the outdated toy cell table is no longer used.
Cleaning steps belong in data/interim/. Raw data, intermediate data, and
results are excluded from Git; source code, configuration, tests,
documentation, and the checksum manifest are committed.

Tests use temporary inputs and do not require the reference dataset.

## Checks

```bash
pytest        # the tests must pass
ruff check .  # the style checks must pass
```

---

The sections below are stubs. Each is filled in during the week named, and together they are what
makes this repository runnable by someone who has never met you. The long-form explanation of what a
runner does with them is the course document *Running a project from GitHub* — do not copy it in
here; a README is answers, not a manual.

## Run in a container *(Week 11)*

There are two separate things here. **A** happens on GitHub without anyone doing anything; **B** is
optional and is the only part that needs software on your computer. Fill in both.

### A. Built and run on GitHub — nothing to install

Every push, the `build-image` workflow builds this project's image from the `Dockerfile`, publishes
it, and runs the analysis inside it. This happens on GitHub's machines, so it works whether or not
anyone in this group — or anyone reading this — has Docker.

* **Image:** `<ghcr.io/your-group/your-project:tag>`
* **A recent successful build:** `<link to the workflow run>`

<**Where to find both.** On your repository page, the *Packages* panel on the right lists the image;
click it and the page shows the full name and every tag that has been published. The *Actions* tab
lists the workflow runs — open the most recent green `image` run and copy its URL.

The name is not a mystery once you have seen it: it is `ghcr.io/` + the owner and repository, in
lower case, exactly as they appear in your GitHub URL. So `github.com/Ada-Group/spatial-decode`
publishes to `ghcr.io/ada-group/spatial-decode`. The **tag** after the colon is added by the
workflow: `sha-3f9a1c2` for the commit it was built from, and the branch name for the latest build
on that branch. Use the commit tag here — it names one exact version, which `main` does not.>

<If the *Packages* panel is empty, the workflow has not published yet: check the *Actions* tab for a
failed run. If it is there but a classmate cannot pull it, the package is private — open it, then
*Package settings → Change visibility → Public*.>

### B. Running the image yourself — needs Docker or Podman

**Optional.** Only do this if you want the image on your own machine; the project is fully usable
without it, through the *Install* and *Run* sections above. Installing a runtime needs administrator
rights, and on macOS or Windows it runs a Linux virtual machine.

```bash
<docker pull ghcr.io/your-group/your-project:tag>
<docker run --rm \
    -v "$PWD/data:/data:ro" \
    -v "$PWD/results:/results" \
    ghcr.io/your-group/your-project:tag --data-dir /data/raw --out-dir /results>
```

<`podman` takes the same arguments. The two `-v` options are what give the container your data and
somewhere to write; without them it can see neither.>

## Run the whole pipeline *(Week 12)*

```bash
<nextflow run main.nf -profile conda ...>      # conda: no container runtime needed
<nextflow run main.nf -profile docker ...>     # docker: what CI uses, if the reader has it
```

## Verifying you got the right answer *(Week 13–14)*

<The result a correct run produces on the reference dataset — a number, or the checksums of the
output files. Say whether agreement should be exact or to a tolerance, and why.>

## Citation and version *(Week 14)*

<The version this README describes, and the DOI. `CITATION.cff` holds the full record.

In this course the deposit is made on sandbox.zenodo.org, so the DOI starts with 10.5072 and the
record is a practice one — say that here for complete clarity.>

## License *(Week 14)*

<The license your group chose, matching the LICENSE file.>
