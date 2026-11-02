## Week 3 — 2026-10-04

Used Codex to understand manifest(), verify(), and stamp(), then to guide
terminal commands and provide code for the project structure, provenance
module, Week 3 toy simulation, and tests. I ran the commands and edited
the files myself.

Copied the supplied reference dataset into data/raw/, excluding Finder
metadata and notebook checkpoints, and removed write permissions.
Created a SHA-256 manifest covering 10 files.

Verified the package installation and read-only permissions. Deliberately
disabled changed-file detection: two tests failed for the expected reasons.
Restored it and confirmed all seven tests pass. ruff check . also passes.
The tests cover changed, missing, and added inputs, refusal to process a
changed input, end-to-end output, and deletion followed by reproducible
rebuilding of the cell table.

The current simulator is a two-domain, two-cell-type toy, not yet the full
Data Contract implementation. The supplied ground truth is given data;
our own future simulation outputs will be generated data.


## Week 4 2026-10-13

- Prompt:  Help split the simulation into modular functions and
  implement a measurement grid, CSV writer, and configuration-driven workflow.
- What it produced: make_grid(), write_coordinates(), a YAML reader,
  a terminal entry point, and tests.
- How was it checked:  Known-answer coordinates copied from the instructor's
  CSV; invalid-input tests; CSV content and overwrite-protection tests;
  workflow test using temporary files. All 18 repository tests passed.
  Ran the workflow and inspected the generated CSV header, coordinates,
  and line count.
- What was wrong: Initially used two-digit rather than three-digit square
  IDs; corrected them to match the reference. Resolved test file location
  and module import issues. Also over-complicated the script for where we are
  in the project right now - asked it to be more in alignment with what we have
  worked on in our .ipynb modules so far. 
- Decision: Accepted after verification. Used an exact whole-square check
  and exclusive file creation to protect existing outputs.

26-10-14
- Prompt: Modify the read_config() function so that the check to ensure the file 
  is read correctly is actually happening when we run the function. Also pass the file
  as an object, so that when we need to reference values from the config file later
  we do not need to type strings. We can just access the values through attributes. 
- What it produced: A modified read_config() function which uses dataclass to define 
  a container which names attributes that we will want to use later when running
  the make_grid() function. 
- How it was checked: there is a configuration reader test which returns validated
  settings from the .yaml file. This test asserts that the function is able to read variables from the 
  file and return the expected numerical values associated with them.
- How it was wrong: This was a great improvement, but the previous commit was far too
  modular and redundant. The previous make_grid.py and write_coordinates.py files were 
  depricated; everything now lives in one file and can be run with one wrapper function main(). 
- Decision: Kept the AI improvements, and also the oranigzational modifications that I made. 

## Weeks 6 and 7 — 2026-11-02

- Request: implement executable documentation, README path checks, the simulation
  conformance checklist and validate_dataset(path), using both weekly handouts
  and Data Contract 1.0; install the supplied broken dataset pack locally.
- Codex produced: two function doctests, clarified units/provenance in touched
  public docstrings, a shorter README, automated path checks, a validator with
  actionable ValueError messages, independently sampled reference cells, full
  contract metadata, validation at export, and generated-fixture and six course
  dataset regression tests. Added the provided contract as a local documentation
  reference. No course dataset is committed.
- How it was checked: pytest --doctest-modules: 170 passed; Ruff passed. The original
  reference dataset validates silently. The five faulty datasets fail for their
  specific defects; dataset_OK is silent. A default end-to-end run preserves the
  earlier 25,600 squares, 20 tissue cells and 613 transcripts, and now adds 80
  independent reference cells per type. All 10 original raw-file hashes remain
  unchanged; the updated manifest additionally fingerprints the 48 pack files.
- How it was wrong: removed trailing whitespace in a doctest's
  expected output; specified zip truncation explicitly in mismatched-order error
  reporting; excluded contract timestamps from byte-for-byte repeat comparisons.
- Decision: changes prepared and verified locally for team review and no commit or push made by this task.


### Week 7 notebook alignment — 2026-11-02

- Prompt: validate_dataset(path)should be based on the supplied
  week07_notebook.ipynb. Reorganized it around the notebook's read_csv,
  check_counts, check_ids_match and check_composition_sums functions, using
  csv.DictReader and explicit dictionary/row iteration instead of pandas.
- What it produced: Kept the full contract extensions for metadata, geometry, independent
  reference tables, sparse composition representation, domains and parameters.
  Required files are checked before inspecting their contents. Invalid input
  continues to raise actionable ValueError messages, preserving the package API. Added one focused regression test reproducing the notebook's partial 3x3
  fixture, explicit zero composition rows, negative count, extra square ID and
  the 12-versus-9 composition mismatch. Helpers accept the valid partial fixture;
  full validation still requires every contract file and full metadata.
- How it was checked: Verified the existing tests against the refactor and the new notebook example.
- Decision: implemented
