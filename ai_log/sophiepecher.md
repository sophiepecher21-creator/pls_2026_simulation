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

