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
