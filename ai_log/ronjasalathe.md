## <date> — <what you were building>

- **Prompt / ask:** <what you asked, or "wrote by hand">
- **What it produced:** <a line or two>
- **How you checked it:** <the specific check you ran — not "I read it">
- **What was wrong:** <bugs, hallucinations, or "nothing found">
- **Decision:** <accepted / edited / rejected — and why>

---

## 27.09.2026 — SPEC.md "What it does"

- **Prompt / ask:** Used ChatGPT and Notebook LM to understand the simulation track specification (review acceptance criteria) and drafting the "What it does" section of SPEC.md.
- **What it produced:** A comprehensive feedback on the acceptance criteria and what could be improved. A short description of the tissue-sim explaining the simulated 2D tissue, stochastic gene expression, ground truth and reference data.
- **How you checked it:** Worked through the documents provided - 00_data_contract.pdf, 01_simulation_track.pdf, 06_git_reference.pdf.
- **What was wrong:** Correction of the index of the grid. 
- **Decision:** Edited and commented the Pull request of Sophie. Submitted a pull request for the "What it does" section for the SPEC.md document.

---

## 14.10.2026 — Organizing and understanding the domain component

- **Prompt / ask:** Asked AI to explain every function, simplify the file structure, and preserve the results.
- **What it produced:** Grouped three main functions and five short helpers in `domains.py`; updated `START_HERE_DOMAINS.md` and illustrated the cell and square outputs.
- **How you checked it:** AI ran 82 tests and Ruff checks, compared outputs across 1,600 positions and four geometries, and checked that output tables and RNG states matched before and after refactoring.
- **What was wrong:** Separate files made the component harder to navigate; putting all validation inside the main functions made them unnecessarily long.
- **Decision:** Accepted the single-module structure with short helpers because it is easier to understand and preserves behavior.
  
---

## 14.10.2026 — Reviewing validation and teammate interfaces

- **Prompt / ask:** Asked AI to review correctness and complexity, check placement/transcript compatibility, and implement the identified validation fixes.
- **What it produced:** Rejected unordered cell-type containers and complex probabilities; added nine test cases and documented teammate handoffs.
- **How you checked it:** AI first reproduced four failing regression cases, then ran all 91 tests with warnings treated as errors. Ruff passed; valid-input outputs and RNG states matched the previous version.
- **What was wrong:** Sets could change the probability-column mapping; complex NumPy probabilities could lose their imaginary parts during conversion. Placement pseudocode did not return actual cell centers, and transcript pseudocode still needed agreement on the shared RNG.
- **Decision:** Accepted the validation fixes. Kept teammate-interface issues documented for coordination rather than changing their code.

---

## 14.10.2026 — Preparing the component for a GitHub pull request

- **Prompt / ask:** Asked AI for terminal/Git instructions and help transferring the component into the group repository.
- **What it produced:** Guided creation of `add-domains` from `origin/main`, preserved the local AI-log commit on a backup branch, and corrected notebook paths and Ruff issues.
- **How you checked it:** In the group’s `pls2026` environment, ran `PYTHONPATH=src python -m pytest -q tests/test_domains.py` — 91 passed. Ran Ruff on the module, tests, and notebook — “All checks passed!” Inspected `git status` and the staged file summary.
- **What was wrong:** Local `main` had diverged; the notebook path was wrong, some `zip()` calls lacked `strict=True`, and a delayed import triggered E402. Adding Python before notebook JSON also made the file invalid.
- **Decision:** Accepted the corrected files and kept the `add-domains` branch. Left unrelated `REAMDE.md` unstaged. Commit/push/merge and successful notebook execution were not confirmed in the chat.

---

## 26.10.2026 — Reviewing requested changes for domains.py code on GitHub

- **Prompt / ask:** Checked domains.py code for correct interface to cell-placement code to have aligned input/output. Before change, not all tests passed.
- **What it produced:** Small corrections on where the code imports the input from cell-placement. Highlighted where x_um has to be changed to center_x_um (same for y_um).
- **How you checked it:** In group's pls2026 environment, ran `PYTHONPATH=src python -m pytest -q tests/test_domains.py` — 92 passed. Inspected `git status` and staged the file summary.
- **What was wrong:** AI earlier suggested to use `x_um` and `y_um` for the cell centers, which conflicted with cell-placement name conventions. Cell centers needed `center_x_um` and `center_y_um`. The tests initially used old names. Some suggested documentary was unnecessary long.
- **Decision:** Edited individual codelines and accepted. Kept the cell-placement column names consistent, shortened explanations, and updated the affected tests before committing.

---

## 26.10.2026 — Change output of domains.py code

- **Prompt / ask:** Change domains.py to add domain and cell-type annotations to the original coordinates.csv, preserving its existing columns instead of creating a new file.
- **What it produced:** Updated the function `domains_from_csv()` to overwrite the input CSV with added domain and cell_type columns. 
- **How you checked it:** Ran `PYTHONPATH=src python -m pytest -q tests/test_domains.py` — 92 passed. The new CSV test `test_coordinates_csv_is_updated_in_place` verified that original columns and row order were preserved, annotations were correct, and no additional CSV was created. 
- **What was wrong:** Intermediate AI suggestion returned square_domains without defining it. The was corrected before the code was accepted and passed the tests after corrections.
- **Decision:** Edited and accepted. Updating the original CSV preserved the existing data instead of copying it, and providing the added annotation for downstream use.

---

## 28.10.2026 — Simplifying domain tests file

- **Prompt / ask:** Review `domains.py` and `test_domains.py`, then simplify the test file. 
- **What it produced:** A shorter file with 11 tests covering domain assignment, cell extraction, sampling, reproducibility, and DataFrame integration.
- **How you checked it:** Ran both `PYTHONPATH=src python -m pytest -q tests/test_domains.py` (11 tests passed, 100 %) and `PYTHONPATH=src python -m pytest -q test` (54 passed, 100 %).
- **What was wrong:** The first proposal removed too many useful checks and included unnecessary  CSV tests. Both points were corrected.
- **Decision:** Accepted the revised proposal because it checks the main functions without extensively repeating upstream input validation.

---

## 03.11.2026 — Reviewing and testing the simulation pipeline

- **Prompt / ask:** Review the connected pipeline and tests, and explain how to test the unmerged changes.
- **What it produced:** Review comments with line numbers and terminal commands for checking the pipeline branch.
- **How you checked it:** Ran ruff check . Successfully and pytest --doctest-modules: 165 passed, 6 skipped. The full simulation stopped because the course datasets were missing.
- **What was wrong:** Some suggestions were optional rather than required. One proposed regex correction was identical to the original and was withdrawn.
- **Decision:** Accepted the testing instructions and assessed review comments individually. Full pipeline verification remains incomplete until the course datasets are installed.

---

## 04.11.2026 — Acceptance criteria test mapping in SPEC.md

- **Prompt / ask:** Helped me add a table to SPEC.md with one row per acceptance criterion, showing the claim, the kind of check and test file and function.
- **What it produced:** A Markdown table mapping the seven Weeks 6-7 appearance criteria to existing tests, with notes about coverage gaps.
- **How you checked it:** Used `git grep` and `sed` to identify dockets functions and supplied test code for mapping. Ran `pytest --doctest-modules`: 84 passed, 6 external-fixture tests skipped. Ran `git diff --check` with no errors and inspected `git diff --SPED.md`.
- **What was wrong:** The first draft lacked the exact doctest function names. Some criteria were only partly covered by the existing tests; passing tests did not establish every claim.
- **Decision:** Edited and accepted. Added the confirmed function and kept coverage gaps explicit. Committed and pushed the SPEC.md change for teammate review

---

## 15.11.2026 — Week 9 Debugging and Regression Testing

- **Prompt / ask:** Help identify and fix five bugs in the spatial_mini pipeline and write regression tests.
- **What it produced:** Explanations of bugs, suggested code fixes, and pytest regression tests.
- **How you checked it:** Ran `pytest -q` (9/9 passed) and compared pipeline output with expected.json.
- **What was wrong:** Five bugs in the original code; AI suggestions needed testing but no further issues found after corrections.
- **Decision:** Accepted and implemented fixes after verifying the tests and expected output.

---

## 23.11.2026 — Week 10: Configuration, CLI and reproducibility

- **Prompt / ask:** Used AI to help me implement the Week 10 requirements in
  the existing simulation pipeline, while keeping
  changes minimal and preserving the existing scientific defaults.

- **What it produced:** Suggestions for moving hard-coded parameters
  into YAML, configuration validation, CLI overrides, run metadata,
  reference sequencing depth, expression-program correlation,
  automated tests and documentation.

- **How you checked it:**
  - Ran pytest after individual implementation steps (96 tests passed).
  - Tested the CLI using --help and a simulation with parameter overrides.
  - Inspected the generated run_metadata.json.
  - Verified that increasing sequencing depth affects reference counts
    without changing spatial counts.
  - Added six test cases for validation, configuration precedence,
    sequencing depth and CLI reproducibility (102 tests passed).
  - Checked test runtime using pytest --durations=6.
  - Ran documentation tests (3 passed) and git diff --check.

- **What was wrong:**
  - The initial configuration validation order caused an existing
    test to fail because it reported missing keys before an invalid seed.
  - The CLI reproducibility test used an unnecessarily large dataset
    and took approximately 15 seconds.
  - The README path checker rejected references to generated output
    files because they did not exist in the repository.

- **Decision:** Edited and accepted. Adjusted the validation order,
  reduced the integration test dataset size and corrected the README
  formatting. Applied changes manually and verified them incrementally.
  Preserved the existing scientific defaults and restricted sequencing
  depth to the independent reference experiment. The optimized Week 10
  tests all passed in 3.57 seconds. Final verification: 104 tests passed, including doctests.
  Ruff and git diff --check reported no issues.

---