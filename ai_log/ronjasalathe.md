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

## 2026-10-14 — Reviewing validation and teammate interfaces

- **Prompt / ask:** Asked AI to review correctness and complexity, check placement/transcript compatibility, and implement the identified validation fixes.
- **What it produced:** Rejected unordered cell-type containers and complex probabilities; added nine test cases and documented teammate handoffs.
- **How you checked it:** AI first reproduced four failing regression cases, then ran all 91 tests with warnings treated as errors. Ruff passed; valid-input outputs and RNG states matched the previous version.
- **What was wrong:** Sets could change the probability-column mapping; complex NumPy probabilities could lose their imaginary parts during conversion. Placement pseudocode did not return actual cell centers, and transcript pseudocode still needed agreement on the shared RNG.
- **Decision:** Accepted the validation fixes. Kept teammate-interface issues documented for coordination rather than changing their code.

---

## 2026-10-14 — Preparing the component for a GitHub pull request

- **Prompt / ask:** Asked AI for terminal/Git instructions and help transferring the component into the group repository.
- **What it produced:** Guided creation of `add-domains` from `origin/main`, preserved the local AI-log commit on a backup branch, and corrected notebook paths and Ruff issues.
- **How you checked it:** In the group’s `pls2026` environment, ran `PYTHONPATH=src python -m pytest -q tests/test_domains.py` — 91 passed. Ran Ruff on the module, tests, and notebook — “All checks passed!” Inspected `git status` and the staged file summary.
- **What was wrong:** Local `main` had diverged; the notebook path was wrong, some `zip()` calls lacked `strict=True`, and a delayed import triggered E402. Adding Python before notebook JSON also made the file invalid.
- **Decision:** Accepted the corrected files and kept the `add-domains` branch. Left unrelated `REAMDE.md` unstaged. Commit/push/merge and successful notebook execution were not confirmed in the chat.

---
