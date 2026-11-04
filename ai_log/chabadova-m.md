## 07-10-2026 — manifest() function
- Prompt: What is the purpose of a SHA-256 checksum manifest? Can raw data be regenerated from the hash?
- What it produced: Explained that SHA-256 creates a fingerprint of each raw file. The manifest stores these hashes so the current raw data can be compared with the original version. It cannot regenerate the original data.
- How I checked it: Implemented manifest() in the notebook. The notebook's known-answer check compared the generated file list and counts.csv hash.
- What was wrong: I used RAW/f, but f was already a full path. Then the manifest stored Path objects instead of strings, so got.get("counts.csv") returned None.
- Decision: Edited/accepted — fixed the path handling and converted relative paths to strings.

## 07-10-2026 — PYTHONPATH meaning
- Prompt: What does PYTHONPATH mean, and do I need it in my own project?
- What it produced: Explained that the notebook uses a temporary project that is not installed, so it temporarily adds src/ to PYTHONPATH to allow from spatial_decode.cleaning import clean. In the real project, we just use pip install -e .
- Decision: Accepted — PYTHONPATH is only a temporary workaround for the notebook.

## 14-10-2026 — read_globals()
- Prompt / question: Correct my implementation for read_globals()
- What it produced: corrected implementation with line ```value = func.__globals__[name]```
- How I checked it: I ran the code but there was an error
- What was wrong: ```func.__globals__[name]``` has no items
- Decision: rejected because ```__globals__``` is not a dict

## 20-10-2026 — expression_programs.py
- Prompt / question: Used AI to assist with implementing and reviewing the expression-program functions.
- What it produced: An implementation with input validation and tests for the expression-program functions.
- How I checked it: Reviewed the implementation and checked the marker-selection and expression-probability behavior. All 24 pytests passed.
- What was wrong: The initial implementation allowed marker genes to overlap between different cell types.
- Decision: Edited the implementation to ensure that marker gene sets do not overlap between cell types.

## 21-10-2026 — inspecting expression probabilities of marker and non-marker genes
- Prompt / question: Asked how to assess whether marker genes have higher expression probabilities than non-marker genes.
- What it produced: Code calculating mean and median probabilities for marker and non-marker genes, their ratio and comparisons between cell-type expression programs.
- How I checked it: Compared the mean probabilities of marker and non-marker genes and examined the highest-ranked genes in each cell type. I also compared marker-gene probabilities between their designated cell type and the other cell type.
- What was wrong: Some highly ranked genes were not markers. These genes had high baseline probabilities and were highly expressed in both cell types.
- Decision: Accepted the results because marker status is determined by the selected marker sets, not by expression level alone. High-probability non-marker genes are consistent with the model because they arise from the shared baseline distribution rather than being specifically boosted for one cell type.

## 04-11-2026 — Raw data verification
- Prompt / question: My simulation failed because the raw-data verification detected unexpected `.DS_Store` files.
- What it produced: AI explained that these files are automatically created by macOS and suggested removing them.
- How I checked it: I removed the `.DS_Store` files from `data/raw` and ran the simulation again.
- What was wrong: The additional files were not included in the raw-data manifest.
- Decision: Accepted the suggestion.

## 04-11-2026 — Reproducibility test
- Prompt / question: I needed to implement a test to check whether the same seed produces identical simulation outputs and a different seed produces different results.
- What it produced: AI suggested a pytest test using `tmp_path`, running the simulation three times and comparing the generated files.
- How I checked it: I ran `python -m pytest tests/test_reproducibility.py -v`. The test passed after I adjusted the comparison.
- What was wrong: The first version failed because `metadata.json` and `run_metadata.json` contain timestamps that change between runs. I also needed to make sure the different-seed test checks actual simulation results, not just the seed stored in metadata.
- Decision: Edited the test to ignore timestamp fields when comparing runs with the same seed. For different seeds, I excluded metadata and parameter files and checked that at least one simulated data file differs. Kept the actual timestamps in the output.

## 04-11-2026 — Changing simulation parameters
- Prompt / question: I wanted to check whether our cell radius, number of cells and number of transcripts were reasonable based on the spatial transcriptomics course material.
- What it produced: AI helped me calculate cell size, tissue coverage and expected transcript counts for different parameter values.
- How I checked it: I tried different parameter values, ran the simulation and checked the generated `counts.csv` using pandas. With 300 cells and `mean_transcripts_per_type: [600.0, 1200.0]`, the simulation produced 264,382 transcripts across 25,600 squares. The mean was 10.33, median 8, and 38.18% of squares were empty.
- What was wrong: The original radius of 15 µm meant a diameter of 30 µm, which was too large compared to the typical cell diameter of 10–20 µm mentioned in the course material. Also, 1,000 cells were too many for the grid with non-overlapping cells, and the original transcript counts were very low compared to the example of approximately 11 median transcripts per 2 µm square.
- Decision: Changed `cell_radius_um` to 7.5, `n_cells` to 300 and `mean_transcripts_per_type` to `[600.0, 1200.0]`. These values seemed like a reasonable starting point for the simulation, although further testing may be needed.

## <date> — <what I was building>
- Prompt / question: <what I asked the AI or wrote myself>
- What it produced: <a line or two>
- How I checked it: <the specific check I ran>
- What was wrong: <bugs, hallucinations, or "nothing found">
- Decision: <accepted / edited / rejected — and why>This is a test, it should not require approval
