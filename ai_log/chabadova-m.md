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

## <date> — <what I was building>
- Prompt / question: <what I asked the AI or wrote myself>
- What it produced: <a line or two>
- How I checked it: <the specific check I ran>
- What was wrong: <bugs, hallucinations, or "nothing found">
- Decision: <accepted / edited / rejected — and why>This is a test, it should not require approval
