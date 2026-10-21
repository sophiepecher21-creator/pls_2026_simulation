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
- Prompt / question: Helped me implement all necessary functions for my python file, validate all parameters. 
- What it produced: Code which a carefully read and edited if not working correctly.
- How I checked it: run several tests for different types of possible problems.
- What was wrong: firstly, marker genes of different cell types could overlap
- Decision: edited, made sure that marker genes don't overlap between cell types.

## 21-10-2026 — check high expression of non-marker genes
- Prompt / question: How can we check if these are significantly higher genes compared to others?
- What it produced: Code calculating mean and median values for each cell type and comparing to non-marker genes using ratio.
- How I checked it: Mean ratio of marker genes was much higher to non-marker genes and also between different cell types.
- What was wrong: Printing top 20 genes for each cell type, some of these top genes were not marker-genes, they had higher expression than marker genes, but these high expressed non-marker genes were same for all cell types.
- Decision: accepted, because although non-marker have high expression, they are equally expressed in all the types, which is not a problem, because then they are not considered as marker-genes. 

## <date> — <what I was building>
- Prompt / question: <what I asked the AI or wrote myself>
- What it produced: <a line or two>
- How I checked it: <the specific check I ran>
- What was wrong: <bugs, hallucinations, or "nothing found">
- Decision: <accepted / edited / rejected — and why>This is a test, it should not require approval
