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

## <date> — <what I was building>
- Prompt / question: <what I asked the AI or wrote myself>
- What it produced: <a line or two>
- How I checked it: <the specific check I ran>
- What was wrong: <bugs, hallucinations, or "nothing found">
- Decision: <accepted / edited / rejected — and why>This is a test, it should not require approval
