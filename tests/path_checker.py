"""README path checker.

Reads README.md, pulls out every file or directory path it mentions, and
asserts each one exists in the repository. This catches the README going
stale: a renamed file, a deleted folder, a typo in a command.

What is deliberately not checked:

* Placeholders inside angle brackets (<...>). These are unfilled stubs for
  later weeks, or text the group still has to replace.
* Paths that are generated or supplied outside Git (GENERATED below). The
  README says raw data, intermediate data and results are never committed.
* Paths owned by a later week (PENDING below). Delete an entry once the
  file is committed, so it is checked from then on.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"

# Not in Git by design: supplied data, intermediate data, program output.
GENERATED = {
    "data/raw",
    "data/interim",
    "results",
    "toy_cells.csv",
    "run_metadata.json",
}

# Mentioned in the README before the week that creates them.
PENDING = {
    "Dockerfile",  # Week 11
}

# A token with no "/" counts as a path only if it looks like a file name.
FILE_SUFFIXES = (
    ".md", ".py", ".yaml", ".yml", ".json", ".csv", ".sha256",
    ".nf", ".cff", ".toml", ".txt", ".cfg",
)
BARE_NAMES = {"Dockerfile", "LICENSE", "Makefile"}


def strip_placeholders(text: str) -> str:
    """Remove <...> spans, including ones that run over several lines."""
    return re.sub(r"<[^<>]*>", " ", text, flags=re.DOTALL)


def extract_paths(text: str) -> set[str]:
    """Return every path-like token in the README outside placeholders."""
    paths = set()
    for token in re.findall(r"[A-Za-z0-9_.\-/]+", strip_placeholders(text)):
        token = token.rstrip(".,;:")
        if not token or token.startswith(("-", "/")) or "//" in token:
            continue
        if "/" in token or token in BARE_NAMES or token.endswith(FILE_SUFFIXES):
            paths.add(token.rstrip("/"))
    return paths


def is_exempt(path: str) -> bool:
    """True for generated or not-yet-created paths."""
    for exempt in GENERATED | PENDING:
        if path == exempt or path.startswith(exempt + "/"):
            return True
    return False


def test_readme_exists():
    assert README.is_file(), "README.md is missing from the repository root"


def test_extractor_finds_paths():
    """Guard against the regex silently matching nothing."""
    sample = "Run `config/default.yaml` and read SPEC.md. <stub/ignored.md>"
    assert extract_paths(sample) == {"config/default.yaml", "SPEC.md"}


def test_every_readme_path_exists():
    found = extract_paths(README.read_text(encoding="utf-8"))
    assert found, "no paths found in README.md; has the extractor broken?"
    missing = sorted(p for p in found if not is_exempt(p) and not (ROOT / p).exists())
    assert not missing, (
        "README.md mentions paths that do not exist: "
        + ", ".join(missing)
        + ". Fix the README or add the file."
    )


@pytest.mark.parametrize("path", sorted(GENERATED))
def test_generated_paths_are_not_committed_data(path):
    """Generated paths are exempt from existence, but must not be tracked."""
    gitignore = ROOT / ".gitignore"
    if not gitignore.is_file():
        pytest.skip("no .gitignore yet")
    if "/" not in path:
        pytest.skip("bare file name; lives under an ignored directory")
    assert path.split("/")[0] in gitignore.read_text(encoding="utf-8"), (
        f"{path} is described as not committed, but .gitignore does not "
        f"mention {path.split('/')[0]}/"
    )