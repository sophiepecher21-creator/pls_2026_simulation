"""Check README file references using the Week 6 notebook's path convention."""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PATHY = re.compile(r"^[\w./-]+\.(py|ya?ml|md|csv|toml|json|pdf)$")


def missing_paths(readme_text, root):
    """Return missing inline, linked, and command-block file paths.

    This mechanical check cannot prove that commands run or replace a classmate
    following the README. Paths are relative to the repository root.
    """
    inline = re.sub(r"```.*?```", "", readme_text, flags=re.DOTALL)
    candidates = set(re.findall(r"`([^`\n]+)`", inline))
    candidates.update(re.findall(r"\[[^\]]+\]\(([^)]+)\)", readme_text))
    for block in re.findall(r"```[^\n]*\n(.*?)```", readme_text, re.DOTALL):
        candidates.update(block.split())
    return sorted(
        token
        for token in candidates
        if PATHY.fullmatch(token) and not (root / token).is_file()
    )


def test_readme_paths_exist():
    assert missing_paths((ROOT / "README.md").read_text(), ROOT) == []


def test_path_checker_known_answer(tmp_path):
    (tmp_path / "present.py").touch()
    assert missing_paths(
        "`present.py` `missing.csv` `python -m tissue_sim`", tmp_path
    ) == ["missing.csv"]


def test_checker_catches_command_and_link_paths(tmp_path):
    text = "[guide](missing.md)\n```bash\npython --config config/absent.yaml\n```"
    assert missing_paths(text, tmp_path) == ["config/absent.yaml", "missing.md"]
