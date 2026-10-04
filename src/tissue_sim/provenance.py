"""Fingerprint input data and record the provenance of a run."""

import datetime
import hashlib
import platform
import subprocess
from pathlib import Path


def manifest(root):
    """Return relative file paths mapped to SHA-256 hashes."""
    root = Path(root)
    hashes = {}

    for path in sorted(root.rglob("*")):
        if path.is_file():
            name = str(path.relative_to(root))
            hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()

    return hashes


def write_manifest(hashes, path):
    """Write hashes in sorted order, using the sha256sum format."""
    path = Path(path)
    path.write_text(
        "".join(
            f"{digest}  {name}\n"
            for name, digest in sorted(hashes.items())
        ),
        encoding="utf-8",
    )


def read_manifest(path):
    """Read a saved manifest into a dictionary."""
    hashes = {}

    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split("  ", 1)
            hashes[name] = digest

    return hashes


def verify(root, manifest_path):
    """Report changed, missing, and added files."""
    expected = read_manifest(manifest_path)
    current = manifest(root)

    expected_files = set(expected)
    current_files = set(current)

    changed = set()
    for name in expected_files & current_files:
        if expected[name] != current[name]:
            changed.add(name)

    return {
        "changed": changed,
        "missing": expected_files - current_files,
        "added": current_files - expected_files,
    }


def require_verified(root, manifest_path):
    """Stop processing if the directory differs from its manifest."""
    report = verify(root, manifest_path)

    if any(report.values()):
        details = "; ".join(
            f"{kind}: {sorted(names)}"
            for kind, names in report.items()
            if names
        )
        raise ValueError(f"Raw-data verification failed: {details}")


def git_commit():
    """Return the code's short Git commit, or None if it is untracked."""
    code_file = Path(__file__).resolve()

    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", code_file.name],
        cwd=code_file.parent,
        capture_output=True,
        text=True,
    )
    if tracked.returncode != 0:
        return None

    result = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        cwd=code_file.parent,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def git_clean():
    """Return True if the code's repository has no uncommitted changes."""
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=Path(__file__).resolve().parent,
        capture_output=True,
        text=True,
    )
    return result.returncode == 0 and result.stdout.strip() == ""


def stamp(config_path, seed, manifest_path):
    """Return a provenance record for one run."""
    return {
        "git_commit": git_commit(),
        "git_clean": git_clean(),
        "config": str(config_path),
        "seed": seed,
        "input_manifest_sha256": hashlib.sha256(
            Path(manifest_path).read_bytes()
        ).hexdigest(),
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "python": platform.python_version(),
    }
