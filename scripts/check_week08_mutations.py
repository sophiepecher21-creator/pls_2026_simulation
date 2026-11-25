"""Run focused Week 8 mutation checks in disposable source copies."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TEST = "tests/test_simulation_invariants.py"

MUTATIONS = [
    (
        "Subtract instead of accumulating transcript counts",
        "src/tissue_sim/aggregation.py",
        "np.add.at(matrix, (rows, cols), 1)",
        "np.add.at(matrix, (rows, cols), -1)",
    ),
    (
        "Reverse the square-to-matrix-row mapping",
        "src/tissue_sim/aggregation.py",
        "enumerate(ids)",
        "enumerate(reversed(ids))",
    ),
    (
        "Ignore the configured seed",
        "src/tissue_sim/simulation.py",
        "np.random.default_rng(config.seed)",
        "np.random.default_rng(42)",
    ),
]


def run_case(directory, mutation=None):
    """Run the focused tests against an isolated copy of the current source."""
    shutil.copytree(ROOT / "src", directory / "src")
    (directory / "tests").mkdir()
    shutil.copy2(ROOT / TEST, directory / TEST)

    if mutation is not None:
        _, filename, before, after = mutation
        path = directory / filename
        source = path.read_text(encoding="utf-8")
        if before not in source:
            raise ValueError(f"Mutation pattern missing in {filename}: {before}")
        path.write_text(source.replace(before, after, 1), encoding="utf-8")

    environment = {
        **os.environ,
        "PYTHONPATH": str(directory / "src"),
        "PYTHONDONTWRITEBYTECODE": "1",
        "NO_COLOR": "1",
    }
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "-p",
            "no:cacheprovider",
            TEST,
        ],
        cwd=directory,
        env=environment,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    print(output)
    return result.returncode, output


def main():
    """Record the baseline and every mutant, including survivors and errors."""
    lines = [
        "# Week 8 mutation record",
        "",
        "Command: `python scripts/check_week08_mutations.py`.",
        "",
        "Scope: the seven cases in "
        "`tests/test_simulation_invariants.py`, run against temporary "
        "copies of the current source. This is a focused mutation experiment, "
        "not a full-suite run or an exhaustive mutation analysis.",
        "",
    ]
    unsuccessful = False

    with tempfile.TemporaryDirectory(prefix="pls-week08-mutations-") as temporary:
        workspace = Path(temporary)
        baseline = workspace / "baseline"
        baseline.mkdir()
        print("BASELINE")
        code, _ = run_case(baseline)
        if code != 0:
            raise SystemExit("Baseline failed; fix it before interpreting mutations.")

        lines.extend(
            [
                "Baseline: PASS (pytest exit code 0).",
                "",
                "| Mutation | Result | Failing tests |",
                "| --- | --- | --- |",
            ]
        )

        for index, mutation in enumerate(MUTATIONS):
            directory = workspace / f"mutation-{index}"
            directory.mkdir()
            print(f"\nMUTATION: {mutation[0]}")
            code, output = run_case(directory, mutation)
            failures = [
                line.split(" - ", 1)[0].removeprefix("FAILED ")
                for line in output.splitlines()
                if line.startswith("FAILED ")
            ]
            if code == 1 and failures:
                status = "CAUGHT"
            elif code == 0:
                status = "SURVIVED"
            else:
                status = f"ERROR (exit {code}; inspect terminal output)"
            unsuccessful |= status != "CAUGHT"
            names = "; ".join(f"`{name}`" for name in failures) or "None"
            lines.append(f"| {mutation[0]} | {status} | {names} |")
            print(status)

    lines.extend(
        [
            "",
            "Every case starts from a fresh copy. The working source is never "
            "mutated. A caught mutant demonstrates protection against that "
            "specific change; it does not establish complete correctness.",
            "",
        ]
    )
    report = ROOT / "docs/week08_mutations.md"
    report.parent.mkdir(exist_ok=True)
    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {report}")
    if unsuccessful:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
