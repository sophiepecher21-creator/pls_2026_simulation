from pathlib import Path
import yaml
import json

from tissue_sim.simulation import run


CONFIG = Path("config/default.yaml")
RAW_DIR = Path("data/raw")


def get_files(directory):
    """Return all generated files and their byte contents."""
    return {
        path.relative_to(directory): path.read_bytes()
        for path in directory.rglob("*")
        if path.is_file()
    }


def comparable_content(filename, content):
    """Ignore only execution timestamps when comparing outputs because these are different for every run."""

    if filename.name not in {"metadata.json", "run_metadata.json"}:
        return content

    metadata = json.loads(content)

    if filename.name == "metadata.json":
        metadata.pop("created_utc", None)

    if filename.name == "run_metadata.json":
        metadata.pop("created", None)

    return json.dumps(metadata, sort_keys=True).encode("utf-8")


def test_reproducibility(tmp_path):
    # 1. Run twice with the same seed
    out1 = tmp_path / "run_001"
    out2 = tmp_path / "run_002"

    run(CONFIG, RAW_DIR, out1)
    run(CONFIG, RAW_DIR, out2)

    files1 = get_files(out1)
    files2 = get_files(out2)

    # Check both runs produced exactly the same files
    assert files1.keys() == files2.keys(), (
        "Same seed produced different sets of output files"
    )

    for filename in files1:
        assert comparable_content(filename, files1[filename]) == comparable_content(
            filename, files2[filename]
        ), f"Same seed produced different contents in {filename}"

    # 2. Create a configuration with a different seed
    config_data = yaml.safe_load(CONFIG.read_text())

    original_seed = config_data["seed"]
    config_data["seed"] = original_seed + 1

    different_config = tmp_path / "different_config.yaml"
    different_config.write_text(yaml.safe_dump(config_data))

    # 3. Run with different seed
    out3 = tmp_path / "run_003"

    run(different_config, RAW_DIR, out3)

    files3 = get_files(out3)

    # Compare all files except metadata that directly records the seed
    simulation_files = [
        filename for filename in files1
        if filename not in {
            Path("metadata.json"),
            Path("run_metadata.json"),
            Path("ground_truth/params.json"),
        }
    ]

    assert any(
        files1[filename] != files3[filename]
        for filename in simulation_files
    ), "Different seeds produced identical simulation outputs"
