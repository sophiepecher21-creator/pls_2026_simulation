"""One regression test per course dataset (external, read-only fixtures)."""

from pathlib import Path

import pytest

from tissue_sim.validate import validate_dataset

PACK = Path(__file__).resolve().parents[1] / "data/raw/broken_datasets"


@pytest.fixture
def pack():
    """Skip external fixtures on clones where the course pack is not installed."""
    if not PACK.is_dir():
        pytest.skip(
            "Install the course pack under data/raw/broken_datasets to run these six tests"
        )
    return PACK


def test_dataset_A(pack):
    with pytest.raises(
        ValueError,
        match=r"counts.csv:.*sq_r017_c023.*GENE_0007='-3'.*non-negative integer",
    ):
        validate_dataset(pack / "dataset_A")


def test_dataset_B(pack):
    with pytest.raises(
        ValueError, match=r"composition.csv:.*sq_r000_c001.*sum=7.*total=2"
    ):
        validate_dataset(pack / "dataset_B")


def test_dataset_C(pack):
    with pytest.raises(
        ValueError, match=r"counts.csv: missing required columns.*square_id"
    ):
        validate_dataset(pack / "dataset_C")


def test_dataset_D(pack):
    with pytest.raises(
        ValueError, match=r"reference_labels.csv: required readable CSV file"
    ):
        validate_dataset(pack / "dataset_D")


def test_dataset_E(pack):
    with pytest.raises(
        ValueError, match=r"counts.csv: square_id mismatch.*sq_r031_c008.*sq_r131_c008"
    ):
        validate_dataset(pack / "dataset_E")


def test_dataset_OK(pack, capsys):
    assert validate_dataset(pack / "dataset_OK") is None
    assert capsys.readouterr() == ("", "")
