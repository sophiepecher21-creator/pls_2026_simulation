"""Safely load simulation settings into a typed configuration object."""

from dataclasses import dataclass
from pathlib import Path
import math
import yaml


@dataclass(frozen=True)
class SimulationConfig:
    """Settings shared by all simulation stages; lengths are in micrometers."""

    seed: int = 20260101
    n_cells: int = 20
    field_width_um: float = 320.0
    field_height_um: float = 320.0
    square_size_um: float = 2.0
    cell_radius_um: float = 15.0
    distribution_shape: str = "uniform"
    cell_types: tuple[str, ...] = ("type_a", "type_b")
    n_domains: int = 2
    domain_probabilities: tuple[tuple[float, ...], ...] = ((0.8, 0.2), (0.2, 0.8))
    n_genes: int = 100
    n_marker_genes: int = 10
    marker_boost: float = 10.0
    baseline_spread: float = 1.0
    mean_transcripts_per_type: tuple[float, ...] = (20.0, 40.0)

    def __post_init__(self):
        """Reject malformed settings before simulation or random sampling."""
        for name in ("seed", "n_cells", "n_domains", "n_genes", "n_marker_genes"):
            value = getattr(self, name)
            if type(value) is not int or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if min(self.n_cells, self.n_domains, self.n_genes) == 0:
            raise ValueError("n_cells, n_domains and n_genes must be positive")
        for name in (
            "field_width_um",
            "field_height_um",
            "square_size_um",
            "cell_radius_um",
            "marker_boost",
            "baseline_spread",
        ):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ValueError(f"{name} must be a finite number")
            if value < 0 or (name != "baseline_spread" and value == 0):
                raise ValueError(f"Invalid {name}")
        if self.marker_boost < 1:
            raise ValueError("marker_boost must be at least 1")
        for length in (self.field_width_um, self.field_height_um):
            if not (length / self.square_size_um).is_integer():
                raise ValueError("Field dimensions must be multiples of square_size_um")
        if (
            not self.cell_types
            or any(not isinstance(t, str) or not t.strip() for t in self.cell_types)
            or len(set(self.cell_types)) != len(self.cell_types)
        ):
            raise ValueError("cell_types must be unique nonblank names")
        if self.distribution_shape not in ("uniform", "gaussian", "clustered"):
            raise ValueError("Unknown distribution_shape")
        if len(self.domain_probabilities) != self.n_domains:
            raise ValueError("domain_probabilities must have n_domains rows")
        for row in self.domain_probabilities:
            if (
                len(row) != len(self.cell_types)
                or any(
                    type(v) not in (int, float) or not math.isfinite(v) or v < 0
                    for v in row
                )
                or not math.isclose(sum(row), 1, abs_tol=1e-8, rel_tol=0)
            ):
                raise ValueError(
                    "Each domain probability row must match cell_types and sum to 1"
                )
        if len(self.mean_transcripts_per_type) != len(self.cell_types) or any(
            type(v) not in (int, float) or not math.isfinite(v) or v < 0
            for v in self.mean_transcripts_per_type
        ):
            raise ValueError(
                "mean_transcripts_per_type must give a nonnegative mean for each type"
            )
        if self.n_marker_genes * len(self.cell_types) > self.n_genes:
            raise ValueError("Not enough genes for disjoint marker sets")


def read_simulation_config(path):
    """Load plain YAML values without evaluating Python code, then validate them."""
    with Path(path).open(encoding="utf-8") as stream:
        values = yaml.safe_load(stream)
    if not isinstance(values, dict):
        raise ValueError("Configuration must be a YAML mapping")
    values = dict(values)
    for name in ("cell_types", "mean_transcripts_per_type"):
        if name in values:
            if not isinstance(values[name], list):
                raise ValueError(f"{name} must be a YAML list")
            values[name] = tuple(values[name])
    if "domain_probabilities" in values:
        rows = values["domain_probabilities"]
        if not isinstance(rows, list) or any(not isinstance(row, list) for row in rows):
            raise ValueError("domain_probabilities must be a list of lists")
        values["domain_probabilities"] = tuple(tuple(row) for row in rows)
    return SimulationConfig(**values)
