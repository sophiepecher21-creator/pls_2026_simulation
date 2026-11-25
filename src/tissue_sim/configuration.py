"""Safely load simulation settings into a typed configuration object."""

from dataclasses import dataclass, fields
from pathlib import Path
import difflib
import math
import yaml


DEFAULT_CONFIG_PATH = Path("config/default.yaml")


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
    gaussian_spread_divisor: float = 6.0
    cluster_spread_divisor: float = 15.0
    min_clusters: int = 2
    cluster_count_divisor: float = 2.0
    max_attempts_per_cell: int = 200
    cell_types: tuple[str, ...] = ("type_a", "type_b")
    n_domains: int = 2
    domain_probabilities: tuple[tuple[float, ...], ...] = ((0.8, 0.2), (0.2, 0.8))
    reference_n_per_type: int = 80
    sequencing_depth: float = 1.0
    n_genes: int = 100
    n_marker_genes: int = 10
    marker_boost: float = 10.0
    baseline_spread: float = 1.0
    mean_transcripts_per_type: tuple[float, ...] = (20.0, 40.0)

    def __post_init__(self):
        """Reject malformed settings before simulation or random sampling."""
        for name in (
            "seed",
            "n_cells",
            "n_domains",
            "n_genes",
            "n_marker_genes",
            "reference_n_per_type",
            "min_clusters",
            "max_attempts_per_cell",
        ):
            value = getattr(self, name)
            if type(value) is not int or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")

        if self.min_clusters < 1 or self.max_attempts_per_cell < 1:
            raise ValueError(
                "min_clusters and max_attempts_per_cell must be positive"
            )

        if (
            min(self.n_cells, self.n_domains, self.n_genes, self.reference_n_per_type)
            == 0
        ):
            raise ValueError(
                "n_cells, n_domains, n_genes and reference_n_per_type must be positive"
            )
        for name in (
            "field_width_um",
            "field_height_um",
            "square_size_um",
            "cell_radius_um",
            "marker_boost",
            "baseline_spread",
            "sequencing_depth",
            "gaussian_spread_divisor",
            "cluster_spread_divisor",
            "cluster_count_divisor",
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
        if self.square_size_um != 2:
            raise ValueError("Data Contract 1.0 requires square_size_um == 2")
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
            raise ValueError(
                f"Unknown distribution_shape: {self.distribution_shape!r}. "
                "Choose uniform, gaussian or clustered."
            )
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

def check_keys(cfg, valid_keys, filename):
    """Reject unknown configuration keys and suggest corrections."""

    for key in cfg:
        if key not in valid_keys:
            suggestions = difflib.get_close_matches(
                str(key), sorted(valid_keys), n=1
            )

            message = f"{filename}: Unknown configuration key {key!r}"

            if suggestions:
                message += f" — did you mean {suggestions[0]!r}?"

            raise ValueError(message)

def read_simulation_config(path, overrides=None):
    """Return validated SimulationConfig settings from a YAML file path.

    Physical lengths are micrometres; means are captured transcripts per cell.
    These are model inputs, not quantities estimated from reference data.
    Only plain YAML values are accepted, never evaluated Python expressions.
    """
    with Path(path).open(encoding="utf-8") as stream:
        values = yaml.safe_load(stream)
    if not isinstance(values, dict):
        raise ValueError(f"{path}: Configuration must be a YAML mapping")
    values = dict(values)
    valid_keys = {field.name for field in fields(SimulationConfig)}

    check_keys(values, valid_keys, path)

    missing_keys = valid_keys - values.keys()
    if missing_keys:
        raise ValueError(
            f"{path}: Missing required configuration keys: {sorted(missing_keys)}"
        )

    if overrides:
        check_keys(overrides, valid_keys, path)
        values.update(overrides)

    if overrides:
        values.update(overrides)

    for name in ("cell_types", "mean_transcripts_per_type"):
        if name in values:
            if not isinstance(values[name], list):
                raise ValueError(f"{path}: {name} must be a YAML list")
            values[name] = tuple(values[name])
    if "domain_probabilities" in values:
        rows = values["domain_probabilities"]
        if not isinstance(rows, list) or any(not isinstance(row, list) for row in rows):
            raise ValueError(f"{path}: domain_probabilities must be a list of lists")
        values["domain_probabilities"] = tuple(tuple(row) for row in rows)
    try:
        config = SimulationConfig(**values)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{path}: {exc}") from exc

    return config
