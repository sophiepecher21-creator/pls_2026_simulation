"""Spatial domains and cell-type assignment.

Input: Cell-placement DataFrames or coordinates.csv.
Use domains_from_csv() to load the CSV and annotate cells and squares.

Square domains:
    build_domains_dataframe() uses square_id, x_um and y_um.

Cell annotations:
    cells_from_placement() extracts unique cells.
    Cell centres remain named center_x_um and center_y_um.
    assign_cell_types() assigns cell domains and cell types.

Empty squares have cell_id=0.
Prefer the complete cell-placement table to include every placed cell.
"""


from collections.abc import Mapping, Sequence
from collections.abc import Set as AbstractSet
from numbers import Real

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

__all__ = [
	"assign_cell_types", 
	"build_domains_dataframe", 
	"get_domain_id",
	"cells_from_placement",
	"domains_from_csv",
]

def cells_from_placement(placement_df: pd.DataFrame) -> pd.DataFrame:
    """Prepare cell-placement data.

    Parameters
    ----------
    placement_df : pandas.DataFrame
        coordinates.csv, already loaded as a DataFrame.
        Alternatively, a complete cell-placement DataFrame.
        Required columns:
        cell_id, center_x_um, center_y_um, radius_um.

    Returns
    -------
    pandas.DataFrame
        One row per cell.
        Columns: cell_id, center_x_um, center_y_um, radius_um.
        Cell-centre column names are preserved.
        Empty squares (cell_id=0) are excluded.
        Pass this result to assign_cell_types().

    Raises
    ------
    KeyError
        If a required column is missing.
    ValueError
        If one cell_id has conflicting centres or radii.
    """   


    cells = (
        placement_df[
            ["cell_id", "center_x_um", "center_y_um", "radius_um"]
        ]
        .loc[lambda df: ~df["cell_id"].isin([0, "0"])]
        .dropna(subset=["cell_id"])
        .drop_duplicates()
    )

    if not cells["cell_id"].is_unique:
        raise ValueError("Conflicting centres or radii for the same cell_id")

    return cells


def domains_from_csv(
    path, K, field_width_um, cell_types, domain_probabilities, rng
):
    """Read Sven's CSV. Return annotated cells and square domains.

    Cell coordinates: center_x_um, center_y_um.
    Square output: square_id, domain.
    """
    grid = pd.read_csv(
        path, dtype={"square_id": "string", "cell_id": "string"}
    )
    square_domains = build_domains_dataframe(grid, K, field_width_um)
    cells = assign_cell_types(
        cells_from_placement(grid),
        K, field_width_um, cell_types, domain_probabilities, rng,
    )

    return cells, square_domains



def _finite_number(value: float, name: str) -> float:
    """Return a finite real value as a float, rejecting invalid scalar inputs."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    try:
        number = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not np.isfinite(number):
        raise ValueError(f"{name} must be a finite real number")
    return number


def _validate_position(
    x: float, y: float, K: int, field_width_um: float
) -> tuple[float, float, float]:
    """Validate geometry and position, returning width, x, and y as floats."""
    if isinstance(K, (bool, np.bool_)) or not isinstance(K, (int, np.integer)):
        raise TypeError("K must be a positive integer")
    if K <= 0:
        raise ValueError("K must be a positive integer")
    width = _finite_number(field_width_um, "field_width_um")
    x = _finite_number(x, "x")
    y = _finite_number(y, "y")
    if width <= 0:
        raise ValueError("field_width_um must be positive")
    if not 0 <= x < width:
        raise ValueError("x must satisfy 0 <= x < field_width_um")
    if y < 0:
        raise ValueError("y must be nonnegative")
    return width, x, y


def get_domain_id(x: float, y: float, K: int, field_width_um: float) -> int:
    """Return the equal-width vertical strip containing a physical position.

    Parameters
    ----------
    x, y : float
        Center coordinates in micrometers. Both must be finite and nonnegative;
        x must be strictly less than field_width_um. No upper y bound is checked
        because this interface has no field height. Only x determines the label.
    K : int
        Positive number of strips, numbered from 0 at the left to K-1.
    field_width_um : float
        Positive, finite tissue width in micrometers.

    Returns
    -------
    int
        Domain ID. An internal boundary belongs to its right-hand strip.

    Raises
    ------
    TypeError
        If an argument has an inappropriate type.
    ValueError
        If geometry or coordinates are invalid or x is outside the field.

    Notes
    -----
    This implements floor(x*K/W), capped at K-1. Floating boundaries are
    represented by W*i/K (multiplication before division). Binary search compares
    x directly with these boundaries, avoiding division rounding a point just
    below a boundary into the next strip. No epsilon snaps nearby points.
    Cell radii and measurement-square extents never influence this label.
    """
    width, x, y = _validate_position(x, y, K, field_width_um)

    # Search implicit boundaries rather than allocating a stored domain map.
    left, right = 0, int(K)
    while right - left > 1:
        middle = (left + right) // 2
        boundary = width * middle / int(K)
        if not np.isfinite(boundary):
            # Avoid intermediate overflow for very large finite field widths.
            boundary = width * (middle / int(K))
        if x < boundary:
            right = middle
        else:
            left = middle
    return left


def _position_domains(
    positions_df: pd.DataFrame,
    K: int,
    field_width_um: float,
    x_column: str,
    y_column: str,
) -> np.ndarray:
    """Validate coordinates and return one domain per position."""
    get_domain_id(0.0, 0.0, K, field_width_um)

    if not isinstance(positions_df, pd.DataFrame):
        raise TypeError("positions must be a pandas DataFrame")
    if not positions_df.columns.is_unique:
        raise ValueError("DataFrame column names must be unique")
    if not {x_column, y_column}.issubset(positions_df.columns):
        raise ValueError(
            f"positions must contain {x_column} and {y_column} columns"
        )

    labels = []
    for row, (x, y) in enumerate(
        zip(positions_df[x_column], positions_df[y_column], strict=True)
    ):
        try:
            labels.append(get_domain_id(x, y, K, field_width_um))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid coordinates at row position {row}: {exc}"
            ) from exc

    return np.asarray(labels, dtype=np.int64)


def _validate_sampling_inputs(
    cell_types: Sequence[str],
    domain_probabilities: ArrayLike,
    rng: np.random.Generator,
    K: int,
) -> tuple[list[str], np.ndarray]:
    """Validate the RNG, type names, and probabilities for categorical sampling."""
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator")
    # Column meanings require an explicit order; mappings are not name sequences.
    if isinstance(cell_types, (str, bytes, Mapping, AbstractSet)):
        raise TypeError(
            "cell_types must be an ordered sequence of unique nonblank strings"
        )
    try:
        names = list(cell_types)
    except TypeError as exc:
        raise ValueError("cell_types must be a sequence of strings") from exc
    if not names or any(not isinstance(n, str) or not n.strip() for n in names):
        raise ValueError("cell_types must contain nonblank strings")
    if len(set(names)) != len(names):
        raise ValueError("cell_types must be unique")
    try:
        raw_probabilities = np.asarray(domain_probabilities)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "domain_probabilities must be a numeric rectangular matrix"
        ) from exc
    # Check before float conversion, which could discard imaginary components.
    if np.iscomplexobj(raw_probabilities):
        raise ValueError("domain_probabilities must be real-valued")
    try:
        probabilities = np.array(raw_probabilities, dtype=float, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "domain_probabilities must be a numeric rectangular matrix"
        ) from exc
    if probabilities.shape != (K, len(names)):
        raise ValueError("domain_probabilities must have shape (K, len(cell_types))")
    if not np.isfinite(probabilities).all() or (probabilities < 0).any():
        raise ValueError("domain_probabilities must be finite and nonnegative")
    sums = probabilities.sum(axis=1)
    if not np.allclose(sums, 1.0, atol=1e-8, rtol=0.0):
        raise ValueError("domain_probabilities must sum to 1 in every domain")
    probabilities /= sums[:, None]

    return names, probabilities


def assign_cell_types(
    cells_df: pd.DataFrame,
    K: int,
    field_width_um: float,
    cell_types: Sequence[str],
    domain_probabilities: ArrayLike,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Copy positioned cells, assign center domains, and sample one type per cell.

    Parameters
    ----------
    cells_df : pandas.DataFrame
        One row per cell.
        Cell centres use center_x_um and center_y_um.
        Prepare cell-placement data with cells_from_placement() first.
        Extra columns and the index are preserved.
        Existing domain or cell_type columns are rejected.
    K : int
        Positive number of equal-width vertical strips.
    field_width_um : float
        Positive, finite tissue width in micrometers.
    cell_types : sequence of str
        Nonempty, unique, nonblank names in an explicit order (for example, a
        list or tuple). Sets and mappings are rejected. Order defines probability
        columns.
    domain_probabilities : array-like of shape (K, len(cell_types))
        Real, finite, nonnegative probabilities. Complex values are rejected.
        Row sums must equal 1 within absolute
        tolerance 1e-8 (no relative tolerance). Accepted rounding error is
        normalized on a private copy before sampling.
    rng : numpy.random.Generator
        Shared generator. One choice is made per cell in input row order;
        no generator or seed is created internally.

    Returns
    -------
    pandas.DataFrame
        Copy with integer domain and string cell_type columns.
        Cell coordinates remain center_x_um and center_y_um.
        Input data are unchanged.
        Cell radii do not affect assignment.


    Raises
    ------
    TypeError
        If an argument has an inappropriate type.
    ValueError
        For invalid geometry, coordinates, names, probabilities, RNG, or columns.
        All validation completes before any random draws are made.
    """
    domains = _position_domains(cells_df, K, field_width_um, "center_x_um", "center_y_um")
    if {"domain", "cell_type"}.intersection(cells_df.columns):
        raise ValueError("cells_df already contains domain or cell_type columns")
    names, probabilities = _validate_sampling_inputs(
        cell_types, domain_probabilities, rng, K
    )

    result = cells_df.copy(deep=True)
    result["domain"] = domains
    result["cell_type"] = pd.array(
        [rng.choice(names, p=probabilities[d]) for d in domains], dtype="string"
    )
    return result


def _square_domains(
    coordinates_df: pd.DataFrame, K: int, field_width_um: float
) -> np.ndarray:
    """Validate square IDs and return domains from square centres."""    
    domains = _position_domains(coordinates_df, K, field_width_um, "x_um", "y_um")

    if "square_id" not in coordinates_df:
        raise ValueError("coordinates_df must contain square_id")

    ids = coordinates_df["square_id"]
    if ids.isna().any() or any(
        not isinstance(s, str) or not s.strip() for s in ids
    ):
        raise ValueError("square_id values must be nonmissing, nonblank strings")
    if not ids.is_unique:
        raise ValueError("square_id values must be unique")

    return domains


def build_domains_dataframe(
    coordinates_df: pd.DataFrame, K: int, field_width_um: float
) -> pd.DataFrame:
    """Label measurement-square centers independently of cell occupancy.

    Parameters
    ----------
    coordinates_df : pandas.DataFrame
        Contains square_id, x_um and y_um. IDs must be unique nonblank strings.
        Coordinates are centers in micrometers. Array indices are optional here:
        grid completeness and 2 micrometer spacing belong to the grid component.
    K : int
        Positive number of equal-width vertical strips.
    field_width_um : float
        Positive, finite field width in micrometers.

    Returns
    -------
    pandas.DataFrame
        Exactly square_id and domain, in input order with the original index.
        Domain labels are strings, as required by Data Contract v1.0. Empty
        inputs with the required columns yield an empty table with this schema.

    Raises
    ------
    TypeError
        If an argument has an inappropriate type.
    ValueError
        For invalid geometry, coordinates, missing columns, or invalid IDs.
        Input data are never modified.

    Notes
    -----
    CSV has no type metadata. Readers should explicitly request string dtypes
    for square_id and domain, even though labels look numeric.
    """
    domains = _square_domains(coordinates_df, K, field_width_um)
    result = coordinates_df[["square_id"]].copy(deep=True)
    result["domain"] = pd.array([str(d) for d in domains], dtype="string")
    return result

