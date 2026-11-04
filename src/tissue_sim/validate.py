"""Validate simulation output using the Week 7 notebook's layered checks.

Based on week07_notebook.ipynb: presence, counts, matching IDs, composition
sums. Extended to the full Data Contract 1.0 checklist, including reference
and aggregate ground truth. Public checks return silently on success.
"""

import argparse
import csv
import datetime
import json
import math
from pathlib import Path
import re

import numpy as np

REQUIRED_FILES = (
    "metadata.json",
    "counts.csv",
    "coordinates.csv",
    "reference/reference_counts.csv",
    "reference/reference_labels.csv",
    "ground_truth/params.json",
    "ground_truth/composition.csv",
    "ground_truth/domains.csv",
)


def _fail(path, detail):
    raise ValueError(f"{path}: {detail}")


def _json(path):
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        _fail(path, f"required readable JSON object; {exc}")
    if not isinstance(value, dict):
        _fail(path, f"expected JSON object, got {value!r}")
    return value


def _read_table(path, required):
    try:
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, strict=True)
            fields = reader.fieldnames
            if not fields:
                _fail(path, "empty file; required CSV header is missing")
            if len(set(fields)) != len(fields) or any(not v.strip() for v in fields):
                _fail(path, f"column names must be unique and nonblank, got {fields!r}")
            missing = set(required) - set(fields)
            if missing:
                _fail(
                    path,
                    f"missing required columns {sorted(missing)!r}; found {fields!r}",
                )
            rows = list(reader)
    except (OSError, csv.Error, UnicodeError) as exc:
        _fail(path, f"required readable CSV file; {exc}")
    for line, row in enumerate(rows, 2):
        if None in row or any(value is None for value in row.values()):
            _fail(path, f"row {line}: expected exactly {len(fields)} fields")
    return fields, rows


def read_csv(path):
    """Return a CSV file as a list of dictionaries, as in the Week 7 notebook.

    path names a UTF-8 CSV with a unique nonblank header. Values remain strings;
    no units or provenance are inferred here. Raises ValueError for malformed CSV.
    """
    _, rows = _read_table(Path(path), [])
    return rows


def _unique_ids(rows, key, path):
    seen = set()
    for row in rows:
        value = row[key]
        if not value.strip() or value in seen:
            _fail(path, f"{key}={value!r}; IDs must be unique and nonblank")
        seen.add(value)


def _integer(value, path, detail, positive=False):
    if not re.fullmatch(r"[0-9]+", value) or (positive and int(value) == 0):
        _fail(
            path,
            f"{detail}={value!r}; expected {'positive' if positive else 'non-negative'} integer",
        )
    return int(value)


def _read_counts(dataset):
    path = Path(dataset) / "counts.csv"
    fields, rows = _read_table(path, ["square_id"])
    if fields[0] != "square_id" or len(fields) < 2:
        _fail(
            path,
            f"columns={fields!r}; expected square_id first and at least one gene column",
        )
    _unique_ids(rows, "square_id", path)
    for line, row in enumerate(rows, 2):
        for gene in fields[1:]:
            row[gene] = _integer(
                row[gene], path, f"row {line}, {row['square_id']!r}, {gene}"
            )
    return fields[1:], rows


def check_counts(dataset):
    """Check counts.csv schema and non-negative integer transcript counts.

    dataset is a directory; rows are simulated measurement squares and columns
    are captured transcripts per gene, not cell counts or inferred fractions.
    Return None silently, or raise ValueError naming the file, square and value.
    """
    _read_counts(dataset)


def _same_order(actual, expected, path, key):
    if actual != expected:
        mismatch = next(
            (
                i
                for i, pair in enumerate(zip(actual, expected, strict=False))
                if pair[0] != pair[1]
            ),
            min(len(actual), len(expected)),
        )
        _fail(
            path,
            f"{key} mismatch at row {mismatch + 2}: got {actual[mismatch : mismatch + 1]!r}, expected {expected[mismatch : mismatch + 1]!r}; IDs and order must match",
        )


def check_ids_match(dataset):
    """Check that counts.csv and coordinates.csv list the same squares in order.

    dataset contains the simulated measurement tables; square_id identifies a
    square, independent of physical units. Return None, or raise ValueError
    naming the differing IDs, including an extra ID at the end of either file.
    """
    root = Path(dataset)
    _, counts = _read_table(root / "counts.csv", ["square_id"])
    _, coords = _read_table(root / "coordinates.csv", ["square_id"])
    _same_order(
        [r["square_id"] for r in counts],
        [r["square_id"] for r in coords],
        root / "counts.csv",
        "square_id",
    )


def check_composition_sums(dataset):
    """Check exact transcript conservation between composition and gene counts.

    dataset contains true per-square/type transcript counts, not fractions or
    inferred estimates. Missing sparse pairs mean zero. Return None silently,
    or raise ValueError naming the square, composition sum and gene-count total.
    This notebook-level check permits explicit zero rows; full validation also
    checks the contract's sparse-only-positive-row representation.
    """
    root = Path(dataset)
    genes, counts = _read_counts(root)
    expected = {row["square_id"]: sum(row[g] for g in genes) for row in counts}
    path = root / "ground_truth/composition.csv"
    _, composition = _read_table(path, ["square_id", "cell_type", "n_transcripts"])
    totals = dict.fromkeys(expected, 0)
    for row in composition:
        square = row["square_id"]
        if square not in totals:
            _fail(
                path,
                f"square_id={square!r}; must occur in coordinates.csv and counts.csv",
            )
        totals[square] += _integer(
            row["n_transcripts"], path, f"square_id={square!r}, n_transcripts"
        )
    for square, total in totals.items():
        if total != expected[square]:
            _fail(
                path,
                f"square_id={square!r}, n_transcripts sum={total}; must equal counts.csv total={expected[square]} (missing pairs mean zero)",
            )


def _metadata(meta, path):
    if meta.get("contract_version") != "1.0":
        _fail(
            path,
            f"contract_version={meta.get('contract_version')!r}; supported version is '1.0'",
        )
    for key in ("n_squares", "n_genes", "random_seed"):
        if type(meta.get(key)) is not int or meta[key] < (
            0 if key == "random_seed" else 1
        ):
            _fail(
                path,
                f"{key}={meta.get(key)!r}; expected {'non-negative' if key == 'random_seed' else 'positive'} integer",
            )
    grid = meta.get("grid")
    if not isinstance(grid, dict) or any(
        type(grid.get(k)) is not int or grid[k] <= 0 for k in ("n_rows", "n_cols")
    ):
        _fail(path, f"grid={grid!r}; n_rows and n_cols must be positive integers")
    if grid["n_rows"] * grid["n_cols"] != meta["n_squares"]:
        _fail(
            path, f"n_squares={meta['n_squares']!r}; must equal grid rows times columns"
        )
    if (
        type(meta.get("square_size_um")) not in (int, float)
        or meta["square_size_um"] != 2
    ):
        _fail(
            path,
            f"square_size_um={meta.get('square_size_um')!r}; contract 1.0 requires 2 micrometres",
        )
    field = meta.get("field_um")
    for key, size in (("width", grid["n_cols"]), ("height", grid["n_rows"])):
        if (
            not isinstance(field, dict)
            or type(field.get(key)) not in (int, float)
            or field[key] != size * meta["square_size_um"]
        ):
            _fail(
                path,
                f"field_um={field!r}; {key} must agree with grid and square_size_um",
            )
    types = meta.get("cell_types")
    if (
        not isinstance(types, list)
        or not types
        or any(not isinstance(t, str) or not t.strip() for t in types)
        or len(set(types)) != len(types)
    ):
        _fail(path, f"cell_types={types!r}; expected unique nonblank type names")
    generator = meta.get("generator")
    if not isinstance(generator, dict) or any(
        not isinstance(generator.get(k), str) or not generator[k].strip()
        for k in ("name", "version", "git_commit")
    ):
        _fail(
            path,
            f"generator={generator!r}; name, version and git_commit must be nonblank strings",
        )
    try:
        timestamp = datetime.datetime.fromisoformat(
            meta.get("created_utc", "").replace("Z", "+00:00")
        )
        if timestamp.utcoffset() != datetime.timedelta(0):
            raise ValueError("timestamp must be UTC")
    except (ValueError, TypeError, AttributeError) as exc:
        _fail(
            path,
            f"created_utc={meta.get('created_utc')!r}; expected ISO-8601 UTC timestamp ({exc})",
        )


def _params(params, meta, path):
    required = {
        "seed",
        "cell_types",
        "expression_programs",
        "marker_genes",
        "n_markers_per_type",
        "marker_gene_boost",
        "reference_n_per_type",
    }
    if required - params.keys():
        _fail(path, f"missing required keys {sorted(required - params.keys())!r}")
    if (
        type(params["seed"]) is not int
        or params["seed"] != meta["random_seed"]
        or params["cell_types"] != meta["cell_types"]
    ):
        _fail(
            path,
            "seed/cell_types must match metadata.json random_seed/cell_types, including order",
        )
    if any("bin_size" in key for key in params):
        _fail(
            path,
            "bin_size is an analysis parameter and must not appear in simulation params",
        )
    for key in ("n_markers_per_type", "reference_n_per_type"):
        if type(params[key]) is not int or params[key] < (
            1 if key == "reference_n_per_type" else 0
        ):
            _fail(
                path,
                f"{key}={params[key]!r}; expected {'positive' if key == 'reference_n_per_type' else 'non-negative'} integer",
            )
    boost = params["marker_gene_boost"]
    if type(boost) not in (int, float) or not math.isfinite(boost) or boost < 1:
        _fail(path, f"marker_gene_boost={boost!r}; expected finite fold-change >= 1")
    try:
        programs = np.array(params["expression_programs"], dtype=float)
    except (ValueError, TypeError):
        _fail(path, "expression_programs must contain numeric probability vectors")
    if (
        programs.shape != (len(meta["cell_types"]), meta["n_genes"])
        or not np.isfinite(programs).all()
        or (programs < 0).any()
        or not np.allclose(programs.sum(axis=1), 1, atol=1e-8, rtol=0)
    ):
        _fail(
            path,
            "expression_programs must have one finite non-negative n_genes vector per type, each summing to 1",
        )
    markers = params["marker_genes"]
    if not isinstance(markers, dict) or set(markers) != set(meta["cell_types"]):
        _fail(path, f"marker_genes={markers!r}; keys must match cell_types")
    for cell_type, indices in markers.items():
        if (
            not isinstance(indices, list)
            or len(indices) != params["n_markers_per_type"]
            or any(type(i) is not int or not 0 <= i < meta["n_genes"] for i in indices)
            or len(set(indices)) != len(indices)
        ):
            _fail(
                path,
                f"marker_genes[{cell_type!r}]={indices!r}; expected n_markers_per_type unique zero-based gene indices",
            )


def _check_geometry(root, meta):
    path = root / "coordinates.csv"
    _, rows = _read_table(path, ["square_id", "array_row", "array_col", "x_um", "y_um"])
    _unique_ids(rows, "square_id", path)
    positions = []
    for row in rows:
        r = _integer(
            row["array_row"], path, f"square_id={row['square_id']!r}, array_row"
        )
        c = _integer(
            row["array_col"], path, f"square_id={row['square_id']!r}, array_col"
        )
        positions.append((r, c))
    expected_positions = {
        (r, c)
        for r in range(meta["grid"]["n_rows"])
        for c in range(meta["grid"]["n_cols"])
    }
    if (
        len(positions) != meta["n_squares"]
        or len(set(positions)) != len(positions)
        or set(positions) != expected_positions
    ):
        _fail(
            path,
            f"grid has {len(positions)} rows and {len(set(positions))} unique positions; expected full grid of {meta['n_squares']} positions without gaps or duplicates",
        )
    for row, (r, c) in zip(rows, positions, strict=True):
        for column, index in (("x_um", c), ("y_um", r)):
            expected = (index + 0.5) * meta["square_size_um"]
            try:
                value = float(row[column])
            except ValueError:
                value = math.nan
            if not math.isfinite(value) or not math.isclose(
                value, expected, abs_tol=1e-8, rel_tol=0
            ):
                _fail(
                    path,
                    f"square_id={row['square_id']!r}, {column}={row[column]!r}; expected {expected} micrometres from grid index and square_size_um",
                )
    return [row["square_id"] for row in rows]


def _check_reference(root, meta, genes, params):
    path = root / "reference/reference_counts.csv"
    fields, counts = _read_table(path, ["cell_id"])
    if fields != ["cell_id", *genes]:
        _fail(
            path, "gene columns must match counts.csv in the same order, after cell_id"
        )
    _unique_ids(counts, "cell_id", path)
    for row in counts:
        for gene in genes:
            _integer(row[gene], path, f"cell_id={row['cell_id']!r}, {gene}")
    path = root / "reference/reference_labels.csv"
    _, labels = _read_table(path, ["cell_id", "cell_type"])
    _same_order(
        [r["cell_id"] for r in labels], [r["cell_id"] for r in counts], path, "cell_id"
    )
    types = [r["cell_type"] for r in labels]
    if set(types) != set(meta["cell_types"]):
        _fail(
            path,
            f"cell_type values={sorted(set(types))!r}; reference must cover metadata cell_types={meta['cell_types']!r}",
        )
    for cell_type in meta["cell_types"]:
        n = types.count(cell_type)
        if n != params["reference_n_per_type"]:
            _fail(
                path,
                f"cell_type={cell_type!r}, reference cells={n}; expected reference_n_per_type={params['reference_n_per_type']}",
            )


def _check_ground_truth(root, meta, square_ids):
    path = root / "ground_truth/domains.csv"
    _, domains = _read_table(path, ["square_id", "domain"])
    _same_order([r["square_id"] for r in domains], square_ids, path, "square_id")
    for row in domains:
        if not row["domain"].strip():
            _fail(
                path,
                f"square_id={row['square_id']!r}, domain=''; every square requires a nonblank label",
            )
    path = root / "ground_truth/composition.csv"
    _, composition = _read_table(path, ["square_id", "cell_type", "n_transcripts"])
    ranks = {square: i for i, square in enumerate(square_ids)}
    previous_rank = -1
    seen = set()
    for row in composition:
        square, cell_type = row["square_id"], row["cell_type"]
        _integer(
            row["n_transcripts"],
            path,
            f"square_id={square!r}, n_transcripts",
            positive=True,
        )
        if square not in ranks:
            _fail(path, f"square_id={square!r}; must occur in coordinates.csv")
        if cell_type not in meta["cell_types"]:
            _fail(
                path, f"cell_type={cell_type!r}; must occur in metadata.json cell_types"
            )
        if (square, cell_type) in seen:
            _fail(
                path,
                f"duplicate (square_id, cell_type) pair={(square, cell_type)!r}; aggregate contributions into one row",
            )
        seen.add((square, cell_type))
        if ranks[square] < previous_rank:
            _fail(
                path,
                f"square_id={square!r}; square_id order must follow coordinates.csv, omitting zero pairs",
            )
        previous_rank = ranks[square]


def validate_dataset(path):
    """Validate a complete simulation dataset, following the Week 7 notebook.

    path is a directory of simulated transcript counts and aggregate ground
    truth, with coordinates in micrometres. Check required files, then counts,
    matching IDs and exact composition sums; extend those notebook checks with
    full metadata, geometry, reference and ground-truth conformance checks.
    Return None silently on success; raise ValueError on the first violation.
    The notebook's partial teaching fixture is not a full contract dataset.
    """
    dataset = Path(path)
    for name in REQUIRED_FILES:
        if not (dataset / name).is_file():
            _fail(
                dataset / name,
                f"required readable {'JSON' if name.endswith('.json') else 'CSV'} file is missing (contract section 2)",
            )
    check_counts(dataset)
    check_ids_match(dataset)
    # Validate representation before totals so duplicate pairs are reported clearly.
    meta = _json(dataset / "metadata.json")
    _metadata(meta, dataset / "metadata.json")
    square_ids = _check_geometry(dataset, meta)
    genes, _ = _read_counts(dataset)
    if len(genes) != meta["n_genes"]:
        _fail(
            dataset / "counts.csv",
            f"gene columns={len(genes)}; expected n_genes={meta['n_genes']}",
        )
    _check_ground_truth(dataset, meta, square_ids)
    check_composition_sums(dataset)
    params = _json(dataset / "ground_truth/params.json")
    _params(params, meta, dataset / "ground_truth/params.json")
    _check_reference(dataset, meta, genes, params)


def main():
    """Validate one dataset directory; successful validation prints nothing."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    try:
        validate_dataset(args.path)
    except ValueError as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
