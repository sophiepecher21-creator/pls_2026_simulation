# Week 8 mutation record

Command: `python scripts/check_week08_mutations.py`.

Scope: the seven cases in `tests/test_simulation_invariants.py`, run against temporary copies of the current source. This is a focused mutation experiment, not a full-suite run or an exhaustive mutation analysis.

Baseline: PASS (pytest exit code 0).

| Mutation | Result | Failing tests |
| --- | --- | --- |
| Subtract instead of accumulating transcript counts | CAUGHT | `tests/test_simulation_invariants.py::test_aggregation_preserves_exact_gene_and_type_counts`; `tests/test_simulation_invariants.py::test_four_8um_bins_sum_to_each_16um_bin_gene_by_gene` |
| Reverse the square-to-matrix-row mapping | CAUGHT | `tests/test_simulation_invariants.py::test_aggregation_preserves_exact_gene_and_type_counts`; `tests/test_simulation_invariants.py::test_four_8um_bins_sum_to_each_16um_bin_gene_by_gene` |
| Ignore the configured seed | CAUGHT | `tests/test_simulation_invariants.py::test_seed_controls_random_sampling` |

Every case starts from a fresh copy. The working source is never mutated. A caught mutant demonstrates protection against that specific change; it does not establish complete correctness.
