#!/usr/bin/env python3
"""Replay all site metrics and 900 protein-only surface draws from offline archives."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.provenance import write_json
from aav_benchmark.site_reproduction import verify_site_experiment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--numeric-tolerance", type=float, default=0.0,
                        help="explicit tolerance for named derived geometry fields; default exact")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        verification = verify_site_experiment(args.root, numeric_tolerance=args.numeric_tolerance)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(1, f"Site experiment replay failed: {error}\n")
    if args.output is not None:
        write_json(args.output, verification)
    maximum = max(row["max_numeric_difference"] for row in verification["tables"])
    print(f"Verified {len(verification['tables'])} site tables and "
          f"{sum(row['replicates_resampled'] for row in verification['surface_nulls'])} surface draws; "
          f"all table values exact: {verification['all_tables_exact']}; "
          f"maximum derived difference: {maximum:.3g}")


if __name__ == "__main__":
    main()
