"""Recompute committed structural and assay summaries from archived input bytes.

Uses temporary working copies and never overwrites the deposited raw results.
Run from the repository root: python scripts/verify_archived_results.py
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.archive import restore
from aav_benchmark.engine import analyze
from aav_benchmark.panel import summarize_panel
from aav_benchmark.provenance import sha256, write_json


def verify(root):
    root = Path(root).resolve()
    original = root / "results/development"
    checks = []
    with tempfile.TemporaryDirectory(prefix="aav-archive-verification-") as directory:
        temporary = Path(directory)
        prepared_root, results_root = temporary / "prepared", temporary / "results"
        for output in sorted(original.iterdir()):
            if not output.is_dir() or not (output / "runs.json").is_file():
                continue
            prepared, scratch = prepared_root / output.name, results_root / output.name
            restore(output, prepared)
            scratch.mkdir(parents=True)
            runs = json.loads((output / "runs.json").read_text())
            files = {"runs.json", "configuration.json"} | {r["pose_file"] for r in runs if r["status"] == "ok"}
            for filename in files:
                if Path(filename).name != filename:
                    raise ValueError("Result filenames must not contain directories")
                shutil.copyfile(output / filename, scratch / filename)
            reference = json.loads((prepared / "reference.json").read_text())
            if reference.get("pose_reference_eligible") is False:
                continue
            analyze(prepared, scratch)
            for filename in ("summary.json", "pose_metrics.csv"):
                expected, actual = sha256(output / filename), sha256(scratch / filename)
                if expected != actual:
                    raise ValueError(f"Recomputed structural table differs: {output.name}/{filename}")
                checks.append({"file": str((output / filename).relative_to(root)), "sha256": actual, "exact_recomputation": True})
        summarize_panel(prepared_root, results_root, root / "benchmarks/galectin3_itc_panel.json")
        panel_name = "gal3_panel_summary.json"
        actual, expected = sha256(results_root / panel_name), sha256(original / panel_name)
        if actual != expected:
            raise ValueError("Recomputed assay ranking differs from committed table")
        checks.append({"file": "results/development/" + panel_name, "sha256": actual, "exact_recomputation": True})
    result = {"scope": "Exact recomputation from archived prepared inputs and raw poses; no docking rerun or independent structural-quality adjudication",
              "all_tables_match": True, "tables": checks, "verifier_sha256": sha256(Path(__file__))}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    verification = verify(args.root)
    write_json(args.root / "results/archive_verification.json", verification)
    print(f"Verified exact recomputation of {len(verification['tables'])} tables")
