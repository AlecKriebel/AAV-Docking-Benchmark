"""Compare derived tables while retaining strict scientific identities/counts."""

import csv
import json
import math
from pathlib import Path

from .provenance import sha256


def compare_tables(expected_path, actual_path, numeric_tolerance=0.0):
    """Byte exact by default; explicit tolerance permits only derived float noise.

    Input/pose hashes are independently checked before table comparison. Boolean
    decisions, integer counts, atom/case identities and collection shapes must
    match exactly even when a cross-platform tolerance is requested.
    """
    if not math.isfinite(numeric_tolerance) or not 0 <= numeric_tolerance <= 1e-8:
        raise ValueError("Numeric tolerance must be finite and between 0 and 1e-8")
    expected_path, actual_path = Path(expected_path), Path(actual_path)
    expected_sha, actual_sha = sha256(expected_path), sha256(actual_path)
    if expected_sha == actual_sha:
        return {"exact_recomputation": True, "max_numeric_difference": 0.0, "recomputed_sha256": actual_sha}
    if numeric_tolerance == 0:
        raise ValueError(f"Recomputed table bytes differ: {expected_path.name}")
    differences = []

    def compare(a, b):
        if type(a) is not type(b):
            raise ValueError("Derived table types differ")
        if isinstance(a, float):
            if not math.isfinite(a) or not math.isfinite(b) or abs(a - b) > numeric_tolerance:
                raise ValueError("Derived floating value differs beyond tolerance")
            differences.append(abs(a - b))
        elif isinstance(a, dict):
            if a.keys() != b.keys():
                raise ValueError("Derived table fields differ")
            for key in a:
                compare(a[key], b[key])
        elif isinstance(a, list):
            if len(a) != len(b):
                raise ValueError("Derived table row counts differ")
            for left, right in zip(a, b):
                compare(left, right)
        elif a != b:
            raise ValueError("Derived table identity, count or decision differs")

    if expected_path.suffix == ".json":
        compare(json.loads(expected_path.read_text()), json.loads(actual_path.read_text()))
    elif expected_path.suffix == ".csv":
        def rows(path):
            with path.open(newline="") as stream:
                data = list(csv.reader(stream))
            if not data:
                raise ValueError("Empty derived CSV")
            result = [data[0]]
            for row in data[1:]:
                converted = []
                for value in row:
                    try:
                        converted.append(int(value))
                    except ValueError:
                        try:
                            converted.append(float(value))
                        except ValueError:
                            converted.append(value)
                result.append(converted)
            return result
        compare(rows(expected_path), rows(actual_path))
    else:
        raise ValueError("Only JSON and CSV derived tables support numeric comparison")
    return {"exact_recomputation": False, "max_numeric_difference": max(differences, default=0.0),
            "recomputed_sha256": actual_sha}
