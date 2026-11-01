"""Controls for exact scientific identities and offline site replay inputs."""
import gzip
import hashlib
import json

import pytest

from aav_benchmark.site_reproduction import compare_site_values, restore_truth_sources, _check_sources


def test_geometry_tolerance_does_not_change_discrete_decisions_or_copied_values():
    original = {"center_to_native_atom_A": 3.9999999999999,
                "top1_success_4A": True, "score": -8.5, "rank": 1,
                "center_A": [1.2, 3.4, 5.6], "coordinates": [[1.0, 2.0, 3.0]],
                "residue_keys": [["A", 14, "", "GAL"]]}
    rounded = {**original, "center_to_native_atom_A": 4.0}
    check = compare_site_values(original, rounded, numeric_tolerance=1e-10)
    assert check["exact_recomputation"] is False
    assert 0 < check["max_numeric_difference"] < 1e-10
    for field, changed in (("top1_success_4A", False), ("score", -8.5 + 1e-12),
                           ("rank", 2), ("center_A", [1.2 + 1e-12, 3.4, 5.6]),
                           ("coordinates", [[1.0 + 1e-12, 2.0, 3.0]]),
                           ("residue_keys", [["B", 14, "", "GAL"]])):
        with pytest.raises(ValueError):
            compare_site_values(original, {**original, field: changed}, numeric_tolerance=1e-10)


def test_tolerance_inherits_matrix_field_and_rejects_types_shapes_and_nonfinite_values():
    original = {"rotation": [[1.0, 0.0], [0.0, 1.0]], "count": 3}
    rounded = {"rotation": [[1.0, 1e-14], [0.0, 1.0]], "count": 3}
    assert compare_site_values(original, rounded, numeric_tolerance=1e-10)["max_numeric_difference"] == 1e-14
    with pytest.raises(ValueError):
        compare_site_values(original, rounded)
    for changed in ({**original, "count": 3.0}, {**original, "rotation": [[1.0, 0.0]]},
                    {**original, "rotation": [[float("nan"), 0.0], [0.0, 1.0]]}):
        with pytest.raises(ValueError):
            compare_site_values(original, changed, numeric_tolerance=1e-10)
    for tolerance in (-1.0, float("inf"), True, 1e-5):
        with pytest.raises(ValueError):
            compare_site_values(original, original, numeric_tolerance=tolerance)


def _offline_fixture(root):
    benchmarks = root / "benchmarks"
    benchmarks.mkdir()
    archive = root / "data/archive"
    archive.mkdir(parents=True)
    sources = []
    for filename, payload in (("7MTZ.cif", b"native capsid source"),
                              ("TEST_r_b.pdb", b"native protein"),
                              ("TEST_l_b.pdb", b"native ligand")):
        compressed = archive / (filename + ".gz")
        compressed.write_bytes(gzip.compress(payload, mtime=0))
        sources.append({"filename": filename, "sha256": hashlib.sha256(payload).hexdigest(),
                        "archive": compressed.relative_to(root).as_posix()})
    (benchmarks / "site_experiment.json").write_text(json.dumps({"cases": [
        {"condition_id": "aav9_holo", "complex_id": "aav9"},
        {"condition_id": "aav9_apo", "complex_id": "aav9"},
        {"condition_id": "TEST", "complex_id": "TEST"}]}))
    (benchmarks / "sources.lock.json").write_text(json.dumps({"sources": sources[:1]}))
    (benchmarks / "external_cohort.json").write_text(json.dumps({"cases": [
        {"case_id": "TEST", "sources": sources[1:]}]}))
    return sources


def test_native_truth_is_restored_offline_and_never_from_existing_raw(tmp_path):
    sources = _offline_fixture(tmp_path)
    existing_raw = tmp_path / "data/raw"
    existing_raw.mkdir()
    (existing_raw / "7MTZ.cif").write_bytes(b"untrusted unrelated local file")
    restored = tmp_path / "restored"
    receipt = restore_truth_sources(tmp_path, restored)
    assert len(receipt) == 3  # Paired apo is not a second native truth source.
    assert (restored / "7MTZ.cif").read_bytes() == b"native capsid source"
    assert (restored / "external/TEST_l_b.pdb").read_bytes() == b"native ligand"
    assert receipt[0]["sha256"] == sources[0]["sha256"]
    with pytest.raises(FileExistsError):
        restore_truth_sources(tmp_path, restored)


def test_native_truth_refuses_missing_or_corrupt_archive_without_fallback(tmp_path):
    sources = _offline_fixture(tmp_path)
    (tmp_path / sources[0]["archive"]).write_bytes(gzip.compress(b"mutated native source"))
    with pytest.raises(ValueError, match="checksum"):
        restore_truth_sources(tmp_path, tmp_path / "corrupt")
    (tmp_path / sources[0]["archive"]).unlink()
    with pytest.raises(ValueError, match="Missing offline"):
        restore_truth_sources(tmp_path, tmp_path / "missing")


def test_source_inventory_is_verified_before_any_numeric_tolerance(tmp_path):
    source = tmp_path / "scores.json"
    source.write_bytes(b'{"score":-8.5}')
    inventory = {"scores.json": hashlib.sha256(source.read_bytes()).hexdigest()}
    _check_sources(tmp_path, inventory)
    source.write_bytes(b'{"score":-8.500000000001}')
    with pytest.raises(ValueError, match="checksum"):
        _check_sources(tmp_path, inventory)
    with pytest.raises(ValueError):
        _check_sources(tmp_path, {"../scores.json": inventory["scores.json"]})
