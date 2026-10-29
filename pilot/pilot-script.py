"""AAV9 local, reference-site-informed pilot. Not blind docking.
Run from project root: python src/pilot.py prepare|dock|analyze
"""

import argparse
import json
import pathlib
import time
import platform
import importlib.metadata
import csv
import numpy as np
import gemmi
from scipy.spatial.distance import cdist

# Paths assume this script is stored in the project's src directory.
PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data/processed"
RESULTS_DIR = PROJECT_ROOT / "results"
PROCESSED_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)


def xyz(residue):
    """Return the non-hydrogen atom coordinates of a residue as a NumPy array."""
    return np.array(
        [list(atom.pos) for atom in residue if not atom.element.is_hydrogen]
    )


def protein(structure):
    """Return protein-containing chains from the first structure model."""
    return {
        chain.name: chain
        for chain in structure[0]
        if any((residue.het_flag == "A" for residue in chain))
    }


def atomdict(chains):
    """Index protein heavy-atom coordinates by chain, residue number/name, and atom name."""
    return {
        (chain_name, residue.seqid.num, residue.name, atom.name): np.array(
            list(atom.pos)
        )
        for chain_name, chain in chains.items()
        for residue in chain
        if residue.het_flag == "A"
        for atom in residue
        if not atom.element.is_hydrogen
    }


def fit(source_coords, target_coords):
    """Find a rigid rotation and translation that align source coordinates to target coordinates."""
    source_center = source_coords.mean(0)
    target_center = target_coords.mean(0)
    u, _, vt = np.linalg.svd(
        (source_coords - source_center).T @ (target_coords - target_center)
    )
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return (rotation, target_center - source_center @ rotation)


# Preparation: receptor alignment, ligand generation, and reference metadata.
def prepare():
    """Align apo/holo receptors, prepare galactose, and save reference data and metadata."""
    holo_structure = gemmi.read_structure(str(PROJECT_ROOT / "data/raw/7MTZ.cif"))
    apo_structure = gemmi.read_structure(str(PROJECT_ROOT / "data/raw/7MT0.cif"))
    ligand = next(
        (
            residue
            for chain in holo_structure[0]
            for residue in chain
            if residue.name == "GAL"
        )
    )
    reference_coords = xyz(ligand)
    center = reference_coords.mean(0)
    holo_chains = protein(holo_structure)
    apo_chains = protein(apo_structure)
    holo_atoms = atomdict(holo_chains)
    apo_atoms = atomdict(apo_chains)
    # Entire subunits touching a 20-A ligand shell; retains neighbors around a 20-A box.
    selected = [
        chain_name
        for chain_name, chain in holo_chains.items()
        if cdist(
            np.vstack([xyz(residue) for residue in chain if residue.het_flag == "A"]),
            reference_coords,
        ).min()
        < 20
    ]
    keys = [
        atom_key
        for atom_key in holo_atoms
        if atom_key[0] in selected and atom_key[3] == "CA" and (atom_key in apo_atoms)
    ]
    rotation, translation = fit(
        np.array([apo_atoms[atom_key] for atom_key in keys]),
        np.array([holo_atoms[atom_key] for atom_key in keys]),
    )
    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (
                        np.array([apo_atoms[atom_key] for atom_key in keys]) @ rotation
                        + translation
                        - np.array([holo_atoms[atom_key] for atom_key in keys])
                    )
                    ** 2,
                    axis=1,
                )
            )
        )
    )
    # Preserve author numbering, remap selected chains to short IDs for legacy PDB tools.
    chainmap = {chain_name: chr(65 + i) for i, chain_name in enumerate(selected)}
    for name, protein_chains in [("holo", holo_chains), ("apo", apo_chains)]:
        output_structure = gemmi.Structure()
        output_structure.add_model(gemmi.Model("1"))
        for chain_name in selected:
            chain = protein_chains[chain_name].clone()
            chain.name = chainmap[chain_name]
            for i in reversed(range(len(chain))):
                if chain[i].het_flag != "A":
                    del chain[i]
            if name == "apo":
                for res in chain:
                    for atom in res:
                        position = np.array(list(atom.pos)) @ rotation + translation
                        atom.pos = gemmi.Position(*position)
            output_structure[0].add_chain(chain)
        output_structure.write_pdb(str(PROCESSED_DIR / f"{name}_receptor.pdb"))
    contacts = []
    for atom_key, position in holo_atoms.items():
        distances = np.linalg.norm(reference_coords - position, axis=1)
        for i in np.where(distances < 4.0)[0]:
            contacts.append(
                {
                    "chain": atom_key[0],
                    "resnum": atom_key[1],
                    "resname": atom_key[2],
                    "atom": atom_key[3],
                    "ligand_atom": ligand[i].name,
                    "distance_A": float(distances[i]),
                }
            )
    # Build ligand from CCD explicit connectivity; never infer bonds from pose proximity.
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    ccd_block = gemmi.cif.read_file(str(PROJECT_ROOT / "data/raw/GAL.cif")).sole_block()
    editable_molecule = Chem.RWMol()
    atom_indices = {}
    for row in ccd_block.find("_chem_comp_atom.", ["atom_id", "type_symbol", "charge"]):
        atom = Chem.Atom(row[1])
        atom.SetFormalCharge(int(row[2]))
        atom.SetProp("ccd_name", row[0])
        atom_indices[row[0]] = editable_molecule.AddAtom(atom)
    for row in ccd_block.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        editable_molecule.AddBond(
            atom_indices[row[0]],
            atom_indices[row[1]],
            {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE}[row[2]],
        )
    mol = editable_molecule.GetMol()
    Chem.SanitizeMol(mol)
    # Set ideal CCD coordinates (including H), assign stereochemistry there, then embed independently.
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for row in ccd_block.find(
        "_chem_comp_atom.",
        [
            "atom_id",
            "pdbx_model_Cartn_x_ideal",
            "pdbx_model_Cartn_y_ideal",
            "pdbx_model_Cartn_z_ideal",
        ],
    ):
        conformer.SetAtomPosition(
            atom_indices[row[0]], tuple(map(float, list(row)[1:]))
        )
    mol.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    initial_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    mol.RemoveAllConformers()
    embedding_parameters = AllChem.ETKDGv3()
    embedding_parameters.randomSeed = 20261024
    if AllChem.EmbedMolecule(mol, embedding_parameters) != 0:
        raise RuntimeError("ligand embedding failed")
    AllChem.MMFFOptimizeMolecule(mol)
    conformer = mol.GetConformer()
    positions = conformer.GetPositions()
    positions += center - positions.mean(0)
    for i, position in enumerate(positions):
        conformer.SetAtomPosition(i, position)
    # Names survive Meeko PDBQT output via PDB residue info.
    for atom in mol.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(atom.GetProp("ccd_name").rjust(4))
        info.SetResidueName("GAL")
        info.SetResidueNumber(1)
        atom.SetMonomerInfo(info)
    setup = MoleculePreparation().prepare(mol)[0]
    pdbqt_text, success, error_message = PDBQTWriterLegacy.write_string(setup)
    if not success:
        raise RuntimeError(error_message)
    (PROCESSED_DIR / "galactose.pdbqt").write_text(pdbqt_text)
    sdf_writer = Chem.SDWriter(str(PROCESSED_DIR / "galactose.sdf"))
    sdf_writer.write(mol)
    sdf_writer.close()
    (PROCESSED_DIR / "reference.json").write_text(
        json.dumps(
            {
                "names": [atom.name for atom in ligand],
                "coordinates": reference_coords.tolist(),
            },
            indent=2,
        )
    )
    result = {
        "status": "structural_preparation_complete",
        "selected_author_chains": selected,
        "pdb_chain_mapping": chainmap,
        "box_center": center.tolist(),
        "box_size": [20, 20, 20],
        "selection_shell_A": 20,
        "apo_to_holo_CA_rmsd_A": rmsd,
        "alignment_CA_count": len(keys),
        "rotation": rotation.tolist(),
        "translation": translation.tolist(),
        "galactose_isomeric_smiles": initial_smiles,
        "contacts_4A": contacts,
        "contact_cutoff_A": 4.0,
        "reference_ligand": "7MTZ author chain A GAL 801",
        "information_regime": "true-site-local; not blind",
        "versions": {
            x: importlib.metadata.version(x)
            for x in ["gemmi", "rdkit", "vina", "meeko", "numpy", "scipy"]
        },
        "platform": platform.platform(),
    }
    (RESULTS_DIR / "preparation.json").write_text(json.dumps(result, indent=2))
    print(
        json.dumps(
            {
                atom_key: v
                for atom_key, v in result.items()
                if atom_key != "contacts_4A"
            },
            indent=2,
        )
    )


# Docking: repeated local searches for both receptor structures.
def dock():
    """Run Vina for each prepared receptor and seed, recording pose files and run metadata."""
    from vina import Vina

    preparation_metadata = json.loads((RESULTS_DIR / "preparation.json").read_text())
    runs = []
    for receptor in ["holo", "apo"]:
        if not (PROCESSED_DIR / f"{receptor}_receptor.pdbqt").exists():
            (RESULTS_DIR / f"{receptor}_blocked.json").write_text(
                json.dumps(
                    {
                        "status": "blocked",
                        "reason": "prepared receptor PDBQT missing; inspect preparation log",
                    }
                )
            )
            continue
        for seed in [104729, 130363, 155921]:
            start = time.perf_counter()
            docking_engine = Vina(sf_name="vina", cpu=2, seed=seed, verbosity=0)
            docking_engine.set_receptor(
                str(PROCESSED_DIR / f"{receptor}_receptor.pdbqt")
            )
            docking_engine.set_ligand_from_file(str(PROCESSED_DIR / "galactose.pdbqt"))
            docking_engine.compute_vina_maps(
                center=preparation_metadata["box_center"],
                box_size=preparation_metadata["box_size"],
            )
            docking_engine.dock(exhaustiveness=16, n_poses=20, min_rmsd=0.5)
            path = RESULTS_DIR / f"{receptor}_{seed}.pdbqt"
            docking_engine.write_poses(
                str(path), n_poses=20, energy_range=10, overwrite=True
            )
            runs.append(
                {
                    "receptor": receptor,
                    "seed": seed,
                    "wall_seconds": time.perf_counter() - start,
                    "exhaustiveness": 16,
                    "cpu": 2,
                    "requested_poses": 20,
                    "pose_file": path.name,
                }
            )
            (RESULTS_DIR / "runs.json").write_text(json.dumps(runs, indent=2))
            print(runs[-1], flush=True)


# Analysis: atom-matched RMSD and recovery of reference contacts.
def analyze():
    """Compare docked poses with the reference and save RMSD and contact metrics."""
    reference = json.loads((PROCESSED_DIR / "reference.json").read_text())
    reference = dict(zip(reference["names"], reference["coordinates"]))
    names = list(reference)
    target = np.array([reference[atom_name] for atom_name in names])
    rows = []
    for run in json.loads((RESULTS_DIR / "runs.json").read_text()):
        receptor = gemmi.read_structure(
            str(PROCESSED_DIR / f"{run['receptor']}_receptor.pdb")
        )
        receptor_atoms = atomdict(protein(receptor))
        receptor_atom_keys = list(receptor_atoms)
        receptor_coords = np.array(list(receptor_atoms.values()))
        holo = atomdict(
            protein(gemmi.read_structure(str(PROCESSED_DIR / "holo_receptor.pdb")))
        )
        holo_coords = np.array([holo[k] for k in receptor_atom_keys])
        native_contacts = set(map(tuple, np.argwhere(cdist(target, holo_coords) < 4)))
        for rank, block in enumerate(
            (RESULTS_DIR / run["pose_file"]).read_text().split("MODEL")[1:], 1
        ):
            atomlines = [
                line
                for line in block.splitlines()
                if line.startswith(("ATOM", "HETATM"))
            ]
            pose_coords = {
                line[12:16].strip(): [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
                for line in atomlines
                if line[12:16].strip() in reference
            }
            if set(pose_coords) != set(reference):
                raise ValueError("Atom identity mismatch; refusing RMSD")
            pose = np.array([pose_coords[atom_name] for atom_name in names])
            rmsd = float(np.sqrt(np.mean(np.sum((pose - target) ** 2, axis=1))))
            predicted_contacts = set(
                map(tuple, np.argwhere(cdist(pose, receptor_coords) < 4))
            )
            # Each receptor compared with the same reference ligand in its aligned frame.
            score = float(
                next(
                    (
                        line
                        for line in block.splitlines()
                        if line.startswith("REMARK VINA RESULT:")
                    )
                ).split()[3]
            )
            rows.append(
                {
                    **run,
                    "rank": rank,
                    "vina_score": score,
                    "rmsd_A": rmsd,
                    "contact_recall": (
                        len(native_contacts & predicted_contacts) / len(native_contacts)
                        if native_contacts
                        else None
                    ),
                    "contact_precision": (
                        len(native_contacts & predicted_contacts)
                        / len(predicted_contacts)
                        if predicted_contacts
                        else 0
                    ),
                    "success_2A": rmsd <= 2,
                    "success_3A": rmsd <= 3,
                }
            )
    with (RESULTS_DIR / "pose_metrics.csv").open("w") as output_file:
        csv_writer = csv.DictWriter(output_file, fieldnames=rows[0])
        csv_writer.writeheader()
        csv_writer.writerows(rows)
    summary = []
    for run in json.loads((RESULTS_DIR / "runs.json").read_text()):
        run_rows = [pose for pose in rows if pose["pose_file"] == run["pose_file"]]
        best = min(run_rows, key=lambda pose: pose["rmsd_A"])
        summary.append(
            {
                **run,
                "returned_poses": len(run_rows),
                "top1_rmsd_A": run_rows[0]["rmsd_A"],
                "top5_min_rmsd_A": min((pose["rmsd_A"] for pose in run_rows[:5])),
                "best_returned_rmsd_A": best["rmsd_A"],
                "best_returned_rank": best["rank"],
                "top1_score": run_rows[0]["vina_score"],
            }
        )
    (RESULTS_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["prepare", "dock", "analyze"])
    globals()[p.parse_args().stage]()
