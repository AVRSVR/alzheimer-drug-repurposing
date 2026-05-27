"""
PHASE 3 - STEP 15: Prepare Hub Proteins for Docking
=====================================================
We prepare all 3 hub proteins for AutoDock Vina.
Same process as Step 8 but for 3 proteins at once.

For each protein:
1. Analyze structure
2. Find center coordinates
3. Convert PDB to PDBQT
4. Save docking parameters
"""

import subprocess      # To run OpenBabel
import os             # For file operations
import json           # For saving parameters
import numpy as np    # For math
from Bio.PDB import PDBParser    # For reading PDB files

os.makedirs("./structures", exist_ok=True)

print("="*70)
print("PREPARING HUB PROTEINS FOR DOCKING")
print("="*70)

# -------------------------------------------------------
# DEFINE OUR 3 HUB PROTEINS
# -------------------------------------------------------

# Each protein needs its PDB file path
hub_proteins = [
    {
        "name": "CDC42",
        "pdb_path": "./structures/CDC42_pdb.pdb"
    },
    {
        "name": "HSPA8",
        "pdb_path": "./structures/HSPA8_pdb.pdb"
    },
    {
        "name": "PGK1",
        "pdb_path": "./structures/PGK1_pdb.pdb"
    }
]

# Store all docking parameters
all_params = {}

# -------------------------------------------------------
# PROCESS EACH PROTEIN
# -------------------------------------------------------

for protein in hub_proteins:
    name = protein["name"]
    pdb_path = protein["pdb_path"]

    print(f"\n{'='*50}")
    print(f"Processing: {name}")
    print(f"{'='*50}")

    # Check file exists
    if not os.path.exists(pdb_path):
        print(f"  ❌ File not found: {pdb_path}")
        continue

    # ── ANALYZE STRUCTURE ──
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, pdb_path)

    # Count atoms and chains
    atom_count = 0
    residue_count = 0
    chain_ids = []

    for model in structure:
        for chain in model:
            chain_ids.append(chain.id)
            for residue in chain:
                residue_count += 1
                for atom in residue:
                    atom_count += 1

    print(f"  Chains: {chain_ids}")
    print(f"  Residues: {residue_count}")
    print(f"  Atoms: {atom_count}")

    # ── FIND CENTER COORDINATES ──
    # Collect all atom coordinates
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())

    coords = np.array(coords)

    # Calculate center
    center_x = float(np.mean(coords[:, 0]))
    center_y = float(np.mean(coords[:, 1]))
    center_z = float(np.mean(coords[:, 2]))

    print(f"  Center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")

    # ── CONVERT TO PDBQT ──
    pdbqt_path = pdb_path.replace(".pdb", ".pdbqt")

    result = subprocess.run([
        "obabel",
        pdb_path,
        "-O", pdbqt_path,
        "-xr",          # Receptor mode
        "-p", "7.4"     # Physiological pH
    ],
    capture_output=True,
    text=True
    )

    if os.path.exists(pdbqt_path):
        size_kb = os.path.getsize(pdbqt_path) / 1024
        print(f"  ✅ PDBQT created: {pdbqt_path}")
        print(f"     Size: {size_kb:.1f} KB")
    else:
        print(f"  ❌ PDBQT conversion failed")
        continue

    # ── SAVE DOCKING PARAMETERS ──
    params = {
        "name": name,
        "pdbqt_path": pdbqt_path,
        "center_x": round(center_x, 3),
        "center_y": round(center_y, 3),
        "center_z": round(center_z, 3),
        "size_x": 80,    # Fixed box size
        "size_y": 80,    # Same as ITPKB
        "size_z": 80     # Works within Vina limits
    }

    all_params[name] = params

# -------------------------------------------------------
# SAVE ALL PARAMETERS TO ONE FILE
# -------------------------------------------------------

params_path = "./structures/hub_docking_params.json"

with open(params_path, "w") as f:
    json.dump(all_params, f, indent=2)

print("\n" + "="*70)
print("SUMMARY:")
print(f"  Proteins prepared: {len(all_params)}")
for name, params in all_params.items():
    print(f"\n  {name}:")
    print(f"    PDBQT: {params['pdbqt_path']}")
    print(f"    Center: ({params['center_x']}, "
          f"{params['center_y']}, {params['center_z']})")
    print(f"    Box: 80 x 80 x 80 Å")

print(f"\nAll parameters saved: {params_path}")
print("\n✅ STEP 15 COMPLETE!")
print("Next: Download 100 drugs and run expanded docking")
print("="*70)