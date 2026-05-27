"""
PHASE 2 - STEP 8: Prepare Protein for Docking
===============================================
We convert our AlphaFold PDB file into PDBQT format.
PDBQT is what AutoDock Vina uses for docking.
We also find the center of the protein to define
the docking search space.
"""

import subprocess        # Lets Python run external programs (like obabel)
import os               # For file and folder operations
from Bio.PDB import PDBParser    # For reading and analyzing PDB files
import numpy as np      # For math calculations

# -------------------------------------------------------
# SECTION 1: VERIFY INPUT FILE EXISTS
# -------------------------------------------------------

input_pdb = "./structures/ITPKB_alphafold.pdb"
output_pdbqt = "./structures/ITPKB_prepared.pdbqt"

print("="*70)
print("PREPARING ITPKB STRUCTURE FOR DOCKING")
print("="*70)

# Check if our PDB file exists
if os.path.exists(input_pdb):
    print(f"\n✅ Found input file: {input_pdb}")
    
    # Get file size
    size_kb = os.path.getsize(input_pdb) / 1024
    print(f"   File size: {size_kb:.1f} KB")
else:
    print(f"\n❌ Input file not found: {input_pdb}")
    print("   Make sure ITPKB_alphafold.pdb is in structures folder")
    exit()

# -------------------------------------------------------
# SECTION 2: ANALYZE THE PDB STRUCTURE
# -------------------------------------------------------

print("\n" + "-"*50)
print("ANALYZING PROTEIN STRUCTURE")
print("-"*50)

# PDBParser reads and parses PDB files
# QUIET=True means don't print warnings
parser = PDBParser(QUIET=True)
structure = parser.get_structure("ITPKB", input_pdb)

# Count atoms and residues
atom_count = 0
residue_count = 0
chain_ids = []

# Loop through structure hierarchy:
# Structure → Model → Chain → Residue → Atom
for model in structure:
    for chain in model:
        chain_ids.append(chain.id)
        for residue in chain:
            residue_count += 1
            for atom in residue:
                atom_count += 1

print(f"Chains found: {chain_ids}")
print(f"Total residues (amino acids): {residue_count}")
print(f"Total atoms: {atom_count}")

# -------------------------------------------------------
# SECTION 3: FIND THE CENTER OF THE PROTEIN
# -------------------------------------------------------

print("\n" + "-"*50)
print("FINDING PROTEIN CENTER (for docking box)")
print("-"*50)

# We need to know where the protein is in 3D space
# The docking box must cover the whole protein
# (since we don't know exact binding site yet)

# Collect all atom coordinates (X, Y, Z positions)
coords = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                coords.append(atom.get_coord())

# Convert to numpy array for easy math
coords = np.array(coords)

# Calculate center of protein
# np.mean() calculates average of each coordinate
center_x = np.mean(coords[:, 0])    # Average X
center_y = np.mean(coords[:, 1])    # Average Y
center_z = np.mean(coords[:, 2])    # Average Z

# Calculate size of protein (max - min in each direction)
size_x = np.max(coords[:, 0]) - np.min(coords[:, 0])
size_y = np.max(coords[:, 1]) - np.min(coords[:, 1])
size_z = np.max(coords[:, 2]) - np.min(coords[:, 2])

# Add padding around protein for docking box
# We add 10 Angstroms on each side
padding = 10
box_size_x = size_x + padding
box_size_y = size_y + padding
box_size_z = size_z + padding

print(f"\nProtein center:")
print(f"  X: {center_x:.2f} Å")
print(f"  Y: {center_y:.2f} Å")
print(f"  Z: {center_z:.2f} Å")

print(f"\nDocking box size:")
print(f"  X: {box_size_x:.2f} Å")
print(f"  Y: {box_size_y:.2f} Å")
print(f"  Z: {box_size_z:.2f} Å")

# Save these coordinates for the docking step
# We will need them in Script 10
import json
docking_params = {
    "center_x": round(float(center_x), 3),
    "center_y": round(float(center_y), 3),
    "center_z": round(float(center_z), 3),
    "size_x": round(float(box_size_x), 3),
    "size_y": round(float(box_size_y), 3),
    "size_z": round(float(box_size_z), 3)
}

with open("./structures/docking_params.json", "w") as f:
    json.dump(docking_params, f, indent=2)

print("\nDocking parameters saved to: structures/docking_params.json")

# -------------------------------------------------------
# SECTION 4: CONVERT PDB TO PDBQT USING OPENBABEL
# -------------------------------------------------------

print("\n" + "-"*50)
print("CONVERTING PDB TO PDBQT FORMAT")
print("-"*50)

print("\nRunning OpenBabel conversion...")
print("This adds charges and atom types needed by Vina")

# subprocess.run() runs an external program from Python
# obabel is OpenBabel's command line tool
# -i pdb = input format is PDB
# -o pdbqt = output format is PDBQT
# -xr = receptor mode (for proteins)
# -p 7.4 = physiological pH (like in human body)

result = subprocess.run([
    "obabel",
    input_pdb,              # Input file
    "-O", output_pdbqt,     # Output file
    "-xr",                  # Receptor mode
    "-p", "7.4"             # Physiological pH
],
capture_output=True,        # Capture any output/errors
text=True                   # Return as text not bytes
)

# Check if conversion worked
if os.path.exists(output_pdbqt):
    size_kb = os.path.getsize(output_pdbqt) / 1024
    print(f"\n✅ Conversion successful!")
    print(f"   Output: {output_pdbqt}")
    print(f"   File size: {size_kb:.1f} KB")
    
    # Count ATOM lines in PDBQT
    with open(output_pdbqt, "r") as f:
        pdbqt_lines = f.readlines()
    atom_lines = [l for l in pdbqt_lines if l.startswith("ATOM")]
    print(f"   Atom entries: {len(atom_lines)}")
    
else:
    print(f"\n❌ Conversion failed")
    print(f"   Error: {result.stderr}")

print("\n" + "="*70)
print("✅ STEP 8 COMPLETE!")
print("Protein is prepared and ready for docking")
print(f"Docking box center: ({center_x:.1f}, {center_y:.1f}, {center_z:.1f})")
print("="*70)