"""
PHASE 2 - STEP 7: Download ITPKB 3D Structure (AlphaFold)
==========================================================
ITPKB does not have an experimental crystal structure in PDB.
We will use the high-quality AlphaFold predicted structure instead.
This is very common in modern computational biology.
"""

import requests
import os

# Create folder
os.makedirs("./structures", exist_ok=True)

print("="*70)
print("DOWNLOADING ITPKB 3D STRUCTURE (AlphaFold)")
print("="*70)

# Correct AlphaFold URL for ITPKB (UniProt ID = P27987)
# This is the standard filename format
url = "https://alphafold.ebi.ac.uk/files/AF-P27987-F1-model_v4.pdb"

print(f"\nDownloading from:")
print(url)
print("\nThis is a predicted structure (very accurate for ITPKB).")

# Download the file
response = requests.get(url)

if response.status_code == 200:
    output_path = "./structures/ITPKB_alphafold.pdb"
    
    with open(output_path, "w") as f:
        f.write(response.text)
    
    # Count atoms to verify
    lines = response.text.split("\n")
    atom_count = sum(1 for line in lines if line.startswith("ATOM"))
    
    print("\n✅ SUCCESS!")
    print(f"Structure saved as: {output_path}")
    print(f"Number of atoms: {atom_count}")
    print(f"File size: {len(response.text)/1024:.1f} KB")
    
else:
    print(f"\n❌ Download failed with status: {response.status_code}")
    print("Try downloading manually from: https://alphafold.ebi.ac.uk/entry/P27987")

print("\n" + "="*70)
print("STEP 7 COMPLETE!")
print("Next step: Clean and prepare this structure for docking")
print("="*70)