"""
PHASE 2 - STEP 10: Convert Drug Files to PDBQT Format
=======================================================
AutoDock Vina needs drugs in PDBQT format.
Our drugs are currently in SDF format.
We use OpenBabel to convert them.

SDF   = Standard chemistry format (from ChEMBL)
PDBQT = AutoDock format (needed for docking)
"""

import subprocess      # To run OpenBabel from Python
import os             # For file operations
import glob           # For finding multiple files at once

print("="*70)
print("CONVERTING DRUG FILES: SDF → PDBQT")
print("="*70)

# -------------------------------------------------------
# SECTION 1: FIND ALL SDF FILES
# -------------------------------------------------------

# glob.glob finds all files matching a pattern
# *.sdf means "any file ending in .sdf"
sdf_files = glob.glob("./drugs/*.sdf")

print(f"\nFound {len(sdf_files)} SDF files to convert:")
for f in sdf_files:
    print(f"  → {os.path.basename(f)}")

# -------------------------------------------------------
# SECTION 2: CONVERT EACH FILE
# -------------------------------------------------------

print("\n" + "-"*50)
print("CONVERTING FILES...")
print("-"*50)

converted = []     # List of successfully converted files
failed = []        # List of files that failed

for sdf_path in sdf_files:

    # Create output path by replacing .sdf with .pdbqt
    pdbqt_path = sdf_path.replace(".sdf", ".pdbqt")

    # Get just the filename for display
    filename = os.path.basename(sdf_path)

    # Run OpenBabel conversion
    # obabel = OpenBabel command
    # sdf_path = input file
    # -O pdbqt_path = output file
    # --gen3d = generate 3D coordinates if missing
    # -p 7.4 = set pH to 7.4 (physiological)
    #          This determines protonation state
    #          (which atoms have H+ attached)
    result = subprocess.run([
        "obabel",
        sdf_path,              # Input SDF file
        "-O", pdbqt_path,      # Output PDBQT file
        "--gen3d",             # Generate 3D coordinates
        "-p", "7.4"            # Physiological pH
    ],
    capture_output=True,       # Capture output silently
    text=True
    )

    # Check if output file was created successfully
    if os.path.exists(pdbqt_path):
        # Get file size to verify it is not empty
        size = os.path.getsize(pdbqt_path)

        if size > 0:
            converted.append(pdbqt_path)
            print(f"  ✅ {filename} → converted ({size} bytes)")
        else:
            # File created but empty = conversion failed
            failed.append(sdf_path)
            print(f"  ❌ {filename} → empty file")
            os.remove(pdbqt_path)    # Remove empty file
    else:
        failed.append(sdf_path)
        print(f"  ❌ {filename} → conversion failed")

# -------------------------------------------------------
# SECTION 3: SUMMARY
# -------------------------------------------------------

print("\n" + "="*70)
print("CONVERSION SUMMARY:")
print(f"  Total SDF files:     {len(sdf_files)}")
print(f"  Successfully converted: {len(converted)}")
print(f"  Failed:              {len(failed)}")

if failed:
    print("\nFailed files:")
    for f in failed:
        print(f"  → {f}")

# -------------------------------------------------------
# SECTION 4: VERIFY PDBQT FILES
# -------------------------------------------------------

print("\n" + "-"*50)
print("VERIFYING PDBQT FILES")
print("-"*50)

# Check a few converted files to make sure they look right
pdbqt_files = glob.glob("./drugs/*.pdbqt")

print(f"\nTotal PDBQT files created: {len(pdbqt_files)}")

# Show first few lines of first file
if pdbqt_files:
    print(f"\nSample from {os.path.basename(pdbqt_files[0])}:")
    with open(pdbqt_files[0], "r") as f:
        lines = f.readlines()[:5]
        for line in lines:
            print(f"  {line.rstrip()}")

print("\n✅ STEP 10 COMPLETE!")
print(f"Ready to dock {len(converted)} drugs against ITPKB")
print("Next: Run AutoDock Vina docking")