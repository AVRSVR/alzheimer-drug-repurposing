"""
PHASE 2 - STEP 11: Run AutoDock Vina Docking
=============================================
We dock all 20 FDA approved drugs against ITPKB.
For each drug, Vina finds the best binding pose
and reports the binding energy in kcal/mol.

More negative = stronger binding = better drug candidate
"""

import subprocess      # To run vina.exe from Python
import os             # For file operations
import glob           # For finding multiple files
import json           # For loading docking parameters
import pandas as pd   # For saving results as table

print("="*70)
print("RUNNING AUTODOCK VINA DOCKING")
print("="*70)

# -------------------------------------------------------
# SECTION 1: LOAD DOCKING PARAMETERS
# -------------------------------------------------------

# Load the box coordinates we calculated in Step 8
# These tell Vina WHERE to search on the protein
with open("./structures/docking_params.json", "r") as f:
    params = json.load(f)

center_x = params["center_x"]
center_y = params["center_y"]
center_z = params["center_z"]
size_x   = params["size_x"]
size_y   = params["size_y"]
size_z   = params["size_z"]

print(f"\nDocking box center: ({center_x}, {center_y}, {center_z})")
print(f"Docking box size:   ({size_x}, {size_y}, {size_z})")

# -------------------------------------------------------
# SECTION 2: SET UP FILE PATHS
# -------------------------------------------------------

# Receptor = our prepared protein
receptor = "./structures/ITPKB_prepared.pdbqt"

# Get all drug PDBQT files
drug_files = glob.glob("./drugs/*.pdbqt")

print(f"\nReceptor: {receptor}")
print(f"Number of drugs to dock: {len(drug_files)}")

# Verify receptor exists
if not os.path.exists(receptor):
    print("❌ Receptor file not found!")
    exit()

# Create output folder for docking results
os.makedirs("./docking_results", exist_ok=True)

# -------------------------------------------------------
# SECTION 3: RUN DOCKING FOR EACH DRUG
# -------------------------------------------------------

print("\n" + "-"*50)
print("STARTING DOCKING RUNS...")
print("-"*50)
print("This may take 10-20 minutes total.")
print("Each drug takes about 1 minute.\n")

all_results = []    # Store results for all drugs

for i, drug_path in enumerate(drug_files):

    # Get drug name from filename
    drug_name = os.path.basename(drug_path).replace(".pdbqt", "")

    print(f"[{i+1}/{len(drug_files)}] Docking: {drug_name}")

    # Output file for this docking run
    output_path = f"./docking_results/{drug_name}_docked.pdbqt"
    log_path    = f"./docking_results/{drug_name}_log.txt"

    # Run AutoDock Vina
    # We call vina.exe directly using subprocess
    result = subprocess.run([
        "./vina.exe",              # Vina executable in our project folder
        "--receptor", receptor,    # Protein file
        "--ligand",   drug_path,   # Drug file
        "--out",      output_path, # Where to save docked poses
        "--log",      log_path,    # Where to save log
        "--center_x", str(center_x),   # Box center X
        "--center_y", str(center_y),   # Box center Y
        "--center_z", str(center_z),   # Box center Z
        "--size_x",   str(size_x),     # Box size X
        "--size_y",   str(size_y),     # Box size Y
        "--size_z",   str(size_z),     # Box size Z
        "--exhaustiveness", "8",   # Search thoroughness
                                   # Higher = more accurate but slower
                                   # 8 is standard for virtual screening
        "--num_modes", "5"         # Number of binding poses to generate
    ],
    capture_output=True,
    text=True
    )

    # -------------------------------------------------------
    # PARSE THE DOCKING RESULTS
    # -------------------------------------------------------

    # Vina writes results to the log file
    # We need to read and extract the best binding energy

    best_affinity = None

    if os.path.exists(log_path):
        with open(log_path, "r") as f:
            log_content = f.readlines()

        # The results table in Vina log looks like:
        # mode | affinity | dist from best mode
        #      | (kcal/mol)| rmsd l.b.| rmsd u.b.
        # -----+----------+----------+----------
        #    1      -7.5       0.000      0.000
        #    2      -7.1       1.234      2.456

        # We want the first mode (best binding energy)
        for line in log_content:
            line = line.strip()
            parts = line.split()

            # Look for lines starting with a number (mode number)
            if (len(parts) >= 4 and
                    parts[0].isdigit() and
                    parts[0] == "1"):        # First mode = best
                try:
                    best_affinity = float(parts[1])
                    break
                except ValueError:
                    continue

    # Store result
    if best_affinity is not None:
        all_results.append({
            "drug": drug_name,
            "best_affinity_kcal_mol": best_affinity
        })
        print(f"   ✅ Best affinity: {best_affinity} kcal/mol")
    else:
        all_results.append({
            "drug": drug_name,
            "best_affinity_kcal_mol": None
        })
        print(f"   ❌ Could not parse result")

# -------------------------------------------------------
# SECTION 4: SAVE AND RANK RESULTS
# -------------------------------------------------------

print("\n" + "="*70)
print("DOCKING COMPLETE - RANKING RESULTS")
print("="*70)

# Convert to DataFrame
results_df = pd.DataFrame(all_results)

# Remove failed runs
results_df = results_df.dropna(subset=["best_affinity_kcal_mol"])

# Sort by binding affinity
# Most negative = strongest binding = best drug
results_df = results_df.sort_values(
    "best_affinity_kcal_mol",
    ascending=True              # Most negative first
)

# Reset index after sorting
results_df = results_df.reset_index(drop=True)
results_df.index = results_df.index + 1    # Start ranking from 1

print("\nRANKED DRUG CANDIDATES:")
print("-"*50)
print(results_df.to_string())

# Save results
results_df.to_csv("./results/docking_results.csv")
print(f"\nResults saved: results/docking_results.csv")

print("\n" + "="*70)
print("✅ STEP 11 COMPLETE!")
print(f"Top drug candidate: {results_df.iloc[0]['drug']}")
print(f"Best binding energy: {results_df.iloc[0]['best_affinity_kcal_mol']} kcal/mol")
print("="*70)