"""
PHASE 4 - STEP 20: Validation Docking of Top Candidates
=========================================================
We download and dock our 6 validated candidates
against all 4 hub proteins.
exhaustiveness=16 for more accurate results.

Tier 1: PROCHLORPERAZINE, FLUPHENAZINE, PERPHENAZINE
Tier 2: ACETOPHENAZINE, DELAVIRDINE
Tier 3: LEVORPHANOL
"""

import requests
import subprocess
import os
import json
import pandas as pd
import time

os.makedirs("./validation", exist_ok=True)
os.makedirs("./validation/structures", exist_ok=True)

print("="*70)
print("VALIDATION DOCKING OF TOP ML CANDIDATES")
print("="*70)

# -------------------------------------------------------
# SECTION 1: DEFINE CANDIDATES
# -------------------------------------------------------

candidates = [
    {"chembl_id": "CHEMBL728",  "name": "PROCHLORPERAZINE", "tier": 1},
    {"chembl_id": "CHEMBL726",  "name": "FLUPHENAZINE",      "tier": 1},
    {"chembl_id": "CHEMBL567",  "name": "PERPHENAZINE",      "tier": 1},
    {"chembl_id": "CHEMBL1085", "name": "ACETOPHENAZINE",    "tier": 2},
    {"chembl_id": "CHEMBL593",  "name": "DELAVIRDINE",       "tier": 2},
    {"chembl_id": "CHEMBL592",  "name": "LEVORPHANOL",       "tier": 3},
]

# -------------------------------------------------------
# SECTION 2: DEFINE TARGET PROTEINS
# -------------------------------------------------------

# Load hub protein docking parameters
with open("./structures/hub_docking_params.json") as f:
    hub_params = json.load(f)

with open("./structures/docking_params.json") as f:
    itpkb_params = json.load(f)

# All 4 targets
targets = [
    {
        "name": "GRIN1",
        "pdbqt": "./structures/GRIN1_clean.pdbqt",
        "cx": -27.787, "cy": 32.504, "cz": -20.379
    },
    {
        "name": "HSPA8",
        "pdbqt": hub_params["HSPA8"]["pdbqt_path"],
        "cx": hub_params["HSPA8"]["center_x"],
        "cy": hub_params["HSPA8"]["center_y"],
        "cz": hub_params["HSPA8"]["center_z"]
    },
    {
        "name": "ITPKB",
        "pdbqt": "./structures/ITPKB_prepared.pdbqt",
        "cx": itpkb_params["center_x"],
        "cy": itpkb_params["center_y"],
        "cz": itpkb_params["center_z"]
    },
    {
        "name": "PGK1",
        "pdbqt": hub_params["PGK1"]["pdbqt_path"],
        "cx": hub_params["PGK1"]["center_x"],
        "cy": hub_params["PGK1"]["center_y"],
        "cz": hub_params["PGK1"]["center_z"]
    }
]

# -------------------------------------------------------
# SECTION 3: DOWNLOAD DRUG STRUCTURES
# -------------------------------------------------------

print("\n1. Downloading candidate structures...")

base_url = "https://www.ebi.ac.uk/chembl/api/data"

for candidate in candidates:
    cid = candidate["chembl_id"]
    name = candidate["name"]
    sdf_path = f"./validation/structures/{cid}.sdf"
    pdbqt_path = f"./validation/structures/{cid}.pdbqt"

    # Skip if already downloaded
    if os.path.exists(pdbqt_path):
        print(f"  ✅ {name} already exists")
        continue

    # Download SDF
    url = f"{base_url}/molecule/{cid}.sdf"
    response = requests.get(url)

    if response.status_code == 200:
        with open(sdf_path, "wb") as f:
            f.write(response.content)

        # Convert to PDBQT
        os.system(
            f"obabel {sdf_path} "
            f"-O {pdbqt_path} "
            f"--gen3d -p 7.4"
        )

        if os.path.exists(pdbqt_path):
            print(f"  ✅ {name} downloaded and converted")
        else:
            print(f"  ❌ {name} conversion failed")
    else:
        print(f"  ❌ {name} download failed")

    time.sleep(0.3)

# -------------------------------------------------------
# SECTION 4: RUN VALIDATION DOCKING
# -------------------------------------------------------

print("\n2. Running validation docking...")
print("Using exhaustiveness=16 (more accurate than screening)")
print("24 total runs (6 drugs × 4 proteins)")
print("Estimated time: ~45-60 minutes\n")

all_results = []

for candidate in candidates:
    cid = candidate["chembl_id"]
    name = candidate["name"]
    tier = candidate["tier"]
    pdbqt_path = f"./validation/structures/{cid}.pdbqt"

    if not os.path.exists(pdbqt_path):
        print(f"  ❌ Skipping {name} - no PDBQT file")
        continue

    print(f"\n[Tier {tier}] {name}:")

    for target in targets:
        t_name = target["name"]
        out = f"./validation/{name}_{t_name}_docked.pdbqt"
        log = f"./validation/{name}_{t_name}_log.txt"

        print(f"  Docking vs {t_name}...", end=" ")

        subprocess.run([
            "./vina.exe",
            "--receptor", target["pdbqt"],
            "--ligand",   pdbqt_path,
            "--out",      out,
            "--log",      log,
            "--center_x", str(target["cx"]),
            "--center_y", str(target["cy"]),
            "--center_z", str(target["cz"]),
            "--size_x",   "80",
            "--size_y",   "80",
            "--size_z",   "80",
            "--exhaustiveness", "16",
            "--num_modes", "9"
        ], capture_output=True)

        # Parse result
        energy = None
        if os.path.exists(log):
            with open(log) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2 and parts[0] == "1":
                        try:
                            energy = float(parts[1])
                            break
                        except:
                            continue

        if energy:
            print(f"✅ {energy} kcal/mol")
            all_results.append({
                "drug": name,
                "tier": tier,
                "target": t_name,
                "energy": energy
            })
        else:
            print(f"❌ failed")

# -------------------------------------------------------
# SECTION 5: ANALYZE AND SAVE RESULTS
# -------------------------------------------------------

print("\n" + "="*70)
print("VALIDATION RESULTS")
print("="*70)

results_df = pd.DataFrame(all_results)

# Pivot for easy reading
pivot = results_df.pivot(
    index='drug',
    columns='target',
    values='energy'
)

# Add average and tier
pivot['avg_energy'] = pivot.mean(axis=1)

# Add tier info
tier_map = {c['name']: c['tier'] for c in candidates}
pivot['tier'] = pivot.index.map(tier_map)

# Sort by average energy
pivot = pivot.sort_values('avg_energy')

print("\nValidation Docking Results:")
print("-"*70)
print(pivot.to_string())

# Save
pivot.to_csv("./results/validation_results.csv")
results_df.to_csv("./results/validation_detailed.csv", index=False)

print("\n" + "="*70)
print("FINAL DRUG RANKING:")
print("="*70)

for i, (drug, row) in enumerate(pivot.iterrows()):
    print(f"\nRank {i+1}: {drug} (Tier {int(row['tier'])})")
    print(f"  Average binding: {row['avg_energy']:.2f} kcal/mol")
    for target in ['GRIN1', 'HSPA8', 'ITPKB', 'PGK1']:
        if target in row:
            print(f"  {target}: {row[target]:.1f} kcal/mol")

print("\n✅ STEP 20 COMPLETE!")
print("Saved: results/validation_results.csv")
print("Next: Literature search + GitHub setup")