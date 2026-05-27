"""
STEP 17: Expanded Multi-Target Docking
100 Drugs x 4 Proteins = 400 Docking Runs.
This is the main "Screening" phase of the project.
"""
import subprocess
import os
import glob
import json
import pandas as pd

# 1. Setup paths
os.makedirs("./docking_results", exist_ok=True)
drugs = glob.glob("./drugs/*.pdbqt")
with open("./structures/hub_docking_params.json", "r") as f: hub_params = json.load(f)
with open("./structures/docking_params.json", "r") as f: itpkb_params = json.load(f)

# Combine all 4 targets
targets = list(hub_params.values())
targets.append({
    "name": "ITPKB", 
    "pdbqt_path": "./structures/ITPKB_prepared.pdbqt",
    "center_x": itpkb_params["center_x"], "center_y": itpkb_params["center_y"], "center_z": itpkb_params["center_z"],
    "size_x": 80, "size_y": 80, "size_z": 80
})

all_results = []

print(f"Starting screen: {len(drugs)} drugs vs {len(targets)} targets (400 runs total)")

for target in targets:
    t_name = target["name"]
    for drug_path in drugs:
        d_name = os.path.basename(drug_path).replace(".pdbqt", "")
        out = f"./docking_results/{t_name}_{d_name}_docked.pdbqt"
        log = f"./docking_results/{t_name}_{d_name}_log.txt"
        
        if os.path.exists(log): continue # Skip if already done
        
        print(f"Docking {d_name} vs {t_name}...")
        subprocess.run([
            "./vina.exe", "--receptor", target["pdbqt_path"], "--ligand", drug_path,
            "--out", out, "--log", log,
            "--center_x", str(target["center_x"]), "--center_y", str(target["center_y"]), "--center_z", str(target["center_z"]),
            "--size_x", "80", "--size_y", "80", "--size_z", "80",
            "--exhaustiveness", "4" # Faster for screening
        ], capture_output=True)
        
        # Parse Energy
        try:
            with open(log, "r") as f:
                for line in f:
                    if line.strip().startswith("1 "): 
                        energy = float(line.split()[1])
                        all_results.append({"drug": d_name, "target": t_name, "energy": energy})
                        break
        except: continue

pd.DataFrame(all_results).to_csv("./results/expanded_docking_results.csv", index=False)
print("\n✅ ALL DOCKING COMPLETE!")