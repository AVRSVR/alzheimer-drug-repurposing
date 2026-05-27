"""
Quick GRIN1 Docking
Docking all 100 drugs against GRIN1
GRIN1 = Glutamate NMDA receptor
Existing AD drug Memantine targets this protein
"""
import subprocess
import os
import glob

receptor = "./structures/GRIN1_clean.pdbqt"
drugs = glob.glob("./drugs/*.pdbqt")

print(f"Docking {len(drugs)} drugs against GRIN1...")
print("This takes ~30-40 minutes...")

for i, drug in enumerate(drugs):
    name = os.path.basename(drug).replace(".pdbqt", "")
    log = f"./docking_results/GRIN1_{name}_log.txt"
    out = f"./docking_results/GRIN1_{name}_docked.pdbqt"

    # Skip if already done
    if os.path.exists(log):
        continue

    print(f"[{i+1}/{len(drugs)}] GRIN1 vs {name}")

    subprocess.run([
        "./vina.exe",
        "--receptor", receptor,
        "--ligand", drug,
        "--out", out,
        "--log", log,
        "--center_x", "-27.787",
        "--center_y", "32.504",
        "--center_z", "-20.379",
        "--size_x", "80",
        "--size_y", "80",
        "--size_z", "80",
        "--exhaustiveness", "4",
        "--num_modes", "5"
    ], capture_output=True)

print("\n✅ GRIN1 DOCKING COMPLETE!")