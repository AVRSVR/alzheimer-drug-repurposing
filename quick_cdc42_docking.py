# quick_cdc42_docking.py
import subprocess, glob, os

drugs = glob.glob("./drugs/*.pdbqt")
receptor = "./structures/CDC42_clean.pdbqt"

for drug in drugs:
    name = os.path.basename(drug).replace(".pdbqt", "")
    out = f"./docking_results/CDC42_{name}_docked.pdbqt"
    log = f"./docking_results/CDC42_{name}_log.txt"
    
    if os.path.exists(log): 
        continue
        
    print(f"CDC42 vs {name}")
    subprocess.run([
        "./vina.exe", "--receptor", receptor, "--ligand", drug,
        "--center_x", "16.161", "--center_y", "-1.96", "--center_z", "15.694",
        "--size_x", "80", "--size_y", "80", "--size_z", "80",
        "--exhaustiveness", "4", "--out", out, "--log", log
    ])
print("CDC42 docking done!")