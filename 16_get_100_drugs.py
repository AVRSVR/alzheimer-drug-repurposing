"""
STEP 16: Download 100 FDA Approved Drugs
We increase our library to 100 molecules to find 
better candidates and provide more data for ML.
"""
import requests
import pandas as pd
import os
import time

os.makedirs("./drugs", exist_ok=True)

print("="*70)
print("DOWNLOADING 100 FDA DRUGS FROM ChEMBL")
print("="*70)

base_url = "https://www.ebi.ac.uk/chembl/api/data"
# We ask for 150 because some will fail the Lipinski filter
drug_url = f"{base_url}/molecule?max_phase=4&molecule_type=Small+molecule&format=json&limit=150"

response = requests.get(drug_url)
molecules = response.json()["molecules"]

drug_list = []
for mol in molecules:
    name = mol.get("pref_name", "Unknown")
    cid = mol.get("molecule_chembl_id", "")
    props = mol.get("molecule_properties", {})
    if not props: continue
    
    try:
        mw = float(props.get("full_mwt", 999))
        logp = float(props.get("alogp", 99))
        hbd = int(props.get("hbd", 99))
        hba = int(props.get("hba", 99))
        
        # Lipinski Filter
        if mw < 500 and logp < 5 and hbd <= 5 and hba <= 10 and mw > 150:
            drug_list.append({"chembl_id": cid, "name": name, "mw": mw, "logp": logp})
    except: continue

# Keep exactly 100
drug_list = drug_list[:100]
print(f"Drugs passing filter: {len(drug_list)}")

# Download SDF and convert to PDBQT immediately
for drug in drug_list:
    cid = drug["chembl_id"]
    if os.path.exists(f"./drugs/{cid}.pdbqt"): continue # Skip if already done
    
    sdf_url = f"{base_url}/molecule/{cid}.sdf"
    res = requests.get(sdf_url)
    if res.status_code == 200:
        sdf_path = f"./drugs/{cid}.sdf"
        with open(sdf_path, "wb") as f: f.write(res.content)
        
        # Convert to PDBQT using OpenBabel
        subprocess_cmd = f"obabel {sdf_path} -O {sdf_path.replace('.sdf', '.pdbqt')} --gen3d -p 7.4"
        os.system(subprocess_cmd)
        print(f"  ✅ Downloaded & Converted: {drug['name']}")
    time.sleep(0.2)

pd.DataFrame(drug_list).to_csv("./results/drug_library_100.csv", index=False)
print("\n✅ STEP 16 COMPLETE!")