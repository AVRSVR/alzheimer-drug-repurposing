"""
PHASE 2 - STEP 9: Get FDA Approved Drug Library
=================================================
We download FDA approved drugs from ChEMBL.
We now use the CORRECT field names from the API.
"""

import requests
import pandas as pd
import os
import time

os.makedirs("./drugs", exist_ok=True)

print("="*70)
print("DOWNLOADING FDA APPROVED DRUG LIBRARY FROM ChEMBL")
print("="*70)

# -------------------------------------------------------
# SECTION 1: DOWNLOAD FDA APPROVED DRUGS
# -------------------------------------------------------

print("\nSearching ChEMBL for FDA approved small molecules...")

base_url = "https://www.ebi.ac.uk/chembl/api/data"

drug_url = (
    f"{base_url}/molecule?"
    f"max_phase=4&"
    f"molecule_type=Small+molecule&"
    f"format=json&"
    f"limit=100"           # Get 100 drugs
)

response = requests.get(drug_url)

if response.status_code == 200:
    data = response.json()
    molecules = data["molecules"]
    print(f"Found {len(molecules)} FDA approved small molecules")
else:
    print(f"Failed: {response.status_code}")
    molecules = []

# -------------------------------------------------------
# SECTION 2: FILTER USING LIPINSKI RULE OF 5
# -------------------------------------------------------

print("\nApplying Lipinski Rule of 5 filter...")
print("This keeps only drug-like molecules that")
print("can be absorbed by the human body.")

drug_list = []

for mol in molecules:

    # Get name from TOP level (not inside molecule_properties)
    pref_name = mol.get("pref_name", "Unknown")
    chembl_id = mol.get("molecule_chembl_id", "")

    # Get properties from inside molecule_properties
    props = mol.get("molecule_properties", {})

    if not props:
        continue

    # Extract values using CORRECT field names
    # We saw these exact names in the debug output
    try:
        mw   = float(props.get("full_mwt", 999))    # Molecular weight
        logp = float(props.get("alogp", 99))          # Lipophilicity
        hbd  = int(props.get("hbd", 99))              # H-bond donors
        hba  = int(props.get("hba", 99))              # H-bond acceptors
    except (TypeError, ValueError):
        continue

    # Lipinski Rule of 5:
    # MW < 500    → not too big to be absorbed
    # LogP < 5    → not too greasy
    # HBD <= 5    → not too many H-bond donors
    # HBA <= 10   → not too many H-bond acceptors
    # MW > 150    → not too small (not just a salt)

    if (mw < 500 and
            logp < 5 and
            hbd <= 5 and
            hba <= 10 and
            mw > 150):

        drug_list.append({
            "chembl_id": chembl_id,
            "name": pref_name,
            "molecular_weight": mw,
            "logp": logp,
            "hbd": hbd,
            "hba": hba
        })

print(f"\nDrugs passing Lipinski filter: {len(drug_list)}")

# -------------------------------------------------------
# SECTION 3: SHOW WHAT WE FOUND
# -------------------------------------------------------

drug_df = pd.DataFrame(drug_list)

print("\nDrug candidates:")
print(drug_df[["name", "molecular_weight", "logp"]].to_string())

# -------------------------------------------------------
# SECTION 4: DOWNLOAD SDF STRUCTURES
# -------------------------------------------------------

print("\n" + "-"*50)
print("DOWNLOADING DRUG 3D STRUCTURES (SDF format)")
print("-"*50)

downloaded = []
failed = []

for i, drug in enumerate(drug_list[:20]):   # Download first 20
    chembl_id = drug["chembl_id"]
    name = drug["name"]

    sdf_url = f"{base_url}/molecule/{chembl_id}.sdf"
    sdf_response = requests.get(sdf_url)

    if sdf_response.status_code == 200:
        output_path = f"./drugs/{chembl_id}.sdf"
        with open(output_path, "wb") as f:
            f.write(sdf_response.content)
        downloaded.append(chembl_id)
        print(f"  ✅ {name} ({chembl_id})")
    else:
        failed.append(chembl_id)
        print(f"  ❌ {name} ({chembl_id}) - failed")

    time.sleep(0.3)    # Small delay to not overload server

# -------------------------------------------------------
# SECTION 5: SAVE DRUG LIST
# -------------------------------------------------------

drug_df.to_csv("./drugs/drug_candidates.csv", index=False)

print("\n" + "="*70)
print("SUMMARY:")
print(f"  Total candidates:          {len(drug_list)}")
print(f"  Successfully downloaded:   {len(downloaded)}")
print(f"  Failed:                    {len(failed)}")
print(f"\nDrug list saved: drugs/drug_candidates.csv")
print("="*70)
print("\n✅ STEP 9 COMPLETE!")
print("Next: Convert SDF to PDBQT and run docking")