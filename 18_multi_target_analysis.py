"""
PHASE 4 - STEP 19: Proper ML with Morgan Fingerprints
=======================================================
We:
1. Re-download drug library WITH smiles strings
2. Generate Morgan fingerprints using RDKit
3. Train proper ML model (2048+ features)
4. Screen 500 new drugs using trained model
5. Find best multi-target candidates
"""

import pandas as pd
import numpy as np
import requests
import time
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, MACCSkeys
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import matplotlib.pyplot as plt
import os

os.makedirs("./results", exist_ok=True)
os.makedirs("./plots", exist_ok=True)

# -------------------------------------------------------
# SECTION 1: RE-DOWNLOAD DRUG LIBRARY WITH SMILES
# -------------------------------------------------------

print("="*70)
print("STEP 1: DOWNLOADING DRUG LIBRARY WITH SMILES")
print("="*70)

base_url = "https://www.ebi.ac.uk/chembl/api/data"

# Download 100 FDA approved drugs
drug_url = (
    f"{base_url}/molecule?"
    f"max_phase=4&"
    f"molecule_type=Small+molecule&"
    f"format=json&"
    f"limit=150"
)

response = requests.get(drug_url)
molecules = response.json()["molecules"]

drug_list = []
for mol in molecules:
    name = mol.get("pref_name", "Unknown")
    cid = mol.get("molecule_chembl_id", "")
    props = mol.get("molecule_properties", {})
    structs = mol.get("molecule_structures", {})

    if not props or not structs:
        continue

    # Get SMILES string
    # SMILES = text representation of molecular structure
    # Example: "CC(=O)Oc1ccccc1C(=O)O" = Aspirin
    smiles = structs.get("canonical_smiles", "")

    if not smiles:
        continue

    try:
        mw   = float(props.get("full_mwt", 999))
        logp = float(props.get("alogp", 99))
        hbd  = int(props.get("hbd", 99))
        hba  = int(props.get("hba", 99))

        # Lipinski filter
        if (mw < 500 and logp < 5 and
                hbd <= 5 and hba <= 10 and mw > 150):

            drug_list.append({
                "chembl_id": cid,
                "name": name,
                "smiles": smiles,
                "mw": mw,
                "logp": logp,
                "hbd": hbd,
                "hba": hba
            })
    except:
        continue

drug_df = pd.DataFrame(drug_list)
print(f"Drugs with SMILES: {len(drug_df)}")
print(f"Columns: {drug_df.columns.tolist()}")

# -------------------------------------------------------
# SECTION 2: GENERATE MOLECULAR FINGERPRINTS
# -------------------------------------------------------

print("\n" + "="*70)
print("STEP 2: GENERATING MOLECULAR FINGERPRINTS")
print("="*70)
print("Morgan fingerprints = 2048 bits encoding molecular structure")
print("MACCS keys = 167 bits encoding chemical features")

def get_fingerprints(smiles):
    """
    Convert SMILES string to molecular fingerprints.
    
    Morgan fingerprint (ECFP4):
    → Looks at each atom
    → Encodes what is around it within 2 bonds
    → Creates 2048 bit vector
    → Each bit = presence/absence of a substructure
    
    MACCS keys:
    → 167 predefined chemical questions
    → "Does molecule have a ring?"
    → "Does it have a carbonyl group?"
    → Each answer = 0 or 1
    """
    try:
        # Convert SMILES to RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            return None

        # Morgan fingerprint (radius=2 = ECFP4)
        # nBits=2048 = 2048 bit vector
        morgan = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=2, nBits=2048
        )

        # MACCS keys
        maccs = MACCSkeys.GenMACCSKeys(mol)

        # Combine both into one array
        fp = np.array(list(morgan) + list(maccs))

        return fp

    except:
        return None

# Generate fingerprints for all drugs
print("\nGenerating fingerprints...")
fingerprints = []
valid_drugs = []

for _, row in drug_df.iterrows():
    fp = get_fingerprints(row["smiles"])
    if fp is not None:
        fingerprints.append(fp)
        valid_drugs.append(row)

fp_array = np.array(fingerprints)
valid_df = pd.DataFrame(valid_drugs).reset_index(drop=True)

print(f"Valid drugs with fingerprints: {len(valid_df)}")
print(f"Fingerprint size per drug: {fp_array.shape[1]} features")
print(f"(2048 Morgan + 167 MACCS = 2215 bits)")

# -------------------------------------------------------
# SECTION 3: MERGE WITH DOCKING RESULTS
# -------------------------------------------------------

print("\n" + "="*70)
print("STEP 3: MERGING FINGERPRINTS WITH DOCKING RESULTS")
print("="*70)

# Load docking results
docking = pd.read_csv("./results/all_docking_results.csv")

# Pivot: rows=drugs, columns=targets
pivot = docking.pivot(
    index='drug',
    columns='target',
    values='energy'
).fillna(docking['energy'].mean())

print(f"Docking results: {pivot.shape}")

# Merge fingerprints with docking results
valid_df_indexed = valid_df.set_index('chembl_id')

# Find drugs present in BOTH fingerprint data AND docking results
common_drugs = [
    d for d in pivot.index
    if d in valid_df_indexed.index
]

print(f"Drugs in both datasets: {len(common_drugs)}")

# Build feature matrix and labels
X_list = []
y_list = {target: [] for target in pivot.columns}

for drug in common_drugs:
    # Get fingerprint index
    drug_idx = valid_df[valid_df['chembl_id'] == drug].index[0]
    X_list.append(fp_array[drug_idx])

    # Get binding energies for each target
    for target in pivot.columns:
        y_list[target].append(pivot.loc[drug, target])

X = np.array(X_list)
print(f"\nFeature matrix shape: {X.shape}")
print(f"(drugs × fingerprint features)")

# -------------------------------------------------------
# SECTION 4: TRAIN ML MODELS
# -------------------------------------------------------

print("\n" + "="*70)
print("STEP 4: TRAINING RANDOM FOREST MODELS")
print("="*70)
print("One model per target protein")
print("Input: 2215 molecular fingerprint features")
print("Output: Predicted binding energy (kcal/mol)")

models = {}
performance = {}

for target in pivot.columns:
    y = np.array(y_list[target])

    # Split into training and test sets
    # 80% training, 20% testing
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.2,
        random_state=42
    )

    print(f"\nTraining model for {target}...")
    print(f"  Training samples: {len(X_train)}")
    print(f"  Test samples: {len(X_test)}")

    # Random Forest Regressor
    model = RandomForestRegressor(
        n_estimators=100,    # 100 decision trees
        max_depth=10,        # Maximum tree depth
        random_state=42,
        n_jobs=-1            # Use all CPU cores
    )
    model.fit(X_train, y_train)

    # Evaluate on test set
    y_pred = model.predict(X_test)
    r2   = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    models[target] = model
    performance[target] = {
        'r2': r2,
        'rmse': rmse,
        'y_test': y_test,
        'y_pred': y_pred
    }

    print(f"  R² Score: {r2:.3f}")
    print(f"  RMSE:     {rmse:.3f} kcal/mol")

# -------------------------------------------------------
# SECTION 5: SCREEN 500 NEW DRUGS
# -------------------------------------------------------

print("\n" + "="*70)
print("STEP 5: SCREENING 500 NEW DRUGS USING ML")
print("="*70)
print("Downloading 500 new FDA approved drugs from ChEMBL...")

# Download new drugs not in our original set
new_drug_url = (
    f"{base_url}/molecule?"
    f"max_phase=4&"
    f"molecule_type=Small+molecule&"
    f"format=json&"
    f"limit=600&"
    f"offset=150"      # Skip first 150 (already used)
)

new_response = requests.get(new_drug_url)
new_molecules = new_response.json()["molecules"]

new_drug_list = []
for mol in new_molecules:
    name = mol.get("pref_name", "Unknown")
    cid = mol.get("molecule_chembl_id", "")
    props = mol.get("molecule_properties", {})
    structs = mol.get("molecule_structures", {})

    if not props or not structs:
        continue

    smiles = structs.get("canonical_smiles", "")
    if not smiles:
        continue

    try:
        mw   = float(props.get("full_mwt", 999))
        logp = float(props.get("alogp", 99))
        hbd  = int(props.get("hbd", 99))
        hba  = int(props.get("hba", 99))

        if (mw < 500 and logp < 5 and
                hbd <= 5 and hba <= 10 and mw > 150):
            new_drug_list.append({
                "chembl_id": cid,
                "name": name,
                "smiles": smiles,
                "mw": mw,
                "logp": logp
            })
    except:
        continue

print(f"New drugs downloaded: {len(new_drug_list)}")

# Generate fingerprints for new drugs
print("Generating fingerprints for new drugs...")
new_fps = []
new_valid = []

for drug in new_drug_list:
    fp = get_fingerprints(drug["smiles"])
    if fp is not None:
        new_fps.append(fp)
        new_valid.append(drug)

new_X = np.array(new_fps)
new_df = pd.DataFrame(new_valid)
print(f"New drugs with fingerprints: {len(new_df)}")

# Predict binding for all new drugs
print("\nPredicting binding energies using ML models...")

predictions = {}
for target, model in models.items():
    preds = model.predict(new_X)
    predictions[target] = preds
    print(f"  {target}: predicted {len(preds)} binding energies")

# Add predictions to dataframe
for target, preds in predictions.items():
    new_df[f'pred_{target}'] = preds

# Calculate average predicted binding
new_df['avg_predicted'] = new_df[
    [f'pred_{t}' for t in pivot.columns]
].mean(axis=1)

# Rank by average predicted binding
new_df = new_df.sort_values('avg_predicted')

print(f"\nTop 10 ML-predicted multi-target candidates:")
print(new_df[['name', 'avg_predicted'] +
             [f'pred_{t}' for t in pivot.columns]
             ].head(10).to_string())

# Save predictions
new_df.to_csv("./results/ml_predicted_candidates.csv", index=False)

# -------------------------------------------------------
# SECTION 6: VISUALIZATIONS
# -------------------------------------------------------

print("\nCreating ML visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Actual vs Predicted for each target
for i, (target, perf) in enumerate(performance.items()):
    ax = axes[i//2, i%2]
    ax.scatter(perf['y_test'], perf['y_pred'],
               alpha=0.6, edgecolors='none', s=30)

    # Perfect prediction line
    mn = min(perf['y_test'].min(), perf['y_pred'].min())
    mx = max(perf['y_test'].max(), perf['y_pred'].max())
    ax.plot([mn, mx], [mn, mx], 'r--', lw=2, label='Perfect')

    ax.set_xlabel('Actual Binding Energy (kcal/mol)')
    ax.set_ylabel('Predicted Binding Energy (kcal/mol)')
    ax.set_title(f'{target}\nR² = {perf["r2"]:.3f}')
    ax.legend(fontsize=9)

plt.suptitle(
    'ML Model Performance: Actual vs Predicted Binding Energy',
    fontsize=14, fontweight='bold'
)
plt.tight_layout()
plt.savefig("./plots/ml_performance.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved: plots/ml_performance.png")

# Plot 2: Top ML predicted candidates
fig, ax = plt.subplots(figsize=(12, 6))

top20 = new_df.head(20)
x = range(len(top20))
ax.bar(x, top20['avg_predicted'], color='steelblue', alpha=0.8)
ax.axhline(y=-7, color='red', linestyle='--', label='Strong binding (-7)')
ax.axhline(y=-6, color='orange', linestyle='--', label='Moderate binding (-6)')
ax.set_xticks(x)
ax.set_xticklabels(top20['name'], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Average Predicted Binding Energy (kcal/mol)')
ax.set_title(
    'Top 20 ML-Predicted Multi-Target Drug Candidates\n'
    'for Alzheimer\'s Disease',
    fontsize=13, fontweight='bold'
)
ax.legend()
plt.tight_layout()
plt.savefig("./plots/ml_top_candidates.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved: plots/ml_top_candidates.png")

# -------------------------------------------------------
# SECTION 7: FINAL SUMMARY
# -------------------------------------------------------

print("\n" + "="*70)
print("COMPLETE SUMMARY")
print("="*70)

avg_r2 = np.mean([p['r2'] for p in performance.values()])

print(f"""
ML MODEL PERFORMANCE:
  Features used: {X.shape[1]} molecular fingerprint bits
  Training drugs: {len(common_drugs)}
  Average R²: {avg_r2:.3f}

  Per target:""")

for target, perf in performance.items():
    print(f"    {target:8s}: R²={perf['r2']:.3f}, RMSE={perf['rmse']:.3f}")

print(f"""
SCREENING RESULTS:
  New drugs screened by ML: {len(new_df)}
  Top candidate: {new_df.iloc[0]['name']}
  Best avg predicted energy: {new_df.iloc[0]['avg_predicted']:.3f} kcal/mol

FILES SAVED:
  results/ml_predicted_candidates.csv
  plots/ml_performance.png
  plots/ml_top_candidates.png
""")

print("✅ STEP 19 COMPLETE!")
print("\nNext: Literature validation of top candidates")