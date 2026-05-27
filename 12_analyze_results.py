"""
PHASE 2 - STEP 12: Analyze and Visualize Docking Results
=========================================================
We merge drug names with docking scores,
create visualizations, and identify top candidates.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

os.makedirs("./plots", exist_ok=True)

print("="*70)
print("ANALYZING DOCKING RESULTS")
print("="*70)

# -------------------------------------------------------
# SECTION 1: LOAD AND MERGE DATA
# -------------------------------------------------------

# Load docking results (has ChEMBL IDs and affinities)
docking = pd.read_csv("./results/docking_results.csv")

# Load drug candidates (has ChEMBL IDs and names)
drugs = pd.read_csv("./drugs/drug_candidates.csv")

# Rename column so we can merge
# docking has "drug" column with ChEMBL ID
# drugs has "chembl_id" column
docking = docking.rename(columns={"drug": "chembl_id"})

# Merge the two tables on ChEMBL ID
# Like VLOOKUP in Excel
merged = docking.merge(
    drugs[["chembl_id", "name", 
           "molecular_weight", "logp"]],
    on="chembl_id",
    how="left"
)

print("\nFull ranked results with drug names:")
print("-"*60)
print(merged[["name", "chembl_id", 
              "best_affinity_kcal_mol"]].to_string())

# -------------------------------------------------------
# SECTION 2: IDENTIFY TOP CANDIDATES
# -------------------------------------------------------

print("\n" + "="*70)
print("TOP 5 DRUG CANDIDATES FOR ITPKB:")
print("="*70)

top5 = merged.head(5)

for i, row in top5.iterrows():
    print(f"\nRank {row['chembl_id']}:")
    print(f"  Drug:     {row['name']}")
    print(f"  Affinity: {row['best_affinity_kcal_mol']} kcal/mol")
    print(f"  MW:       {row['molecular_weight']} g/mol")
    print(f"  LogP:     {row['logp']}")

# -------------------------------------------------------
# SECTION 3: BAR CHART OF ALL RESULTS
# -------------------------------------------------------

print("\nCreating visualization...")

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# ── Plot 1: Bar chart of binding affinities ──
# Sort for plotting (most negative at top)
plot_data = merged.sort_values(
    "best_affinity_kcal_mol",
    ascending=True
)

# Color top 3 differently
colors = []
for i, val in enumerate(plot_data["best_affinity_kcal_mol"]):
    if i < 2:
        colors.append("#E74C3C")    # Red for top 2
    elif i < 5:
        colors.append("#E67E22")    # Orange for top 3-5
    else:
        colors.append("#3498DB")    # Blue for rest

# Create horizontal bar chart
axes[0].barh(
    plot_data["name"],                        # Drug names on Y axis
    abs(plot_data["best_affinity_kcal_mol"]), # Absolute affinity on X
    color=colors,
    edgecolor="white",
    linewidth=0.5
)

axes[0].set_xlabel(
    "Binding Affinity |kcal/mol|\n(larger = stronger binding)",
    fontsize=11
)
axes[0].set_title(
    "Drug Binding Affinity to ITPKB\n(Alzheimer's Target)",
    fontsize=13,
    fontweight="bold"
)

# Add value labels on bars
for i, (val, name) in enumerate(zip(
    plot_data["best_affinity_kcal_mol"],
    plot_data["name"]
)):
    axes[0].text(
        abs(val) + 0.05,
        i,
        f"{val}",
        va="center",
        fontsize=9
    )

# ── Plot 2: Scatter plot (MW vs Affinity) ──
# This shows if molecular weight affects binding

scatter_colors = [
    "#E74C3C" if v <= -8.0
    else "#E67E22" if v <= -7.0
    else "#3498DB"
    for v in merged["best_affinity_kcal_mol"]
]

axes[1].scatter(
    merged["molecular_weight"],
    merged["best_affinity_kcal_mol"],
    c=scatter_colors,
    s=100,
    alpha=0.8,
    edgecolors="white",
    linewidth=0.5
)

# Label each point with drug name
for _, row in merged.iterrows():
    axes[1].annotate(
        row["name"],
        xy=(row["molecular_weight"],
            row["best_affinity_kcal_mol"]),
        fontsize=7,
        ha="left",
        va="bottom"
    )

axes[1].set_xlabel("Molecular Weight (g/mol)", fontsize=11)
axes[1].set_ylabel("Binding Affinity (kcal/mol)", fontsize=11)
axes[1].set_title(
    "Molecular Weight vs Binding Affinity\n"
    "Red = Strong binders | Blue = Weak binders",
    fontsize=13,
    fontweight="bold"
)

# Add horizontal reference line at -7.0
axes[1].axhline(
    y=-7.0,
    color="black",
    linestyle="--",
    linewidth=1,
    alpha=0.5,
    label="Threshold (-7.0 kcal/mol)"
)
axes[1].legend(fontsize=9)

plt.tight_layout()
plt.savefig(
    "./plots/docking_results.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Plot saved: plots/docking_results.png")

# -------------------------------------------------------
# SECTION 4: SAVE FINAL RESULTS
# -------------------------------------------------------

merged.to_csv(
    "./results/final_docking_results.csv",
    index=False
)

print("\nFinal results saved: results/final_docking_results.csv")

# -------------------------------------------------------
# SECTION 5: SUMMARY FOR REPORT
# -------------------------------------------------------

print("\n" + "="*70)
print("SUMMARY FOR YOUR PROJECT REPORT:")
print("="*70)

best_drug = merged.iloc[0]
print(f"""
Target Protein:  ITPKB (Inositol-trisphosphate 3-kinase B)
                 Upregulated in Alzheimer's Disease brain

Drugs Screened:  {len(merged)} FDA approved small molecules

Top Candidate:   {best_drug['name']}
                 Binding affinity: {best_drug['best_affinity_kcal_mol']} kcal/mol
                 Molecular weight: {best_drug['molecular_weight']} g/mol
                 LogP: {best_drug['logp']}

Interpretation:
  Binding affinity of {best_drug['best_affinity_kcal_mol']} kcal/mol indicates
  strong predicted binding to ITPKB.
  This suggests {best_drug['name']} may be a
  repurposing candidate for Alzheimer's Disease
  via ITPKB inhibition.

Next Steps:
  1. Visual confirmation in PyRx
  2. Network analysis for multi-target screening
  3. Literature validation
""")

print("✅ STEP 12 COMPLETE!")
print("Next: Open top 3 drugs in PyRx for visual confirmation")