"""
PHASE 1 - STEP 4: Volcano Plot
================================
We visualize all 54,613 genes at once.
Each dot is one gene.
Position shows fold change and significance.
Color shows direction of change.
"""

# matplotlib is the main plotting library in Python
import matplotlib.pyplot as plt

# numpy for math operations
import numpy as np

# pandas for loading our results table
import pandas as pd

import os
os.makedirs("./plots", exist_ok=True)    # Create plots folder if it doesn't exist

# -------------------------------------------------------
# LOAD RESULTS FROM PREVIOUS STEP
# -------------------------------------------------------

# Load the full results table (all 54,613 genes)
results = pd.read_csv("./results/DEG_all_results.csv")

print("Loaded results for", len(results), "genes")

# -------------------------------------------------------
# CALCULATE Y AXIS VALUES
# -------------------------------------------------------

# Y axis = -log10(adj_p_value)
# Why negative log10?
# Small p-value (0.0001) = very significant
# -log10(0.0001) = 4  → appears HIGH on plot
# Large p-value (0.9) = not significant
# -log10(0.9) = 0.05 → appears LOW on plot
# So significant genes float to the TOP

# We add a tiny number (1e-300) to avoid log10(0) = error
results["neg_log10_p"] = -np.log10(results["adj_p_value"] + 1e-300)

# -------------------------------------------------------
# ASSIGN COLORS TO EACH GENE
# -------------------------------------------------------

# Create a color list for all genes
# Default color is grey (not significant)
colors = []

for index, row in results.iterrows():

    # Check if gene is significantly UPREGULATED
    # Must pass BOTH filters: p-value AND fold change
    if row["adj_p_value"] < 0.05 and row["log2FoldChange"] > 0.2:
        colors.append("red")       # Upregulated = red

    # Check if gene is significantly DOWNREGULATED
    elif row["adj_p_value"] < 0.05 and row["log2FoldChange"] < -0.2:
        colors.append("blue")      # Downregulated = blue

    # Not significant
    else:
        colors.append("grey")      # Not significant = grey

# Count how many of each color
red_count = colors.count("red")
blue_count = colors.count("blue")
grey_count = colors.count("grey")

print(f"Upregulated (red):   {red_count}")
print(f"Downregulated (blue): {blue_count}")
print(f"Not significant (grey): {grey_count}")

# -------------------------------------------------------
# CREATE THE VOLCANO PLOT
# -------------------------------------------------------

# Create a figure with specific size (width=12, height=8 inches)
plt.figure(figsize=(12, 8))

# Draw all dots
# x = fold change, y = significance, c = colors we defined
plt.scatter(
    results["log2FoldChange"],    # X axis values
    results["neg_log10_p"],       # Y axis values
    c=colors,                     # Colors we assigned
    alpha=0.5,                    # Transparency (0=invisible, 1=solid)
    s=10                          # Size of each dot
)

# -------------------------------------------------------
# ADD THRESHOLD LINES
# -------------------------------------------------------

# Horizontal line at p=0.05
# -log10(0.05) = 1.3
# Genes above this line are statistically significant
plt.axhline(
    y=-np.log10(0.05),      # Y position of line
    color="black",           # Line color
    linestyle="--",          # Dashed line
    linewidth=1,             # Line thickness
    label="p = 0.05"         # Label for legend
)

# Vertical line at log2FC = +0.2 (upregulation threshold)
plt.axvline(
    x=0.2,
    color="red",
    linestyle="--",
    linewidth=1,
    alpha=0.5
)

# Vertical line at log2FC = -0.2 (downregulation threshold)
plt.axvline(
    x=-0.2,
    color="blue",
    linestyle="--",
    linewidth=1,
    alpha=0.5
)

# -------------------------------------------------------
# ADD LABELS AND TITLE
# -------------------------------------------------------

plt.xlabel("log2 Fold Change (AD vs Normal)", fontsize=13)
plt.ylabel("-log10 (adjusted p-value)", fontsize=13)
plt.title(
    "Volcano Plot: Alzheimer's Disease vs Normal Brain\n"
    f"Red = Upregulated ({red_count} genes)  |  "
    f"Blue = Downregulated ({blue_count} genes)",
    fontsize=14,
    fontweight="bold"
)

# -------------------------------------------------------
# SAVE AND SHOW THE PLOT
# -------------------------------------------------------

plt.tight_layout()     # Automatically adjusts spacing

# Save as high quality image
plt.savefig("./plots/volcano_plot.png", dpi=300, bbox_inches="tight")
print("\nPlot saved to: plots/volcano_plot.png")

plt.show()             # Display the plot on screen

print("\n✅ VOLCANO PLOT COMPLETE!")