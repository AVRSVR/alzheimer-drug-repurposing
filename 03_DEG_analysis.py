"""
PHASE 1 - STEP 3: Differential Gene Expression Analysis
We compare AD vs Normal expression levels gene by gene.
"""

# Import libraries
import pandas as pd                                        # For working with tables
import numpy as np                                         # For math operations
from scipy import stats                                    # For t-test
from statsmodels.stats.multitest import multipletests      # For FDR correction

# -------------------------------------------------------
# LOAD THE DATA
# -------------------------------------------------------

# Load normalized expression matrix (genes x samples)
expr = pd.read_csv("./results/normalized_expression.csv", index_col=0)

# Load labeled sample file (which sample is AD or Normal)
labels = pd.read_csv("./results/labeled_samples.csv")

print("Expression matrix shape:", expr.shape)
print("Number of labeled samples:", len(labels))

# -------------------------------------------------------
# SEPARATE SAMPLES INTO AD AND NORMAL GROUPS
# -------------------------------------------------------

# Get list of AD sample IDs
ad_samples = labels[labels["Condition"] == "AD"]["Sample_ID"].tolist()

# Get list of Normal sample IDs
normal_samples = labels[labels["Condition"] == "Normal"]["Sample_ID"].tolist()

print(f"\nAD samples found: {len(ad_samples)}")
print(f"Normal samples found: {len(normal_samples)}")

# Keep only samples that exist in expression matrix
ad_samples = [s for s in ad_samples if s in expr.columns]
normal_samples = [s for s in normal_samples if s in expr.columns]

print(f"AD samples confirmed in matrix: {len(ad_samples)}")
print(f"Normal samples confirmed in matrix: {len(normal_samples)}")

# -------------------------------------------------------
# REMOVE CONTROL PROBES BEFORE ANALYSIS
# -------------------------------------------------------

# AFFX probes are technical controls, not real human genes
# We remove them before doing any statistics
# ~ means NOT, startswith() checks the beginning of text
expr = expr[~expr.index.str.startswith("AFFX")]
print(f"\nGenes after removing AFFX controls: {len(expr)}")

# -------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANALYSIS
# -------------------------------------------------------

print("\nStarting differential expression analysis...")
print("This takes 2-3 minutes...")

results = []    # Empty list to store results

# Loop through every gene
for gene in expr.index:

    # Get expression values for this gene in each group
    ad_values = expr.loc[gene, ad_samples]
    normal_values = expr.loc[gene, normal_samples]

    # Calculate mean expression in each group
    mean_ad = ad_values.mean()
    mean_normal = normal_values.mean()

    # Calculate log2 fold change
    # Positive = higher in AD (upregulated)
    # Negative = lower in AD (downregulated)
    log2fc = mean_ad - mean_normal

    # Welch's t-test
    # Tests if the difference is statistically real
    # equal_var=False means we don't assume equal variance
    t_stat, p_value = stats.ttest_ind(
        ad_values,
        normal_values,
        equal_var=False
    )

    # Store this gene's results
    results.append([gene, mean_ad, mean_normal, log2fc, p_value])

print("Finished testing all genes.")

# Convert list into a proper table
results_df = pd.DataFrame(results, columns=[
    "Gene",
    "Mean_AD",
    "Mean_Normal",
    "log2FoldChange",
    "p_value"
])

# -------------------------------------------------------
# MULTIPLE TESTING CORRECTION
# -------------------------------------------------------

# Adjust all p-values together using FDR correction
_, adj_pvalues, _, _ = multipletests(
    results_df["p_value"],
    method="fdr_bh"
)

# Add adjusted p-values as new column
results_df["adj_p_value"] = adj_pvalues

# -------------------------------------------------------
# FILTER SIGNIFICANT GENES
# -------------------------------------------------------

# Keep only genes that pass BOTH filters
sig_genes = results_df[
    (results_df["adj_p_value"] < 0.05) &
    (abs(results_df["log2FoldChange"]) > 0.2)
]

print(f"\nTotal genes tested: {len(results_df)}")
print(f"Significant DEGs found: {len(sig_genes)}")

print(f"\nTop 10 upregulated in AD:")
print(sig_genes.nlargest(10, "log2FoldChange")[
    ["Gene", "log2FoldChange", "adj_p_value"]
])

print(f"\nTop 10 downregulated in AD:")
print(sig_genes.nsmallest(10, "log2FoldChange")[
    ["Gene", "log2FoldChange", "adj_p_value"]
])

# -------------------------------------------------------
# SAVE RESULTS
# -------------------------------------------------------

results_df.to_csv("./results/DEG_all_results.csv", index=False)
sig_genes.to_csv("./results/DEG_significant.csv", index=False)

print("\nSaved: results/DEG_all_results.csv")
print("Saved: results/DEG_significant.csv")
print("\n✅ DEG ANALYSIS COMPLETE!")