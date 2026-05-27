"""
PHASE 1 - NORMALIZATION STEP
We will:
1. Log2 transform the raw expression data
2. Apply quantile normalization
3. Save the cleaned data
"""

# Import libraries
import pandas as pd
import numpy as np
from sklearn.preprocessing import quantile_transform

# -------------------------------------------------------
# LOAD RAW EXPRESSION DATA
# -------------------------------------------------------

# Load the raw matrix (genes x samples)
expr = pd.read_csv("./data/raw_expression_matrix.csv", index_col=0)

print("Raw data shape:", expr.shape)
print("Max value before log transform:", expr.values.max())
print("Min value before log transform:", expr.values.min())
# -------------------------------------------------------
# STEP 1: LOG2 TRANSFORMATION
# -------------------------------------------------------

# We add 1 before log2 to avoid log2(0) = -infinity
# log2(0 + 1) = log2(1) = 0
# So zero expression stays as zero after transformation
expr_log = np.log2(expr + 1)

print("\nAfter log2 transformation:")
print("Max value:", round(expr_log.values.max(), 2))
print("Min value:", round(expr_log.values.min(), 2))
print("Range is now manageable for statistics")
# -------------------------------------------------------
# STEP 2: QUANTILE NORMALIZATION
# -------------------------------------------------------

# quantile_transform makes all samples have
# identical distribution
# axis=0 means normalize across samples (columns)
# copy=True means don't modify original data

expr_norm = pd.DataFrame(
    quantile_transform(expr_log, axis=0, copy=True),
    index=expr_log.index,        # Keep same gene names as rows
    columns=expr_log.columns     # Keep same sample names as columns
)

print("\nAfter quantile normalization:")
print("Shape:", expr_norm.shape)
print("Max value:", round(expr_norm.values.max(), 2))
print("Min value:", round(expr_norm.values.min(), 2))

# -------------------------------------------------------
# STEP 3: SAVE NORMALIZED DATA
# -------------------------------------------------------

# Save to results folder for use in next steps
expr_norm.to_csv("./results/normalized_expression.csv")
print("\nNormalized data saved to: results/normalized_expression.csv")

print("\n✅ NORMALIZATION COMPLETE!")
print("Now our data is ready for proper statistical analysis.")