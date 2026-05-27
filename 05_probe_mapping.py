"""
PHASE 1 - STEP 5: Convert Probe IDs to Gene Names
===================================================
Our DEG results use probe IDs like "236439_at"
We convert these to real gene names like "APP"
using the platform annotation file (GPL96)
"""

import pandas as pd
import GEOparse
import os

# -------------------------------------------------------
# LOAD THE PLATFORM ANNOTATION FILE
# -------------------------------------------------------

# The platform file tells us:
# probe ID → gene symbol → gene description
# It was downloaded when we got GSE5281

print("Loading platform annotation...")
print("This tells us what each probe measures...")

# Load the already downloaded dataset
# This does NOT download again, reads from disk
gse = GEOparse.get_GEO(
    geo="GSE5281",
    destdir="./data/",
    silent=True          # Don't print download messages
)

# Get the platform (GPL) object
# GPL96 is the Affymetrix Human Genome U133A array
gpl_name = list(gse.gpls.keys())[0]
gpl = gse.gpls[gpl_name]

print(f"Platform found: {gpl_name}")
print(f"Annotation table shape: {gpl.table.shape}")
print(f"\nFirst 3 rows of annotation:")
print(gpl.table.head(3))

# -------------------------------------------------------
# BUILD PROBE TO GENE MAPPING DICTIONARY
# -------------------------------------------------------

# Look at what columns are available
print(f"\nAvailable columns: {gpl.table.columns.tolist()}")
# -------------------------------------------------------
# EXTRACT PROBE TO GENE MAPPING
# -------------------------------------------------------

# We only need two columns:
# ID = probe ID
# Gene Symbol = real gene name
annotation = gpl.table[["ID", "Gene Symbol", "Gene Title"]].copy()

# Rename columns for clarity
annotation.columns = ["probe_id", "gene_symbol", "gene_title"]

# Remove rows where gene symbol is empty
# Some probes don't map to any known gene
annotation = annotation[
    annotation["gene_symbol"].notna() &      # Not NaN
    (annotation["gene_symbol"] != "")        # Not empty string
]

print(f"\nProbes with gene symbols: {len(annotation)}")
print("\nExample mappings:")
print(annotation.head(10).to_string())

# -------------------------------------------------------
# LOAD OUR DEG RESULTS
# -------------------------------------------------------

# Load significant genes from previous step
sig_genes = pd.read_csv("./results/DEG_significant.csv")
all_results = pd.read_csv("./results/DEG_all_results.csv")

print(f"\nSignificant DEGs to map: {len(sig_genes)}")

# -------------------------------------------------------
# MERGE DEG RESULTS WITH GENE NAMES
# -------------------------------------------------------

# merge() combines two tables based on a common column
# Like VLOOKUP in Excel
# left_on = column name in sig_genes table
# right_on = column name in annotation table

sig_mapped = sig_genes.merge(
    annotation,
    left_on="Gene",        # probe ID column in our results
    right_on="probe_id",   # probe ID column in annotation
    how="left"             # keep all sig_genes rows
)

all_mapped = all_results.merge(
    annotation,
    left_on="Gene",
    right_on="probe_id",
    how="left"
)

# -------------------------------------------------------
# SHOW RESULTS WITH REAL GENE NAMES
# -------------------------------------------------------

print("\n" + "="*60)
print("TOP 10 UPREGULATED GENES IN ALZHEIMER'S:")
print("="*60)
top_up = sig_mapped.nlargest(10, "log2FoldChange")[
    ["gene_symbol", "gene_title", "log2FoldChange", "adj_p_value"]
]
print(top_up.to_string())

print("\n" + "="*60)
print("TOP 10 DOWNREGULATED GENES IN ALZHEIMER'S:")
print("="*60)
top_down = sig_mapped.nsmallest(10, "log2FoldChange")[
    ["gene_symbol", "gene_title", "log2FoldChange", "adj_p_value"]
]
print(top_down.to_string())

# -------------------------------------------------------
# SAVE MAPPED RESULTS
# -------------------------------------------------------

sig_mapped.to_csv("./results/DEG_significant_named.csv", index=False)
all_mapped.to_csv("./results/DEG_all_named.csv", index=False)

print("\nSaved: results/DEG_significant_named.csv")
print("Saved: results/DEG_all_named.csv")
print("\n✅ PROBE MAPPING COMPLETE!")