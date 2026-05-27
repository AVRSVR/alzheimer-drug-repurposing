"""
PHASE 3 - STEP 13: Protein Interaction Network Analysis
========================================================
We build a network of AD-related proteins using STRING DB.
We find hub proteins (most connected = best drug targets).
We identify which proteins to use for expanded docking.
"""

import requests          # For talking to STRING API
import pandas as pd      # For handling data tables
import numpy as np       # For math
import networkx as nx    # For building and analyzing networks
import matplotlib.pyplot as plt   # For plotting
import os

os.makedirs("./plots", exist_ok=True)
os.makedirs("./results", exist_ok=True)

print("="*70)
print("PROTEIN-PROTEIN INTERACTION NETWORK ANALYSIS")
print("="*70)

# -------------------------------------------------------
# SECTION 1: LOAD OUR SIGNIFICANT AD GENES
# -------------------------------------------------------

# Load the named DEG results from Phase 1
sig_genes = pd.read_csv("./results/DEG_significant_named.csv")

print(f"\nLoaded {len(sig_genes)} significant AD genes")

# Get gene symbols (removing NaN values)
# We need real gene names for STRING database
gene_symbols = sig_genes["gene_symbol"].dropna().unique().tolist()

# Remove empty strings
gene_symbols = [g for g in gene_symbols if g != ""]

# Some probes map to multiple genes (separated by ///)
# We need to split these and take first gene only
clean_genes = []
for gene in gene_symbols:
    if "///" in gene:
        # Take only first gene name
        first_gene = gene.split("///")[0].strip()
        clean_genes.append(first_gene)
    else:
        clean_genes.append(gene.strip())

# Remove duplicates
clean_genes = list(set(clean_genes))

print(f"Unique gene symbols: {len(clean_genes)}")
print(f"First 10 genes: {clean_genes[:10]}")

# -------------------------------------------------------
# SECTION 2: GET INTERACTIONS FROM STRING DATABASE
# -------------------------------------------------------

print("\n" + "-"*50)
print("QUERYING STRING DATABASE FOR INTERACTIONS")
print("-"*50)
print("This connects to string-db.org...")

# STRING API endpoint
string_url = "https://string-db.org/api/json/network"

# We can only send a limited number of genes at once
# Take top 50 genes for network analysis
genes_to_query = clean_genes    # Use ALL genes

# Join gene names with newline for API
genes_str = "%0d".join(genes_to_query)

# Query parameters
params = {
    "identifiers": genes_str,    # Our gene list
    "species": 9606,             # 9606 = Homo sapiens
    "required_score": 300,       # Minimum confidence score
                                 # 400 = medium confidence
                                 # 700 = high confidence
    "caller_identity": "alzheimer_pipeline"
}

response = requests.post(string_url, data=params)

if response.status_code == 200:
    interactions = response.json()
    print(f"Found {len(interactions)} protein interactions")
else:
    print(f"STRING query failed: {response.status_code}")
    interactions = []

# -------------------------------------------------------
# SECTION 3: BUILD THE NETWORK
# -------------------------------------------------------

print("\n" + "-"*50)
print("BUILDING PROTEIN INTERACTION NETWORK")
print("-"*50)

# networkx builds and analyzes networks (graphs)
# A network has:
# → Nodes = proteins
# → Edges = interactions between proteins

G = nx.Graph()    # G = our network graph

# Add each interaction as an edge between two proteins
for interaction in interactions:
    protein_a = interaction["preferredName_A"]
    protein_b = interaction["preferredName_B"]
    score = interaction["score"]        # Confidence score

    # Add edge (connection) between the two proteins
    # weight = confidence of interaction
    G.add_edge(protein_a, protein_b, weight=score)

print(f"Network built:")
print(f"  Nodes (proteins): {G.number_of_nodes()}")
print(f"  Edges (interactions): {G.number_of_edges()}")

# -------------------------------------------------------
# SECTION 4: FIND HUB PROTEINS
# -------------------------------------------------------

print("\n" + "-"*50)
print("FINDING HUB PROTEINS")
print("-"*50)

# Degree centrality = how many connections each protein has
# Higher degree = more connected = more important hub

degree_centrality = nx.degree_centrality(G)

# Betweenness centrality = how often a protein sits
# on the shortest path between other proteins
# Higher = more important bridge in network
betweenness_centrality = nx.betweenness_centrality(G)

# Combine into a DataFrame for easy viewing
centrality_df = pd.DataFrame({
    "protein": list(degree_centrality.keys()),
    "degree_centrality": list(degree_centrality.values()),
    "betweenness_centrality": [
        betweenness_centrality.get(p, 0)
        for p in degree_centrality.keys()
    ]
})

# Calculate combined hub score
# Average of both centrality measures
centrality_df["hub_score"] = (
    centrality_df["degree_centrality"] +
    centrality_df["betweenness_centrality"]
) / 2

# Sort by hub score
centrality_df = centrality_df.sort_values(
    "hub_score",
    ascending=False
).reset_index(drop=True)

print("\nTop 10 Hub Proteins in AD Network:")
print("-"*50)
print(centrality_df.head(10).to_string())

# Save hub proteins
centrality_df.to_csv(
    "./results/hub_proteins.csv",
    index=False
)

# -------------------------------------------------------
# SECTION 5: VISUALIZE THE NETWORK
# -------------------------------------------------------

print("\nCreating network visualization...")

plt.figure(figsize=(16, 12))

# Layout = how nodes are positioned
# spring_layout = nodes repel each other
# connected nodes attract each other
pos = nx.spring_layout(G, seed=42, k=2)

# Node size based on degree centrality
# More connected = bigger node
node_sizes = [
    degree_centrality[node] * 3000 + 100
    for node in G.nodes()
]

# Node color based on hub score
node_colors = [
    centrality_df[
        centrality_df["protein"] == node
    ]["hub_score"].values[0]
    if node in centrality_df["protein"].values
    else 0
    for node in G.nodes()
]

# Draw network
nx.draw_networkx_nodes(
    G, pos,
    node_size=node_sizes,
    node_color=node_colors,
    cmap=plt.cm.RdYlBu_r,    # Red = hub, Blue = peripheral
    alpha=0.9
)

nx.draw_networkx_edges(
    G, pos,
    alpha=0.3,
    edge_color="grey",
    width=0.8
)

nx.draw_networkx_labels(
    G, pos,
    font_size=8,
    font_weight="bold"
)

plt.title(
    "Alzheimer's Disease Protein Interaction Network\n"
    "Node size & color = Hub importance "
    "(Red = Hub, Blue = Peripheral)",
    fontsize=14,
    fontweight="bold"
)
plt.axis("off")
plt.tight_layout()
plt.savefig(
    "./plots/PPI_network.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Network plot saved: plots/PPI_network.png")

# -------------------------------------------------------
# SECTION 6: SELECT TOP TARGETS FOR EXPANDED DOCKING
# -------------------------------------------------------

print("\n" + "="*70)
print("TOP HUB PROTEINS SELECTED FOR EXPANDED DOCKING:")
print("="*70)

top_hubs = centrality_df.head(3)

for i, row in top_hubs.iterrows():
    print(f"\n#{i+1}: {row['protein']}")
    print(f"     Hub Score: {row['hub_score']:.4f}")
    print(f"     Degree Centrality: {row['degree_centrality']:.4f}")
    print(f"     Betweenness: {row['betweenness_centrality']:.4f}")

# Save top hubs for next step
top_hub_names = top_hubs["protein"].tolist()
with open("./results/top_hub_proteins.txt", "w") as f:
    f.write("\n".join(top_hub_names))

print(f"\nTop hub proteins saved: results/top_hub_proteins.txt")
print("\n✅ STEP 13 COMPLETE!")
print("Next: Download structures of hub proteins and dock 100 drugs")