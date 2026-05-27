"""
PHASE 2 - STEP 6: Get Protein Sequence from NCBI
================================================
We found ITPKB is upregulated in Alzheimer's.
Now we download its official human protein sequence.
This helps us confirm we are docking against the right target.
"""

from Bio import Entrez, SeqIO     # Biopython tools to talk to NCBI
import os

# Create a folder for our protein data
os.makedirs("protein_data", exist_ok=True)

# -------------------------------------------------------
# SECTION 1: SET UP NCBI ACCESS
# -------------------------------------------------------

# NCBI requires an email address so they know who is downloading
Entrez.email = "your_email@example.com"  # Change this to your email

print("Searching NCBI for ITPKB protein sequence...")

# -------------------------------------------------------
# SECTION 2: SEARCH FOR THE GENE
# -------------------------------------------------------

# Search for ITPKB in the protein database for Humans (Homo sapiens)
# esearch finds the ID of the record
search_handle = Entrez.esearch(
    db="protein", 
    term="ITPKB[Gene Name] AND Homo sapiens[Organism]", 
    retmax=1
)
search_results = Entrez.read(search_handle)
search_handle.close()

# Get the ID of the first result
if search_results["IdList"]:
    protein_id = search_results["IdList"][0]
    print(f"Found Protein ID: {protein_id}")
else:
    print("Could not find protein. Check gene name.")
    exit()

# -------------------------------------------------------
# SECTION 3: DOWNLOAD THE SEQUENCE
# -------------------------------------------------------

# efetch downloads the actual data
# rettype="fasta" means we want it in FASTA format (plain text sequence)
fetch_handle = Entrez.efetch(
    db="protein", 
    id=protein_id, 
    rettype="fasta", 
    retmode="text"
)
protein_record = SeqIO.read(fetch_handle, "fasta")
fetch_handle.close()

# -------------------------------------------------------
# SECTION 4: SAVE AND PRINT
# -------------------------------------------------------

# Save to a file
output_path = "./protein_data/ITPKB_sequence.fasta"
SeqIO.write(protein_record, output_path, "fasta")

print("\n" + "="*60)
print(f"PROTEIN: {protein_record.description}")
print(f"LENGTH:  {len(protein_record.seq)} amino acids")
print("="*60)
print("\nFirst 50 Amino Acids:")
print(protein_record.seq[:50])

print(f"\nSequence saved to: {output_path}")
print("\n✅ STEP 6 COMPLETE!")