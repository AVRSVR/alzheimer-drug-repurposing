# This is a multi-line comment. It explains what the whole file does.
# We use triple quotes """ at the beginning of scripts to describe them.
"""
PHASE 1 - STEP 1: Download Alzheimer's Disease Brain Data
This is our very first script.
We will download real human brain gene expression data.
"""

# ===================================================================
# SECTION 1: IMPORTING LIBRARIES (Tools)
# ===================================================================

import GEOparse          # This library helps us download data from GEO (NCBI database)
import pandas as pd      # pandas helps us work with tables (like Excel). We use 'pd' as short name
import os                # os library helps us create folders and work with files on computer


# ===================================================================
# SECTION 2: CREATING FOLDERS TO ORGANIZE OUR PROJECT
# ===================================================================

# os.makedirs() creates a new folder. exist_ok=True means "don't show error if folder already exists"
os.makedirs("data", exist_ok=True)        # Creates a folder called "data" to store raw downloaded files
os.makedirs("results", exist_ok=True)     # Creates a folder called "results" to store our analysis outputs


# ===================================================================
# SECTION 3: PRINTING MESSAGES TO USER (For learning and tracking)
# ===================================================================

# print() is used to show messages on the terminal screen
print("="*70)                                   # Prints a line of = signs to make output look clean
print("STARTING DOWNLOAD OF GSE5281 (Alzheimer's Brain Data)")
print("="*70)
print("\nThis may take 3-10 minutes depending on your internet...\n")


# ===================================================================
# SECTION 4: DOWNLOADING THE ACTUAL DATA
# ===================================================================

# GEOparse.get_GEO() connects to NCBI website and downloads the dataset
# geo="GSE5281" is the ID of the Alzheimer's dataset we want
# destdir="./data/" tells it to save the downloaded file in our data folder
gse = GEOparse.get_GEO(geo="GSE5281", destdir="./data/")


# ===================================================================
# SECTION 5: SHOWING SUCCESS MESSAGE
# ===================================================================

print("\n✅ Download completed successfully!")
print(f"Number of samples: {len(gse.gsms)}")                    # len() tells us how many samples (patients) are in the data
print(f"Number of genes/probes: {len(list(gse.gpls.values())[0].table)}")


# ===================================================================
# SECTION 6: BUILDING THE EXPRESSION MATRIX (Main Data Table)
# ===================================================================

print("\nCreating expression matrix...")

# We create an empty dictionary to store data temporarily
expression_data = {}

# This loop goes through every sample in the downloaded data
for gsm_name, gsm in gse.gsms.items():
    # For each sample, we take the gene ID and expression value and store it
    expression_data[gsm_name] = gsm.table.set_index("ID_REF")["VALUE"]

# Convert the dictionary into a proper table (DataFrame) using pandas
df = pd.DataFrame(expression_data)

# Show the shape of our table (rows = genes, columns = samples)
print(f"Expression matrix shape: {df.shape} (genes x samples)")


# ===================================================================
# SECTION 7: SAVING OUR DATA TO FILES
# ===================================================================

# Save the big table as a CSV file so we can use it later
df.to_csv("./data/raw_expression_matrix.csv")
print("\nRaw data saved to: data/raw_expression_matrix.csv")

# Also save basic information about each sample (title, source, etc.)
sample_info = []
for gsm_name, gsm in gse.gsms.items():
    title = gsm.metadata.get('title', ['Unknown'])[0]
    source = gsm.metadata.get('source_name_ch1', ['Unknown'])[0]
    sample_info.append([gsm_name, title, source])

# Convert sample info into a table and save it
sample_df = pd.DataFrame(sample_info, columns=['Sample_ID', 'Title', 'Source'])
sample_df.to_csv("./results/sample_info.csv", index=False)
print("Sample information saved to: results/sample_info.csv")


# ===================================================================
# SECTION 8: FINAL MESSAGE
# ===================================================================

print("\n" + "="*70)
print("STEP 1 COMPLETE!")
print("Next step: We will label which samples are AD and which are Normal.")
print("="*70)