"""
PHASE 3 - STEP 14 (v2): Download Hub Protein Structures from PDB
==================================================================
We search the Protein Data Bank for experimental structures
of our 3 hub proteins using UniProt IDs.

These proteins have known crystal structures.
We download the best available one.
"""

import requests
import os
import json

os.makedirs("./structures", exist_ok=True)

print("="*70)
print("DOWNLOADING HUB PROTEIN STRUCTURES FROM PDB")
print("="*70)

# -------------------------------------------------------
# DEFINE TARGETS WITH UNIPROT IDs
# -------------------------------------------------------

targets = [
    {
        "name": "CDC42",
        "uniprot_id": "P60953",
        "description": "Cell division cycle 42"
    },
    {
        "name": "HSPA8",
        "uniprot_id": "P11142",
        "description": "Heat shock protein A8"
    },
    {
        "name": "PGK1",
        "uniprot_id": "P00558",
        "description": "Phosphoglycerate kinase 1"
    }
]

# -------------------------------------------------------
# SEARCH PDB FOR EACH TARGET
# -------------------------------------------------------

for target in targets:
    name = target["name"]
    uniprot = target["uniprot_id"]
    
    print(f"\nSearching PDB for {name} ({uniprot})...")
    
    # Search RCSB PDB by UniProt accession
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 5},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]
        }
    }
    
    response = requests.post(search_url, json=query)
    
    if response.status_code == 200:
        data = response.json()
        result_set = data.get("result_set", [])
        
        if result_set:
            pdb_id = result_set[0]["identifier"]
            print(f"  Found PDB ID: {pdb_id}")
            
            # Download PDB file
            download_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            pdb_response = requests.get(download_url)
            
            if pdb_response.status_code == 200:
                output_path = f"./structures/{name}_pdb.pdb"
                with open(output_path, "w") as f:
                    f.write(pdb_response.text)
                
                atom_count = sum(1 for line in pdb_response.text.split("\n") if line.startswith("ATOM"))
                print(f"  ✅ {name} saved: {output_path}")
                print(f"     Atoms: {atom_count}")
            else:
                print(f"  ❌ Download failed (status {pdb_response.status_code})")
        else:
            print(f"  ❌ No PDB structures found for {name}")
    else:
        print(f"  ❌ Search failed (status {response.status_code})")

print("\n✅ STEP 14 COMPLETE!")