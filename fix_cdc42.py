# fix_cdc42.py
print("Cleaning CDC42 PDB file...")

with open("./structures/CDC42_pdb.pdb", "r") as f:
    lines = f.readlines()

# Keep only ATOM records (protein atoms)
# Remove HETATM (GDP, MG) and water
clean_lines = []
for line in lines:
    if line.startswith("ATOM"):
        clean_lines.append(line)
    elif line.startswith("END"):
        clean_lines.append(line)

# Save cleaned PDB
with open("./structures/CDC42_clean.pdb", "w") as f:
    f.writelines(clean_lines)

print(f"Original lines: {len(lines)}")
print(f"Clean lines: {len(clean_lines)}")
print("Saved: CDC42_clean.pdb")