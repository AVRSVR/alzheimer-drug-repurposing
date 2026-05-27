"""
PHASE 1 - STEP 2: Label Samples as AD or Normal
=================================================
We look at sample titles and decide:
- "control"  = Normal healthy brain
- "affected" = Alzheimer's disease brain
"""

# pandas helps us work with tables (like Excel in Python)
import pandas as pd

# =====================================================
# SECTION 1: LOAD THE SAMPLE INFO FILE
# =====================================================

# pd.read_csv() opens a CSV file and turns it into a table
# We saved this file in the previous step
sample_df = pd.read_csv("./results/sample_info.csv")

# Print the first 5 rows so we can see what we are working with
print("="*60)
print("FIRST 5 ROWS OF SAMPLE INFO:")
print("="*60)
print(sample_df.head())           # .head() shows first 5 rows by default
print(f"\nTotal samples: {len(sample_df)}")   # len() counts total rows


# =====================================================
# SECTION 2: CREATE A LABELING FUNCTION
# =====================================================

# A function is a reusable block of code
# We define it once, then use it for every row
def assign_label(title):
    """
    This function takes a sample title as input
    and returns either "AD", "Normal", or "Unknown"
    """

    # .lower() converts text to lowercase
    # so "Control" and "control" both match
    title_lower = title.lower()

    # Check if the word "affected" is in the title
    # "affected" = Alzheimer's sample
    if "affected" in title_lower:
        return "AD"

    # Check if the word "control" is in the title
    # "control" = Normal healthy sample
    elif "control" in title_lower:
        return "Normal"

    # If neither word found, label as Unknown
    else:
        return "Unknown"


# =====================================================
# SECTION 3: APPLY THE FUNCTION TO EVERY ROW
# =====================================================

# .apply() runs our function on every single row
# It takes the Title column and passes each value to assign_label()
sample_df["Condition"] = sample_df["Title"].apply(assign_label)

# =====================================================
# SECTION 4: CHECK OUR RESULTS
# =====================================================

print("\n" + "="*60)
print("LABELING RESULTS:")
print("="*60)

# .value_counts() counts how many of each label we have
print(sample_df["Condition"].value_counts())

# Show any samples we could not label
unknown = sample_df[sample_df["Condition"] == "Unknown"]
print(f"\nUnknown samples: {len(unknown)}")
if len(unknown) > 0:
    print(unknown)

# Show a few examples of each group
print("\nExample AD samples:")
print(sample_df[sample_df["Condition"] == "AD"].head(3))

print("\nExample Normal samples:")
print(sample_df[sample_df["Condition"] == "Normal"].head(3))


# =====================================================
# SECTION 5: SAVE THE LABELED FILE
# =====================================================

# Save with the new Condition column added
sample_df.to_csv("./results/labeled_samples.csv", index=False)

# index=False means don't save the row numbers (0,1,2,3...)
# as an extra column in the file

print("\n" + "="*60)
print("✅ STEP 2 COMPLETE!")
print("Labeled file saved to: results/labeled_samples.csv")
print("="*60)