
#!/usr/bin/env python3
#==================================================================================================
# stepThree -- starcode to remove barcode sequence errors
#
# Runs starcode on the text files generated in Step 2 with specified distance to collapse barcodes
# Must download and make starcode before running this script:
#     git clone https://github.com/gui11aume/starcode
#     make -C starcode
# Input: stepTwoCellIDUMIBarcodes.txt
# Output: text file that contains collapsed barcodes
#
# Usage:
#   python stepThree.py \
#       --input /path/to/stepTwo/ \
#       --outdir /path/to/stepThree/
#==================================================================================================

import os
import sys
import csv
import numpy as np
import pandas as pd

csv.field_size_limit(100000000)

inputDirectory = sys.argv[1]
outputDirectory = sys.argv[2]

inputFile = os.path.join(inputDirectory, "stepTwoCellIDUMIBarcodes.txt")
if not os.path.exists(inputFile):
    sys.exit("Error: input file not found: " + inputFile)

# Read only the three core columns
df = pd.read_csv(inputFile, sep="\t", usecols=["cellID", "UMI", "BC"], dtype=str)
nrows = len(df)
print("Loaded " + str(nrows) + " reads. Applying starcode replacements...")

# Initialize target arrays with substring defaults from BC
# Using numpy object arrays for fast index-based assignment
bc = df["BC"].values
bc50 = np.array([s[:50] for s in bc], dtype=object)
bc40 = np.array([s[:40] for s in bc], dtype=object)
bc30_d8 = np.array([s[:30] for s in bc], dtype=object)
bc30_d6 = np.array([s[:30] for s in bc], dtype=object)

target_arrays = {
    "BC50StarcodeD8": bc50,
    "BC40StarcodeD8": bc40,
    "BC30StarcodeD8": bc30_d8,
    "BC30StarcodeD6": bc30_d6,
}

starcode_files = {
    "stepThreeBarcodes50_d8": "BC50StarcodeD8",
    "stepThreeBarcodes40_d8": "BC40StarcodeD8",
    "stepThreeBarcodes30_d8": "BC30StarcodeD8",
    "stepThreeBarcodes30_d6": "BC30StarcodeD6",
}

for filename, col in starcode_files.items():
    filepath = os.path.join(outputDirectory, filename)
    if not os.path.exists(filepath):
        print("Warning: starcode output not found, skipping: " + filepath)
        continue

    arr = target_arrays[col]
    with open(filepath, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            sequence = row[0]
            # Parse indices directly as integer array; starcode uses 1-based indexing
            indices = np.array(row[2].split(","), dtype=np.intp)
            indices -= 1
            arr[indices] = sequence

print("Creating starcode replaced table...")

# Assign completed arrays back to dataframe in one operation per column
df["BC50StarcodeD8"] = bc50
df["BC40StarcodeD8"] = bc40
df["BC30StarcodeD8"] = bc30_d8
df["BC30StarcodeD6"] = bc30_d6

col_order = ["cellID", "UMI", "BC", "BC50StarcodeD8", "BC40StarcodeD8",
             "BC30StarcodeD8", "BC30StarcodeD6"]

outputFile = os.path.join(outputDirectory, "stepThreeStarcodeShavedReads.txt")
df[col_order].to_csv(outputFile, sep="\t", index=False)

print("Done. Output written to " + outputFile)
