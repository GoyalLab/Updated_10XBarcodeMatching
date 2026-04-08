#!/usr/bin/env python3
#==================================================================================================
# stepFour -- identify singlets from barcode extraction output
#
# Wraps the singletCode package (check_sample_sheet + get_singlets) as a CLI script.
# Input: stepThree starcode-shaved reads
# Output: singlet lists, stats, and UMI plots
#
# Usage:
#   python stepFour.py \
#       --input /path/to/stepThreeStarcodeShavedReads.txt \
#       --sample sample_1 \
#       --outdir /path/to/stepFour/
#==================================================================================================

import argparse
import sys

import pandas as pd
from singletCode import check_sample_sheet, get_singlets


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify singlets from stepThree starcode output using singletCode."
    )
    parser.add_argument("--input", required=True,
                        help="Path to stepThreeStarcodeShavedReads.txt (tab-separated)")
    parser.add_argument("--sample", required=True,
                        help="Sample/dataset name (e.g. GFP_barcode_1)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for singletCode results")
    parser.add_argument("--barcode-col", default="BC50StarcodeD8",
                        help="Column name to rename to 'barcode' (default: BC50StarcodeD8)")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load and format input
    print(f"Loading input: {args.input}", flush=True)
    df = pd.read_csv(args.input, sep="\t")

    if args.barcode_col not in df.columns:
        sys.exit(f"ERROR: column '{args.barcode_col}' not found in input. "
                 f"Available columns: {list(df.columns)}")

    df = df.rename(columns={args.barcode_col: "barcode"})
    df["sample"] = args.sample
    df = df[["cellID", "barcode", "sample"]]

    print(f"Sample: {args.sample} | Reads: {len(df):,}", flush=True)

    # Validate
    check_sample_sheet(df)

    # Run singlet identification
    print(f"Running get_singlets...", flush=True)
    cellLabelList, stats = get_singlets(
        df,
        dataset_name=args.sample,
        output_path=args.outdir,
        save_all_singlet_categories=True,
        save_plot_umi=True
    )

    print(f"Stats:\n{stats.to_string()}")
    print(f"Output written to: {args.outdir}")


if __name__ == "__main__":
    main()
