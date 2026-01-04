#!/usr/bin/env python3
"""
Author: Davie Slocum
Purpose: Preprocess TCGA Pan-Cancer expression data for clustering + supervised learning:
         - load expression + phenotype
         - keep tumor samples only
         - create cancer-type labels
         - select top variable genes
         - standardize
         - save processed matrix + labels

NOTE: We intentionally DO NOT log-transform this dataset because EB++ adjusted
expression can contain negative values, and log1p would produce NaNs/infs.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

DATA_DIR = Path("data")

EXPR_IN = DATA_DIR / "expression_pancan.tsv.gz"
PHENO_IN = DATA_DIR / "phenotype_tcga.tsv.gz"

EXPR_OUT = DATA_DIR / "expression_processed.csv.gz"
LABELS_OUT = DATA_DIR / "labels.csv"


def main():
    if not EXPR_IN.exists():
        raise FileNotFoundError(f"Missing {EXPR_IN}. Run: python3 scripts/download_xena_data.py")
    if not PHENO_IN.exists():
        raise FileNotFoundError(f"Missing {PHENO_IN}. Run: python3 scripts/download_xena_data.py")

    # 1) Load phenotype table
    pheno = pd.read_csv(PHENO_IN, sep="\t", compression="gzip", low_memory=False)
    print("Phenotype columns:", pheno.columns.tolist())

    # ['sample', 'sample_type_id', 'sample_type', '_primary_disease']
    sample_col = "sample"
    sample_type_col = "sample_type"
    label_col = "_primary_disease"

    # Keep tumor samples only
    pheno_tumor = pheno[pheno[sample_type_col].str.contains("tumor", case=False, na=False)].copy()

    # Labels
    labels = (
        pheno_tumor[[sample_col, label_col]]
        .dropna()
        .drop_duplicates()
        .set_index(sample_col)[label_col]
    )
    print(f"Tumor samples with labels: {labels.shape[0]}")

    # 2) Load expression matrix (genes x samples)
    expr = pd.read_csv(EXPR_IN, sep="\t", compression="gzip", index_col=0)
    print("Expression shape (genes x samples):", expr.shape)

    # 3) Subset expression to labeled tumor samples
    keep_samples = expr.columns.intersection(labels.index)
    expr = expr[keep_samples]
    labels = labels.loc[keep_samples]
    print("After subsetting (genes x samples):", expr.shape)
    print("Labels length:", labels.shape[0])

    # 4) Transpose to samples x genes
    X = expr.T

    # 5) Keep top variable genes
    top_n = 2000
    variances = X.var(axis=0)
    top_genes = variances.sort_values(ascending=False).head(top_n).index
    X = X[top_genes]

    # 6) Convert to numeric safely and clean non-finite values (just in case)
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)

    if X.isna().any().any():
        nan_count = int(X.isna().sum().sum())
        print(f"Found {nan_count} NaN values; filling with 0.")
        X = X.fillna(0.0)

    if not np.isfinite(X.to_numpy()).all():
        raise ValueError("Still found non-finite values in X after cleanup.")

    # 7) Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    X_out = pd.DataFrame(X_scaled, index=X.index, columns=X.columns)

    # 8) Save outputs
    X_out.to_csv(EXPR_OUT, compression="gzip")
    labels.to_csv(LABELS_OUT)

    print(f"Saved:\n  {EXPR_OUT}\n  {LABELS_OUT}")
    print("Preprocessing complete.")


if __name__ == "__main__":
    main()
