#!/usr/bin/env python3
"""
Author: Davie Slocum
Purpose: Download TCGA Pan-Cancer gene expression + phenotype from UCSC Xena
         (Pan-Cancer Atlas hub) and save locally.
"""

from pathlib import Path
import urllib.request

DATA_DIR = Path("data")

EXPR_URL = "https://pancanatlas.xenahubs.net/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
PHENO_URL = "https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"

EXPR_OUT = DATA_DIR / "expression_pancan.tsv.gz"
PHENO_OUT = DATA_DIR / "phenotype_tcga.tsv.gz"


def download(url: str, outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    if outpath.exists() and outpath.stat().st_size > 0:
        print(f"Already exists, skipping: {outpath}")
        return
    print(f"Downloading:\n  {url}\n→ {outpath}")
    urllib.request.urlretrieve(url, outpath)
    print(f"Done: {outpath} ({outpath.stat().st_size/1e6:.1f} MB)")


def main():
    download(EXPR_URL, EXPR_OUT)
    download(PHENO_URL, PHENO_OUT)
    print("All downloads complete.")


if __name__ == "__main__":
    main()

