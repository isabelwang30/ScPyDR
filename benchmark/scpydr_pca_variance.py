import scanpy as sc
import argparse
import os
import sys
import numpy as np
import pandas as pd
from scPyDR import utils as utils
from scPyDR import __version__


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run Scanpy PCA on 10x Genomics scRNA-seq data")
    parser.add_argument("datadir", help="Directory containing 10x Genomics scRNA-seq data files.", metavar="DIR", type=str)
    parser.add_argument("-n", "--nComp", help="Number of principal components.", type=int)
    args = parser.parse_args()

    # Load data
    try:
        adata = utils.load(args.datadir, prefix="", cache=True)
    except Exception:
        utils.ERROR(">> Failed to load 10x Genomics data into an AnnData object. Please check: \n 1) path to directory \n 2) directory contains properly formatted 10x Genomics files. \n Read more on 10x Genomics files here: https://www.10xgenomics.com/support/software/space-ranger/latest/advanced/hdf5-feature-barcode-matrix-format")
    try:
        adata = utils.preprocess(adata)
    except Exception:
        utils.ERROR(">> Failed to preprocess data. Check arguments and try again.")
    try:
        df = utils.convert(adata)
    except Exception:
        utils.ERROR(">> Failed to convert AnnData object to DataFrame. This should not happen; please make a pull request to scPyDR.")
        
    
    # Fit data
    pca = utils.scpydrPCA(nComp=args.nComp)  # create pca object
    pca.fit(df)  # compute new PCs
    explained_var_perc_n = round(pca.perc_explained_var, 2)
    print(f"Total explained variance for the first {args.nComp} components: {explained_var_perc_n}%")
    

if __name__ == "__main__":
    main()
