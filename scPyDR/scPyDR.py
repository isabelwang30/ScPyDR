#!/usr/bin/env python

"""
Command-line script to reduce the dimensionality and visualize 10x Genomics scRNA-seq data

Similar to scanpy's `scanpy.pp.pca`, `scanpy.tl.tsne`, and `scanpy.tl.umap`.
"""

import os
import sys
import matplotlib.pyplot as plt
import argparse
import numpy as np
import pandas as pd
from . import utils as utils
from scPyDR import __version__

def main():
    sys.stdout.write("Welcome to scPyDR! Starting analysis.\n\n")

    # -------------------- parse command-line arguments into program --------------------

    parser = argparse.ArgumentParser(
        prog="scPyDR",
        description="Command-line script to reduce the dimensionality and visualize 10x Genomics scRNA-seq data"
    )

    # input (required)
    parser.add_argument("datadir", \
                        help="Directory containing 10x Genomics scRNA-seq data files.", \
                        metavar="DIR", \
                        type=str)

    # output (optional)
    parser.add_argument("-o", "--output", \
                        help="Output directory to store results. Default: working directory.", \
                        type=str, \
                        required=False)
    
    # other options (optional)
    parser.add_argument("-g", "--min_genes", \
                        help="Minimum number of genes expressed per cell. Default: 200.", \
                        type=int, \
                        default=200,  # Assign default value
                        required=False)
    parser.add_argument("-c", "--min_cells", \
                        help="Minimum number of cells expressing a gene. Default: 5.", \
                        type=int, \
                        default=5,  # Assign default value
                        required=False)
    parser.add_argument("-cr", "--min_cell_reads", \
                        help="Minimum number of reads per cell. Default: None.", \
                        type=int, \
                        required=False)
    parser.add_argument("-gc", "--min_gene_counts", \
                        help="Minimum number of counts per gene. Default: None.", \
                        type=int, \
                        required=False)
    parser.add_argument("-ntop", "--n_top_genes", \
                        help="Number of highly variable genes to keep. Default: 500.", \
                        type=int, \
                        default=500,  # Assign default value
                        required=False)
    parser.add_argument("-t", "--target_sum", \
                        help="Number of reads per cell for normalization. Default: 1e4.", \
                        type=float, \
                        default=1e4,  # Assign default value
                        required=False)
    parser.add_argument("-n", "--nComp", \
                        help="Number of principal components. Default: for n data points and m features, there are min(n-1,m) PCs.", \
                        type=int, \
                        required=False)
    parser.add_argument("--version", \
                        help="Print the version of scPyDR.", \
                        action="version", \
                        version='{version}'.format(version=__version__))
    parser.add_argument("-u", "--umap", 
                        help="Run UMAP for dimensionality reduction and visualization.", 
                        action="store_true", 
                        required=False)

    # parse args
    args = parser.parse_args()

    # validate args
    if not os.path.isdir(args.datadir):
        utils.ERROR("Invalid directory path to 10x Genomics data files.")

    if args.output and not os.path.isdir(args.output):
        utils.ERROR("Invalid output directory to store results.")

    if args.min_genes is not None and args.min_genes < 0:
        utils.ERROR("Minimum number of genes expressed per cell (-g, --min_genes) must be greater than 0.")

    if args.min_cells is not None and args.min_cells < 0:
        utils.ERROR("Minimum number of cells expressing a gene (-c, --min_cells) must be greater than 0.")

    if args.min_cell_reads is not None and args.min_cell_reads < 0:
        utils.ERROR("Minimum number of reads per cell (-cr, --min_cell_reads) must be greater than 0.")

    if args.min_gene_counts is not None and args.min_gene_counts < 0:
        utils.ERROR("Minimum number of counts per gene (-gc, --min_gene_counts) must be greater than 0.")

    if args.n_top_genes is not None and args.n_top_genes < 0:
        utils.ERROR("Number of highly variable genes to keep (-ntop, --n_top_genes) must be greater than 0.")

    if args.target_sum is not None and args.target_sum < 0:
        utils.ERROR("Number of reads per cell for normalization (-t, --target_sum) must be greater than 0.")

    if args.nComp is not None and args.nComp < 0:
        utils.ERROR("Number of principal components (-n, --nComp) must be greater than 0.")

    # -------------------- set up output and start analysis --------------------

    if args.output is None:
        outdir = os.getcwd()
    else:
        outdir = args.output

    if outdir[-1] != "/":
        outdir += "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # -------------------- load and preprocess data --------------------

    sys.stdout.write("Loading 10x Genomics data files...\n")

    try:
        adata = utils.load(args.datadir, prefix="", cache=True)
    except Exception:
        utils.ERROR(">> Failed to load 10x Genomics data into an AnnData object. Please check: \n 1) path to directory \n 2) directory contains properly formatted 10x Genomics files. \n Read more on 10x Genomics files here: https://www.10xgenomics.com/support/software/space-ranger/latest/advanced/hdf5-feature-barcode-matrix-format")
    try:
        adatac = utils.preprocess(adata, args.min_genes, args.min_cells, args.min_cell_reads,
                                 args.min_gene_counts, args.n_top_genes, args.target_sum)
    except Exception:
        utils.ERROR(">> Failed to preprocess data. Check arguments and try again.")
    try:
        df = utils.convert(adatac)
    except Exception:
        utils.ERROR(">> Failed to convert AnnData object to DataFrame. This should not happen; please make a pull request to scPyDR.")
    sys.stdout.write("Data loaded successfully.\n\n")    

    # -------------------- initialize, fit and transform data using PCA --------------------

    sys.stdout.write("Running PCA for dimensionality reduction... \n")

    nComp = None
    if args.nComp is None:
        # for n data points and m features, there are min(n-1,m) PCs.
        num_pts = df.shape[0] * df.shape[1]
        num_features = df.shape[1]
        nComp = min(num_pts - 1, num_features)
    else:
        nComp = args.nComp

    pca = utils.scpydrPCA(nComp=nComp)  # create pca object
    sys.stdout.write("PCA object created. \n")
    pca.fit(df)  # compute new PCs
    sys.stdout.write("New PCs computed. \n")
    pca_results = pca.transform(df)  # fit data to new PCs
    sys.stdout.write("Original data successfully projected onto the new PCs! \n")

    # -------------------- save and plot PCA results --------------------

    filename_prefix = os.path.basename(args.datadir)

    # Save PCA results to file
    utils.save_pca_results(outdir, filename_prefix, pca_results)

    # Plot PCA results and save the plot to a file
    sys.stdout.write("Plotting top two PCs... \n")
    utils.plot_pca_results(outdir, filename_prefix, pca_results)

    # -------------------- compute and plot UMAP embedding --------------------

    if args.umap:
        sys.stdout.write("Running UMAP for dimensionality reduction and visualization... \n")
        umap_embedding, cluster_labels = utils.umap_embedding(adata)
        # Further actions with umap_embedding if needed
        sys.stdout.write("UMAP computation completed! \n\n")

        # Plot the UMAP embedding
        sys.stdout.write("Plotting UMAP embedding... \n")
        utils.plot_umap_results(outdir, filename_prefix, umap_embedding, cluster_labels)
        sys.stdout.write("Congrats! scPyDR PCA and UMAP were successfully run!\n")
    else:
        sys.stdout.write("Congrats! scPyDR PCA was successfully run!\n")

    sys.stdout.close()
    sys.exit(0)

if __name__ == "__main__":
    main()
