#!/usr/bin/env python

"""
Command-line script to reduce the dimensionality and visualize 10x Genomics scRNA-seq data

Similar to scanpy's `scanpy.pp.pca`, `scanpy.tl.tsne`, and `scanpy.tl.umap`.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from . import utils as utils
from scPyDR import __version__

def main():

    # -------------------- parse command-line arguments into program --------------------

    parser = argparse.ArgumentParser(
        prog = "scPyDR", \
        description = "Command-line script to reduce the dimensionality and visualize 10x Genomics scRNA-seq data"
    )

    # input
    parser.add_argument("datadir", \
                        help = "Directory containing 10x Genomics scRNA-seq data files.", \
                        metavar="DIR", \
                        type = str)

    # output
    parser.add_argument("-o", "--output", \
                        help = "Output directory to store results. Default: working directory.", \
                        type = str, \
                        required = False)
    
    # other options
    parser.add_argument("-g", "--min_genes", \
                        help = "Minimum number of genes expressed per cell. Default: 200.", \
                        type = int, \
                        required = False)
    parser.add_argument("-c", "--min_cells", \
                        help = "Minimum number of cells expressing a gene. Default: 5.", \
                        type = int, \
                        required = False)
    parser.add_argument("-cr", "--min_cell_reads", \
                        help = "Minimum number of reads per cell. Default: None.", \
                        type = int, \
                        required = False)
    parser.add_argument("-gc", "--min_gene_counts", \
                        help = "Minimum number of counts per gene. Default: None.", \
                        type = int, \
                        required = False)
    parser.add_argument("-ntop", "--n_top_genes", \
                        help = "Number of highly variable genes to keep. Default: 500.", \
                        type = int, \
                        required = False)
    parser.add_argument("-t", "--target_sum", \
                        help = "Number of reads per cell for normalization. Default: 1e4.", \
                        type = float, \
                        required = False)
    parser.add_argument("-n", "--nComp", \
                        help = "Number of principal componenets. Default: for n data points and m features, there are min(n-1,m) PCs.", \
                        type = int, \
                        required = False)
    parser.add_argument("--version", \
                        help = "Print the version of scPyDR.", \
                        action = "version", \
                        version = '{version}'.format(version = __version__))
    parser.add_argument("-u", "--umap", 
                        help = "Run UMAP for dimensionality reduction and visualization.", 
                        action="store_true", 
                        required = False)
    # TO DO!
    # parser.add_argument("-v", "--visualize", \
    #                     help = "Visualize the results of scPyDR.", \
    #                     required = False)

    # parse args
    args = parser.parse_args()
    
    # validate args
    if not os.path.isdir(args.datadir):
        print("Invalid directory path to 10x Genomics data files.")
        parser.print_help()
    if not os.path.isdir(args.output):
        print("Invalid output directory to store results.")
        parser.print_help()
    if args.min_genes is not None and args.min_genes < 0:
        print("Minimum number of genes expressed per cell (-g, --min_genes) must be greater than 0.")
    if args.min_cells is not None and args.min_cells < 0:
        print("Minimum number of cells expressing a gene (-c, --min_cells) must be greater than 0.")
    if args.min_cell_reads is not None and args.min_cell_reads < 0:
        print("Minimum number of reads per cell (-cr, --min_cell_reads) must be greater than 0.")
    if args.min_gene_counts is not None and args.min_gene_counts < 0:
        print("Minimum number of counts per gene (-gc, --min_gene_counts) must be greater than 0.")
    if args.n_top_genes is not None and args.n_top_genes < 0:
        print("Number of reads per cell (-t, --target_sum) must be greater than 0.")
    if args.target_sum is not None and args.target_sum < 0:
        print("Number of reads per cell (-t, --target_sum) must be greater than 0.")
    if args.nComp is not None and args.nComp < 0:
        print("Number of PCs (-n, --nComp) must be greater than 0.")

    # -------------------- set up output and start analysis --------------------
 
    if args.output is None:
        outdir = os.getcwd()
    else:
        outdir = args.output
    
    if outdir[-1] != "/":
        outdir += "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outf = open(outdir + "scpydr_results.txt", "w")
    
    sys.stdout.write("Welcome to scPyDR! Starting analysis.\n\n")

    # -------------------- load and preprocess data --------------------

    sys.stdout.write("Loading 10x Genomics data files...\n")

    """
    NOTES: 
    - metadata not used. should we have an option to call metadata into the function?
    - errors are very vague. try to make more user friendly/offer next steps to user.
    """

    try:
        adata = utils.load(args.datadir)
    except Exception:
        utils.ERROR(">> Failed to load 10x Genomics data into an AnnData object. Please check 1) path to directory \
                    and 2) directory contains properly formatted 10x Genomics files. \
                    Read more on 10x Genomics files here: https://www.10xgenomics.com/support/software/space-ranger/latest/advanced/hdf5-feature-barcode-matrix-format")
    try:
        adata = utils.preprocess(adata, args.min_genes, args.min_cells, args.min_cell_reads, \
                                    args.min_gene_counts, args.n_top_genes, args.target_sum)
    except Exception:
        utils.ERROR(">> Failed to preprocess data. Check arguments and try again.")
    try:
        df = utils.convert(adata)
    except Exception:
        utils.ERROR(">> Failed to convert AnnData object to DataFrame. This should not happen; please make a pull request to scPyDR.")

    sys.stdout.write("Data files successfully prepared! \n\n")    

    # -------------------- initialize, fit and transform data using pca --------------------

    sys.stdout.write("Running PCA for dimensionality reduction... \n")
    if nComp is None:
        #for n data points and m features, there are min(n-1,m) PCs.
        num_pts = df.shape[0] * df.shape[1]
        num_features = df.shape[1]
        nComp = min(num_pts - 1, num_features)
    pca = scpydrPCA(nComp=nComp) #create pca object
    sys.stdout.write("PCA object created. \n")
    pca.fit(df) #compute new PCs
    sys.stdout.write("New PCs computed. \n")
    df = pca.transform(df) #fit data to new PCs
    sys.stdout.write("Original data successfully projected onto the new PCs! \n\n")
    # TRY LATER: 
    # sys.stdout.write("Would you like to visualize the data in 2D? [y/n] \n")
    # input()

    # -------------------- analysis and conclusion --------------------

    sys.stdout.write("Summarizing results...")

    """ 
    TO DO: benchmarking, analysis, conclusions
    """
    sys.exit(0)

if __name__ == "__main__":
    main()
