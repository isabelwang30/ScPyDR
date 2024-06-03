import scanpy as sc
import argparse

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run UMAP using Scanpy on 10x Genomics scRNA-seq data")
    parser.add_argument("datadir", help="Directory containing 10x Genomics scRNA-seq data files.", metavar="DIR", type=str)
    args = parser.parse_args()

    # Load data
    adata = sc.read_10x_mtx(args.datadir)

    # Compute neighborhood graphs
    sc.pp.neighbors(adata) 

    # Cluster cells based on expression profiles
    sc.tl.leiden(adata)

    # Compute UMAP embedding
    sc.tl.umap(adata) 

    # Make a UMAP plot
    sc.pl.umap(adata)

if __name__ == "__main__":
    main()
