import scanpy as sc
import argparse

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run Scanpy PCA on 10x Genomics scRNA-seq data")
    parser.add_argument("datadir", help="Directory containing 10x Genomics scRNA-seq data files.", metavar="DIR", type=str)
    parser.add_argument("-o", "--output", help="Output directory to store results.", type=str, default="scanpy_results")
    args = parser.parse_args()

    # Load data
    adata = sc.read_10x_mtx(args.datadir)

    # Run PCA with 10 PCs
    sc.pp.pca(adata, n_comps=10)
    # Plot top two PCs
    sc.pl.pca(adata, save=True, show=False)

    # Save results
    adata.write(args.output)

if __name__ == "__main__":
    main()