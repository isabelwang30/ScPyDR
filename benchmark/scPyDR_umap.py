import scanpy as sc
import argparse
from scPyDR.utils import umap_embedding, plot_umap_results

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run UMAP using scPyDR on 10x Genomics scRNA-seq data")
    parser.add_argument("datadir", help="Directory containing 10x Genomics scRNA-seq data files.", metavar="DIR", type=str)
    parser.add_argument("-o", "--outdir", help="Output directory to store results.", type=str, default="scpydr_results")
    args = parser.parse_args()

    # Load data
    adata = sc.read_10x_mtx(args.datadir)

    # Call scPyDR UMAP functions
    umap_embedding_result, cluster_labels = umap_embedding(adata)
    plot_umap_results(args.outdir, "data", umap_embedding_result, cluster_labels)

if __name__ == "__main__":
    main()