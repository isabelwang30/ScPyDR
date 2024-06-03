import scanpy as sc
import argparse


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run Scanpy PCA on 10x Genomics scRNA-seq data")
    parser.add_argument("datadir", help="Directory containing 10x Genomics scRNA-seq data files.", metavar="DIR", type=str)
    parser.add_argument("-n", "--nComp", help="Number of principal components.", type=int)
    args = parser.parse_args()
    

    # Load data
    adata = sc.read_10x_mtx(args.datadir)
    
    # Run PCA
    sc.pp.pca(adata)
    
    # Calculate percent variance
    explained_var_ratio = adata.uns['pca']['variance_ratio']
    explained_var_ratio_n = explained_var_ratio[:args.nComp].sum()
    explained_var_perc_n = round(explained_var_ratio_n*100, 2)
    print(f"Total explained variance for the first {args.nComp} components: {explained_var_perc_n}%")
    

if __name__ == "__main__":
    main()
