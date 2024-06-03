"""
Utilities for scpydr
"""

# -------------------- set up --------------------
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from anndata import AnnData
import scanpy as sc
import umap
import os
import leidenalg

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# -------------------- error handling --------------------
def ERROR(msg):
    """
    Prints error message and exits.

    Parameters
    ----------
    msg : str
       Error message to print to terminal.
    """
    sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg))
    sys.exit(1)

# -------------------- load and preprocess data --------------------
import shutil

def load(datadir, prefix="", cache=True):
    # List the files in the directory
    files = os.listdir(datadir)
    
    # Debug: print out the found files
    print("Files in directory:", files)
    
    # Find the files with optional prefix
    barcodes_file = next((f for f in files if 'barcodes' in f and f.endswith('.tsv.gz')), None)
    features_file = next((f for f in files if 'features' in f and f.endswith('.tsv.gz')), None)
    matrix_file = next((f for f in files if 'matrix' in f and f.endswith('.mtx.gz')), None)
    
    # Debug: print the identified files
    print("Identified barcodes file:", barcodes_file)
    print("Identified features file:", features_file)
    print("Identified matrix file:", matrix_file)
    
    # Check if all necessary files are found
    if not all([barcodes_file, features_file, matrix_file]):
        ERROR("Missing required files in the directory. Ensure 'barcodes', 'features', and 'matrix' files are present.")
    
    barcodes_path = os.path.join(datadir, barcodes_file)
    features_path = os.path.join(datadir, features_file)
    matrix_path = os.path.join(datadir, matrix_file)
    
    # Create temporary copies of the files with the expected names
    temp_dir = os.path.join(datadir, "temp_10x_files")
    
    try:
        os.makedirs(temp_dir, exist_ok=True)
        
        temp_barcodes_path = os.path.join(temp_dir, "barcodes.tsv.gz")
        temp_features_path = os.path.join(temp_dir, "features.tsv.gz")
        temp_matrix_path = os.path.join(temp_dir, "matrix.mtx.gz")
        
        shutil.copy(barcodes_path, temp_barcodes_path)
        shutil.copy(features_path, temp_features_path)
        shutil.copy(matrix_path, temp_matrix_path)
        
        # Debug: print the temporary paths
        print("Temporary barcodes path:", temp_barcodes_path)
        print("Temporary features path:", temp_features_path)
        print("Temporary matrix path:", temp_matrix_path)
    
    except Exception as e:
        ERROR(f"An error occurred during file copying: {e}")
    
    # Read the data using scanpy's read_10x_mtx function
    try:
        adata = sc.read_10x_mtx(
            temp_dir,
            var_names='gene_symbols' if 'features.tsv.gz' in features_file else 'gene_ids',
            cache=cache
        )
        print("Data loaded successfully.")
    except Exception as e:
        ERROR(f"An error occurred while reading the data: {e}")
    finally:
        # Clean up the temporary directory
        shutil.rmtree(temp_dir)
    
    return adata


def preprocess(adata, min_genes=200, min_cells=5,
                min_cell_reads=None, min_gene_counts=None,
                n_top_genes=500, target_sum=1e4):
    """
    Preprocess an AnnData object for downstream analysis. Filters, normalizes, and log transforms the data, then keeps only the highly variable genes.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to preprocess.
    min_genes : int, optional
        Minimum number of genes expressed per cell (default is 200).
    min_cells : int, optional
        Minimum number of cells expressing a gene (default is 5).
    min_cell_reads : int, optional
        Minimum number of reads per cell (default is None).
    min_gene_counts : int, optional
        Minimum number of counts per gene (default is None).
    n_top_genes : int, optional
        Number of highly variable genes to keep (default is 500).
    target_sum : float, optional
        Number of reads per cell for normalization (default is 1e4).

    Returns
    -------
    AnnData
        Preprocessed AnnData object.
    """
    sc.pp.filter_cells(adata, min_genes=min_genes, inplace=True)
    sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)
    
    # Additional filtering based on counts if specified
    if min_cell_reads is not None:
        sc.pp.filter_cells(adata, min_counts=min_cell_reads, inplace=True)
    if min_gene_counts is not None:
        sc.pp.filter_genes(adata, min_counts=min_gene_counts, inplace=True)

    # Filter out cells with a high percentage of counts (>25%) from mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata = adata[adata.obs.pct_counts_mt <= 25, :]

    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
    # Identify and keep only highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data to have a mean of 0 and a variance of 1 for each gene, and limit extreme values to 10
    sc.pp.scale(adata, max_value=10)
    
    return adata

def convert(adata, metadata_cols=None):
    """
    Convert an AnnData object to a pandas DataFrame.

    Parameters
    ----------
    adata : AnnData
        AnnData object to be converted.
    metadata_cols : list of str, optional
        List of columns from adata.obs to include in the DataFrame. If None, only the expression matrix is included.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the expression matrix and optionally the metadata.
    """
    # Convert the expression matrix (.X) to a df
    if isinstance(adata.X, np.ndarray):  
        df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
    else:
        df = pd.DataFrame.sparse.from_spmatrix(adata.X, index=adata.obs.index, columns=adata.var.index)
    
    if metadata_cols is not None:
        # Select specified metadata columns from .obs
        df_metadata = adata.obs[metadata_cols]
        # Concatenate metadata with the expression df along the columns
        df_combined = pd.concat([df_metadata, df], axis=1)
        return df_combined
    else:
        return df

# -------------------- pca class: initialize, fit and transform --------------------    
class scpydrPCA:
    """PCA class for single-cell RNA-seq data analysis."""

    def __init__(self, nComp):
        """
        Constructor for the PCA class.

        Parameters
        ----------
        nComp : int
            Number of principal components to compute.
        """
        self.nComp = nComp
        self.mean = None
        self.normalize = None
        self.components = None
        self.perc_explained_var = None  # nComp PCs explain this amount of the total variance

    def fit(self, X):
        """
        Compute new principal components.

        Parameters
        ----------
        X : np.ndarray
            Data matrix to compute principal components for.

        Returns
        -------
        scpydrPCA
            The PCA object with fitted principal component axes.
        """
        # Standardize and center data
        self.mean = np.mean(X, axis=0)
        self.normalize = np.std(X, axis=0)
        X_std = (X - self.mean) / self.normalize

        # Extract eigenvalues and eigenvectors through covariance matrix
        cov = np.cov(X.T)
        eigenvalues, eigenvectors = np.linalg.eig(cov)
        sort_idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[sort_idx]
        eigenvectors = eigenvectors[:, sort_idx]  # Column i is the i'th eigenvector
        self.components = eigenvectors[:self.nComp]  # Store subset of eigenvectors as the PCs of our data
        # Explained variance ratio
        self.perc_explained_var = (np.sum(eigenvalues[:self.nComp]) / np.sum(eigenvalues)) * 100  # For analysis later

        return self

    def transform(self, X):
        """
        Project data onto new principal components.

        Parameters
        ----------
        X : np.ndarray
            Data matrix with raw counts.

        Returns
        -------
        np.ndarray
            Transformed data matrix made by projecting raw counts onto the new principal component axes.
        """
        X_std = (X - self.mean) / self.normalize 
        return np.dot(X_std, self.components.T)
    
def save_pca_results(outdir, filename_prefix, pca_results):
    """
    Save PCA results to a file.

    Parameters
    ----------
    outdir : str
        Output directory to save the file.
    filename_prefix : str
        Prefix for the output filename.
    pca_results : np.ndarray
        PCA results to save.

    Returns
    -------
    None
    """
    output_file = os.path.join(outdir, f"{filename_prefix}_pca.txt")
    np.savetxt(output_file, pca_results, delimiter="\t")
    print(f"PCA results saved to {output_file}\n")

def plot_pca_results(outdir, filename_prefix, pca_results):
    """
    Plot PCA results and save the plot to a file.

    Parameters
    ----------
    outdir : str
        Output directory to save the plot.
    filename_prefix : str
        Prefix for the output filename.
    pca_results : np.ndarray
        PCA results to plot.

    Returns
    -------
    None
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_results[:, 0], pca_results[:, 1], s=20, c='b', alpha=0.5)
    plt.title('PCA Plot', fontsize=16)
    plt.xlabel('PC1', fontsize=14)
    plt.ylabel('PC2', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    output_plot = os.path.join(outdir, f"{filename_prefix}_pca_plot.png")
    plt.savefig(output_plot)
    plt.close()
    print(f"PCA plot saved to {output_plot}\n")

def umap_embedding(adata, min_dist=0.1, n_components=2, n_epochs=200, learning_rate=1.0, n_neighbors=None):
    """
    Compute UMAP embedding for visualization.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing the input data.
    min_dist : float, optional
        Minimum distance between points in the UMAP embedding (default is 0.1).
    n_components : int, optional
        Number of dimensions of the UMAP embedding (default is 2).
    n_epochs : int, optional
        Number of epochs for optimizing the UMAP embedding (default is 200).
    learning_rate : float, optional
        Learning rate for optimizing the UMAP embedding (default is 1.0).
    n_neighbors : int or None, optional
        Number of nearest neighbors to use for constructing the UMAP graph. If None, it will be determined automatically based on the size of the data.

    Returns
    -------
    np.ndarray
        UMAP embedding of the input data.
    """
    adatac = adata.copy()  # Make a copy to avoid modifying the original object
    
    # Compute default number of neighbors if not specified by the user
    if n_neighbors is None:
        # Use a heuristic based on the size of the data: >10000 genes has different number of neighbors
        n_neighbors = 15 if adatac.shape[0] > 10000 else 10
    
    # Compute nearest neighbors using scanpy's algorithm
    sc.pp.neighbors(adatac, n_neighbors=n_neighbors)
    
    # Create UMAP object
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        metric='euclidean',
        random_state=None,  # You can set a random state if desired
        n_components=n_components,
        min_dist=min_dist,
        n_epochs=n_epochs,
        learning_rate=learning_rate
    )
    
    # Fit and transform the data
    embedding = reducer.fit_transform(adata.X)
    return embedding

def plot_umap_results(outdir, filename_prefix, umap_embedding):
    """
    Plot UMAP results and save the plot to a file.

    Parameters
    ----------
    outdir : str
        Output directory to save the plot.
    filename_prefix : str
        Prefix for the output filename.
    umap_embedding : np.ndarray
        UMAP results to plot.

    Returns
    -------
    None
    """
    fig, ax = plt.subplots(figsize=(8, 6))  # Set the figure size
    ax.scatter(umap_embedding[:, 0], umap_embedding[:, 1], s=20, c='b', alpha=0.5)  # Adjust marker size, color, and transparency
    ax.set_title('UMAP Embedding', fontsize=16)  # Set title and adjust font size
    ax.set_xlabel('UMAP 1', fontsize=14)  # Set x-axis label and adjust font size
    ax.set_ylabel('UMAP 2', fontsize=14)  # Set y-axis label and adjust font size
    ax.grid(True, linestyle='--', alpha=0.5)  # Add grid with dashed lines and transparency
    plt.tight_layout()  # Adjust layout to prevent overlapping labels
    
    output_plot = os.path.join(outdir, f"{filename_prefix}_umap_plot.png")
    plt.savefig(output_plot)
    plt.close()
    print(f"UMAP plot saved to {output_plot}\n")
