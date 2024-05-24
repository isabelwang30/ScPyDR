"""
utilities for scpydr
"""

# -------------------- set up --------------------
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from anndata import AnnData
import scanpy as sc

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
	prints error message and exits

	Parameters
	----------
	msg : str
	   Error message to print to terminal
	"""
	sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
	sys.exit(1)

# -------------------- preprocess data --------------------
"""filter, normalize, and log transform AnnData to keep highly variable genes"""
def preprocess(adata, min_genes=200, min_cells=5,
				min_cell_reads=None, min_gene_counts=None,
				n_top_genes=500, target_sum=1e4):
	"""
	Preprocess an AnnData object for downstream analysis.
	
	Parameters:
	- adata: AnnData object
	- min_genes: Min number of genes expressed per cell
	- min_cells: Min number of cells expressing a gene
	- min_cell_reads: Min number of reads per cell (optional)
	- min_gene_counts: Min number of counts per gene (optional)
	- n_top_genes: Number of highly variable genes to keep
	- target_sum: Number of reads per cell for normalization
	
	Returns:
	- Preprocessed AnnData object
	"""
	sc.pp.filter_cells(adata, min_genes=min_genes, inplace=True)
	sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)
	
	# additional filtering based on counts if specified
	if min_cell_reads is not None:
		sc.pp.filter_cells(adata, min_counts=min_cell_reads, inplace=True)
	if min_gene_counts is not None:
		sc.pp.filter_genes(adata, min_counts=min_gene_counts, inplace=True)

	# filter out cells with a high percentage of counts (>25%) from mitochondrial genes
	adata.var["mt"] = adata.var_names.str.startswith("MT-")
	sc.pp.calculate_qc_metrics(
		adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
	)
	adata = adata[adata.obs.pct_counts_mt <= 25, :]

	# normalize and log transform
	sc.pp.normalize_total(adata, target_sum=target_sum)
	sc.pp.log1p(adata)
	
	# identify and keep only highly variable genes
	sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
	adata = adata[:, adata.var.highly_variable]
	
	# scale data to have a mean of 0 and a variance of 1 for each gene, and limit extreme values to 10
	sc.pp.scale(adata, max_value=10)
	
	return adata

"""convert AnnData object to dataframe"""
def convert(adata, metadata_cols=None):
	"""
    Convert an AnnData object to a pandas DataFrame.
    
    Parameters:
    - adata: AnnData object to be converted.
    - metadata_columns: List of columns from adata.obs to include in the DataFrame. If None, only the expression matrix is included.
                        
    Returns:
    - A pandas DataFrame containing the expression matrix and optionally the metadata.
    """
    # convert the expression matrix (.X) to a df
	# if .X is a dense matrix
	if isinstance(adata.X, np.ndarray):  
		df = pd.DataFrame(adata.X, index=adata.obs.index, 
						  columns=adata.var.index)
	# if .X is a sparse matrix
	else:
		df = pd.DataFrame.sparse.from_spmatrix(adata.X,
											   index=adata.obs.index, columns=adata.var.index)
	
	if metadata_cols is not None:
        # select specified metadata columns from .obs
		df_metadata = adata.obs[metadata_cols]
        # concatenate metadata with the expression df along the columns
		df_combined = pd.concat([df_metadata, df], axis=1)
		return df_combined
	else:
		return df

# -------------------- pca class: initialize, fit and transform --------------------	
class scpydrPCA:
	
	"""constructor"""
	def __init__(self, nComp):
		self.nComp = nComp
		self.mean = None
		self.normalize = None
		self.components = None
		self.perc_explained_var = None # nComp PCs explain this amount of the total variance

	"""compute new PCs"""
	def fit(self, X):
		# standardize and center data
		self.mean = np.mean(X, axis=0)
		self.normalize = np.std(X, axis=0)
		X_std = (X - self.mean) / self.normalize

		# extract eigenvalues and eigenvectors through covariance matrix
		cov = np.cov(X.T)
		eigenvalues, eigenvectors = np.linalg.eig(cov)
		sort_idx = np.argsort(eigenvalues)[::-1]
		eigenvalues = eigenvalues[sort_idx]
		eigenvectors = eigenvectors[:, sort_idx] # column i is the i'th eigenvector
		self.components = eigenvectors[:self.nComp] # store subset of eigenvectors as the PCs of our data

		# explained variance ratio
		self.perc_explained_var = (np.sum(eigenvalues[:self.nComp])/np.sum(eigenvalues))*100 # for analysis later

		return self

	"""project data onto new pcs"""
	def transform(self, X):
		X_std = (X - self.mean / self.normalize) 
		return np.dot(X_std, self.components.T)