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