import os
import sys
import shutil
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import umap.umap_ as umap
from anndata import AnnData
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

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

def ERROR(msg):
    """
    Prints error message and exits.
    """
    sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + f"{msg}\n")
    sys.exit(1)

def load(datadir, prefix="", cache=True):
    files = os.listdir(datadir)
    barcodes_file = next((f for f in files if 'barcodes' in f and f.endswith('.tsv.gz')), None)
    features_file = next((f for f in files if 'features' in f and f.endswith('.tsv.gz')), None)
    matrix_file = next((f for f in files if 'matrix' in f and f.endswith('.mtx.gz')), None)
    
    if not all([barcodes_file, features_file, matrix_file]):
        ERROR("Missing required files in the directory. Ensure 'barcodes', 'features', and 'matrix' files are present.")
    
    barcodes_path = os.path.join(datadir, barcodes_file)
    features_path = os.path.join(datadir, features_file)
    matrix_path = os.path.join(datadir, matrix_file)
    
    temp_dir = os.path.join(datadir, "temp_10x_files")
    
    try:
        os.makedirs(temp_dir, exist_ok=True)
        
        temp_barcodes_path = os.path.join(temp_dir, "barcodes.tsv.gz")
        temp_features_path = os.path.join(temp_dir, "features.tsv.gz")
        temp_matrix_path = os.path.join(temp_dir, "matrix.mtx.gz")
        
        shutil.copy(barcodes_path, temp_barcodes_path)
        shutil.copy(features_path, temp_features_path)
        shutil.copy(matrix_path, temp_matrix_path)
    
    except Exception as e:
        ERROR(f"An error occurred during file copying: {e}")
    
    try:
        adata = sc.read_10x_mtx(
            temp_dir,
            var_names='gene_symbols' if 'features.tsv.gz' in features_file else 'gene_ids',
            cache=cache
        )
    except Exception as e:
        ERROR(f"An error occurred while reading the data: {e}")
    finally:
        shutil.rmtree(temp_dir)
    
    return adata

def preprocess(adata, min_genes=200, min_cells=5,
               min_cell_reads=None, min_gene_counts=None,
               n_top_genes=500, target_sum=1e4):
    adatac = adata.copy()
    sc.pp.filter_cells(adatac, min_genes=min_genes, inplace=True)
    sc.pp.filter_genes(adatac, min_cells=min_cells, inplace=True)
    
    if min_cell_reads is not None:
        sc.pp.filter_cells(adatac, min_counts=min_cell_reads, inplace=True)
    if min_gene_counts is not None:
        sc.pp.filter_genes(adatac, min_counts=min_gene_counts, inplace=True)

    adatac.var["mt"] = adatac.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adatac, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adatac = adatac[adatac.obs.pct_counts_mt <= 25, :].copy()

    sc.pp.normalize_total(adatac, target_sum=target_sum)
    sc.pp.log1p(adatac)

    # marker genes to keep
    # Enterocytes = Slc10a2, Alpi, Krt20, Klf5, Ccl25
    # Paneth cells = Atg16l1, Dll4, Il4ra, Tlr9, Ctnnb1, Tm4sf20
    # Goblet cells = Muc2, Muc5ac, Clca1, Krt7, Muc13
    # Enteroendocrine cells = Chga, Gfra3
    # T cells = TRBC2, CD3D, CD3G, CD3E, LTB, IL7R
    # B cells = Pxk, Ms4a1, Cd19, Cd74, Cd79a, Cd79b
    # Macrophages = Cd68, Fcgr1, Naaa, Lyz2, Ccl12, Sepp1, March1
    # Fibroblasts = Vim, Pdgfrb, Lum, Col6a2, Vtn, Mfap5, Col1a2
    # Endothelial cells = Cd93, Vwf, Emcn, Egfl7, Flt1, Id3, Kdr

    genes = ['Slc10a2', 'Alpi', 'Krt20', 'Klf5', 'Ccl25',
             'Atg16l1', 'Dll4', 'Il4ra', 'Tlr9', 'Ctnnb1', 'Tm4sf20',
             'Muc2', 'Clca1', 'Krt7', 'Muc13',
             'Chga', 'Gfra3',
             'Trbc2', 'Cd3d', 'Cd3g', 'Cd3e', 'Ltb', 'Il7r',
             'Pxk', 'Ms4a1', 'Cd19', 'Cd74', 'Cd79a', 'Cd79b',
             'Cd68', 'Fcgr1', 'Naaa', 'Lyz2', 'Ccl12', 'March1',
             'Vim', 'Pdgfrb', 'Lum', 'Col6a2', 'Vtn', 'Mfap5', 'Col1a2',
             'Cd93', 'Vwf', 'Emcn', 'Egfl7', 'Flt1', 'Id3', 'Kdr']
    
    # Check if genes are present in adatac
    for gene in genes:
        if gene not in adatac.var_names:
            print(f"Gene {gene} not found in the dataset.")

    # filter top highly variable genes
    sc.pp.highly_variable_genes(adatac, n_top_genes=n_top_genes)

    # keep highly variable genes and marker genes
    adata_var = adatac[:, adatac.var.index.isin(genes) | adatac.var.highly_variable].copy()
    sc.pp.scale(adata_var, max_value=10)
    
    return adata_var

def umap_embedding(adata, min_dist=0.1, n_components=2, n_epochs=200, learning_rate=1.0, n_neighbors=None):
    adatac = adata.copy()
    
    if n_neighbors is None:
        n_neighbors = 15 if adatac.shape[0] > 10000 else 10
    
    sc.pp.pca(adatac, n_comps=50)
    sc.pp.neighbors(adatac, n_neighbors=n_neighbors, use_rep='X_pca')
    sc.tl.leiden(adatac, flavor="igraph", n_iterations=2, directed=False)

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        metric='euclidean',
        random_state=None,
        n_components=n_components,
        min_dist=min_dist,
        n_epochs=n_epochs,
        learning_rate=learning_rate
    )
    
    embedding = reducer.fit_transform(adatac.obsm['X_pca'])
    cluster_labels = adatac.obs['leiden']
    
    return embedding, cluster_labels, adatac

def plot_annotated_results(outdir, filename_prefix, umap_embedding, cluster_labels, adata, genes_to_plot, gene_to_cell_type):
    cluster_labels = cluster_labels.values.astype(int)
    num_clusters = len(np.unique(cluster_labels))
    cmap = ListedColormap(plt.cm.get_cmap('viridis', num_clusters).colors)
    
    num_genes = len(genes_to_plot)
    fig, axes = plt.subplots(1, num_genes + 1, figsize=(8 * (num_genes + 1), 6))
    
    # Plot by cluster labels
    scatter1 = axes[0].scatter(umap_embedding[:, 0], umap_embedding[:, 1], s=20, c=cluster_labels, cmap=cmap, alpha=0.5)
    axes[0].set_title('UMAP Embedding by Clusters', fontsize=16)
    axes[0].set_xlabel('UMAP 1', fontsize=14)
    axes[0].set_ylabel('UMAP 2', fontsize=14)
    axes[0].grid(True, linestyle='--', alpha=0.5)
    
    # Create legend for clusters
    handles1 = []
    for cluster in np.unique(cluster_labels):
        handles1.append(axes[0].scatter([], [], s=50, c=[cmap(cluster)], alpha=0.5, label=cluster))
    axes[0].legend(handles=handles1, title="Clusters", bbox_to_anchor=(1.01, 1.02), loc='upper left')
    
    # Plot expression levels of each gene
    for i, gene in enumerate(genes_to_plot):
        if gene in adata.var_names:
            expression = adata[:, gene].X.toarray().flatten()
        else:
            print(f"Gene {gene} not found in the dataset.")
            continue
        
        # Choose a colormap for the gene expression levels
        cmap_gene = plt.cm.get_cmap('Reds')  # You can change 'Reds' to 'Blues', 'Greens', etc.

        sc = axes[i + 1].scatter(umap_embedding[:, 0], umap_embedding[:, 1], s=20, c=expression, cmap=cmap_gene, alpha=0.5)
        if gene in gene_to_cell_type and gene_to_cell_type[gene]:
            cell_type = gene_to_cell_type[gene]
            axes[i + 1].set_title(f'UMAP Embedding by {gene} ({cell_type}) Expression', fontsize=16)
        else:
            axes[i + 1].set_title(f'UMAP Embedding by {gene} Expression', fontsize=16)
        axes[i + 1].set_xlabel('UMAP 1', fontsize=14)
        axes[i + 1].set_ylabel('UMAP 2', fontsize=14)
        axes[i + 1].grid(True, linestyle='--', alpha=0.5)
        fig.colorbar(sc, ax=axes[i + 1], fraction=0.046, pad=0.04, label=f'{gene} expression')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_plot = os.path.join(outdir, f"{filename_prefix}_umap_expression.png")
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()
    print(f"Annotated UMAP plot saved to {output_plot}\n")

def main():
    parser = argparse.ArgumentParser(
        prog="run_umap",
        description="UMAP Dimensionality Reduction for 10x Genomics scRNA-seq Data"
    )

    parser.add_argument("datadir", help="Directory containing 10x Genomics scRNA-seq data files.", metavar="DIR", type=str)
    parser.add_argument("-o", "--output", help="Output directory to store results. Default: working directory.", type=str, required=False)
    parser.add_argument("-g", "--min_genes", help="Minimum number of genes expressed per cell. Default: 200.", type=int, default=200, required=False)
    parser.add_argument("-c", "--min_cells", help="Minimum number of cells expressing a gene. Default: 5.", type=int, default=5, required=False)
    parser.add_argument("-cr", "--min_cell_reads", help="Minimum number of reads per cell. Default: None.", type=int, required=False)
    parser.add_argument("-gc", "--min_gene_counts", help="Minimum number of counts per gene. Default: None.", type=int, required=False)
    parser.add_argument("-ntop", "--n_top_genes", help="Number of highly variable genes to keep. Default: 500.", type=int, default=500, required=False)
    parser.add_argument("-t", "--target_sum", help="Number of reads per cell for normalization. Default: 1e4.", type=float, default=1e4, required=False)
    parser.add_argument("--min_dist", help="Minimum distance between points in UMAP embedding. Default: 0.1.", type=float, default=0.1, required=False)
    parser.add_argument("--n_components", help="Number of dimensions of UMAP embedding. Default: 2.", type=int, default=2, required=False)
    parser.add_argument("--n_epochs", help="Number of epochs for optimizing UMAP embedding. Default: 200.", type=int, default=200, required=False)
    parser.add_argument("--learning_rate", help="Learning rate for optimizing UMAP embedding. Default: 1.0.", type=float, default=1.0, required=False)
    parser.add_argument("--n_neighbors", help="Number of nearest neighbors to use for constructing UMAP graph. Default: determined automatically.", type=int, required=False)

    args = parser.parse_args()

    if not os.path.isdir(args.datadir):
        ERROR("Invalid directory path to 10x Genomics data files.")
    if args.output and not os.path.isdir(args.output):
        ERROR("Invalid output directory path.")

    datadir = args.datadir
    outdir = args.output if args.output else os.getcwd()
    filename_prefix = os.path.basename(os.path.normpath(datadir))

    print(f"Loading data from {datadir}")
    adata = load(datadir)

    print("Preprocessing data...")
    adata = preprocess(adata, min_genes=args.min_genes, min_cells=args.min_cells,
                       min_cell_reads=args.min_cell_reads, min_gene_counts=args.min_gene_counts,
                       n_top_genes=args.n_top_genes, target_sum=args.target_sum)

    print("Running UMAP...")
    embedding, cluster_labels, adatac = umap_embedding(adata)

    # in mice:
    # Enterocytes = Slc10a2, Alpi, Krt20, Klf5, Ccl25
    # Paneth cells = Atg16l1, Dll4, Il4ra, Tlr9, Ctnnb1, Tm4sf20
    # Goblet cells = Muc2, Muc5ac, Clca1, Krt7, Muc13
    # Enteroendocrine cells = Chga, Gfra3

    # T cells = TRBC2, CD3D, CD3G, CD3E, LTB, IL7R
    # B cells = Pxk, Ms4a1, Cd19, Cd74, Cd79a, Cd79b
    # Macrophages = Cd68, Fcgr1, Naaa, Lyz2, Ccl12, Sepp1, March1

    # Fibroblasts = Vim, Pdgfrb, Lum, Col6a2, Vtn, Mfap5, Col1a2

    # Endothelial cells = Cd93, Vwf, Emcn, Egfl7, Flt1, Id3, Kdr
    
    genes_to_plot = ['Klf5',
             'Tm4sf20',
             'Muc13',
             'Gfra3',
             'Trbc2',
             'Cd79a',
             'Lyz2',
             'Col6a2',
             'Flt1']
    
    # Define the mapping of genes to cell types
    gene_to_cell_type = {
        'Klf5': 'Enterocytes',
        'Tm4sf20': 'Paneth cells',
        'Muc13': 'Goblet cells',
        'Gfra3': 'Enteroendocrine cells',
        'Trbc2': 'T cells',
        'Cd79a': 'B cells',
        'Lyz2': 'Macrophages',
        'Col6a2': 'Fibroblasts',
        'Flt1': 'Endothelial cells'
    }
    print("Plotting annotated UMAP results...")
    plot_annotated_results(outdir, filename_prefix, embedding, cluster_labels, adatac, genes_to_plot, gene_to_cell_type)

    print("UMAP dimensionality reduction completed successfully.")

if __name__ == "__main__":
    main()

"""
def plot_annotated_results(outdir, filename_prefix, umap_embedding, cluster_labels, cell_types):
    # Convert cluster labels to integer type if not already
    cluster_labels = cluster_labels.values.astype(int)
    num_clusters = len(np.unique(cluster_labels))
    cmap = ListedColormap(plt.cm.get_cmap('viridis', num_clusters).colors)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot by cluster labels
    scatter1 = axes[0].scatter(umap_embedding[:, 0], umap_embedding[:, 1], s=20, c=cluster_labels, cmap=cmap, alpha=0.5)
    axes[0].set_title('UMAP Embedding by Clusters', fontsize=16)
    axes[0].set_xlabel('UMAP 1', fontsize=14)
    axes[0].set_ylabel('UMAP 2', fontsize=14)
    axes[0].grid(True, linestyle='--', alpha=0.5)
    
    # Create legend for clusters
    handles1 = []
    for cluster in np.unique(cluster_labels):
        handles1.append(axes[0].scatter([], [], s=50, c=[cmap(cluster)], alpha=0.5, label=cluster))
    axes[0].legend(handles=handles1, title="Clusters", bbox_to_anchor=(1.01, 1.02), loc='upper left')
    
    # Plot by cell types
    unique_cell_types = np.unique(cell_types)
    num_cell_types = len(unique_cell_types)
    cell_type_cmap = ListedColormap(plt.cm.get_cmap('tab20', num_cell_types).colors)
    cell_type_dict = {cell_type: i for i, cell_type in enumerate(unique_cell_types)}
    cell_type_indices = np.vectorize(cell_type_dict.get)(cell_types)
    
    scatter2 = axes[1].scatter(umap_embedding[:, 0], umap_embedding[:, 1], s=20, c=cell_type_indices, cmap=cell_type_cmap, alpha=0.5)
    axes[1].set_title('UMAP Embedding by Cell Types', fontsize=16)
    axes[1].set_xlabel('UMAP 1', fontsize=14)
    axes[1].set_ylabel('UMAP 2', fontsize=14)
    axes[1].grid(True, linestyle='--', alpha=0.5)
    
    # Create legend for cell types
    handles2 = []
    for cell_type in unique_cell_types:
        handles2.append(axes[1].scatter([], [], s=50, c=[cell_type_cmap(cell_type_dict[cell_type])], alpha=0.5, label=cell_type))
    axes[1].legend(handles=handles2, title="Cell Types", bbox_to_anchor=(1.01, 1.02), loc='upper left')
    
    # Adjust layout
    fig.subplots_adjust(right=0.85)
    plt.tight_layout()
    
    # Save the plot
    output_plot = os.path.join(outdir, f"{filename_prefix}_umap_annotated.png")
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()
    print(f"Annotated UMAP plot saved to {output_plot}\n")
"""
# in main
"""
    # Identify marker genes that distinguish the clusters of cells using wilcoxon method to better handle outliers
    sc.tl.rank_genes_groups(adatac, groupby='leiden', method='wilcoxon')

    # Get the marker genes for each cluster
    marker_genes = adatac.uns['rank_genes_groups']['names']

    # Display top 5 marker genes for each cluster
    for cluster in marker_genes.dtype.names:
        print(f"Cluster {cluster}: {marker_genes[cluster][:5]}")


    
    
    common marker genes- PangloaDB mouse 
Dcn (4) = fibroblasts, endocrine cells (atlas)

Tnc (3)- fibroblasts, osteoblasts, astrocytes, Scwann cells- mostly glial/cns related
Col3a1 (3)- fibroblasts, pericytes, schwann cells

* Agr2 (3)- paneth cells and enterocyte (both small intestine)
* Dmbt1 (3)- paneth cells and enterocyte (both small intestine)

Rgs1 (3)- macrophages
Il2rb (3)- T cells
Tpm2 (3) - fibroblasts
Rgs5 (3) - pericytes

Cd74 (2)- macrophages
Iglc2 (2)
Ces1d (2)
Meg3 (2)
Edn1 (2)
Hopx (2)
Igha (2)
Jchain (2)
Guca2a (2)
Cd28 (2)
Ramp3 (2)
Trbc1 (2)
Ano1 (2)
Myh11 (2)
Fcer1g (2)
Fabp4 (2)

    # Manually annotate cell types
    cell_type_annotations = {
    '0': 'Macrophages',
    '1': 'Fibroblasts',
    '2': 'Fibroblasts',
    '3': 'Enterocytes and Panneths',
    '4': 'B cells',
    '5': 'Enterocytes and Panneths',
    '6': 'T cells',
    '7': 'T cells',
    '8': 'Enterocytes',
    '9': 'B cells',
    '10': 'Fibroblasts',
    '11': 'Enterocytes and Panneths',
    '12': '',
    '13': '',
    '14': '',
    '15': '',
    '16': '',
    '17': '',
    '18': '',
    '19': '',
    '20': '',
    '21': '',
    '22': '',
    '23': '',
    '24': '',
    '25': '',
    '26': '',
    '27': ''
    }
    adatac.obs['cell_type'] = adatac.obs['leiden'].map(cell_type_annotations).astype('category')

    print("Plotting annotated UMAP results...")
    plot_annotated_results(outdir, filename_prefix, embedding, cluster_labels, cell_type_annotations)
"""