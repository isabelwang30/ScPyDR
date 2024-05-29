# scPyDR
### Single-Cell Python Dimensionality Reduction
`scPyDR` is a Python package containing tools for the dimensionality reduction and visualization of single-cell RNA sequencing data. The three tools are simpler versions of Scanpy's [`scanpy.pp.pca`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.pca.html), [`scanpy.tl.tsne`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.tsne.html), and [`scanpy.tl.umap`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html).

[Prerequisites](#Prerequisites) | [Installation](#Installation) | [Basic Usage](#Usage) | [scPyDR Options](#Options) | [File Formats](#Formats) | 
[Contributors](#Contributors)

## Prerequisites<a name="Prerequisites"></a>
`scPyDR` requires the following python libraries to be installed:
* numpy
* pandas
* matplotlib
* anndata
* scanpy
* umap

These can be installed with the following `pip` command:
```
pip install numpy pandas matplotlib anndata scanpy umap
```

*Note: if you do not have root access, the packages can be installed locally with the following command:*
```
pip install --user numpy pandas matplotlib anndata scanpy umap
```

## Installation<a name="Installation"></a>
Once the required libraries are installed, `scPyDR` can be installed with the following commands:
```
git clone https://github.com/isabelwang30/scPyDR.git
cd scPyDR
python setup.py install
```

*Note: if you do not have root access, `scPyDR` can be installed locally with the following commands:*
```
git clone https://github.com/isabelwang30/scPyDR.git
cd scPyDR
python setup.py install --user
```

If the install was successful, the command `scpydr --help` should show a useful help message.

*Note: if you get an error that says the `scpydr` command was not found, you may need to include the script installation path in your `$PATH` variable before calling `scpydr`. You can do this with the following command:*
```
export PATH=$PATH:/home/<user>/.local/bin
```

## Basic Usage<a name="Usage"></a>


## scPyDR Options<a name="Options"></a>


## File Formats<a name="Formats"></a>


## Contributors<a name="Contributors"></a>

