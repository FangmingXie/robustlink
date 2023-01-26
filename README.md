# Robust enhancer-gene association prediction using single cell transcriptomes and epigenomes


This repository hosts the package `robustlink`, a tool that integrates single-cell transcriptomes (scRNA-seq) and epigenomes (snATAC-seq and snmC-seq) and identify robust associations between cis-regulatory elements (enhancers) and genes. 

Reference:
- [Xie, Armand et al. 2021; Robust enhancer-gene regulation identified by single-cell transcriptomes and epigenomes](https://www.biorxiv.org/content/10.1101/2021.10.25.465795v1)

Correspondence: [Fangming Xie](mailto:f7xie@ucsd.edu) and [Eran A. Mukamel](mailto:emukamel@ucsd.edu)

# Getting started
### System requirements ###

This package is tested on a Ubuntu 18.04.6 LTS (Bionic Beaver) server. However, we expect it can be operated under a wide range of systems.
We recommend users to use a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) to install dependencies. This requires users to pre-install [Anaconda](https://www.anaconda.com/products/individual).

### Installation ###

Run the following command to clone this repo and install dependencies.
```bash
# clone this repo
git clone https://github.com/FangmingXie/robustlink.git
# create a conda env and install dependancies using the provided `env.yml` file.
cd ./robustlink
conda env create -f env.yml
conda activate env_robustlink 
pip install .
```

### Demo ###

The demo data, which includes mC, ATAC, and RNA profiles for >70,000 neurons from mouse primary motor cortex, can be downloaded using [this](https://drive.google.com/file/d/1FHjSn4MuNz7nxE7h_Ib8oiJcjm8XZaOs/view?usp=sharing) link.

Decompress `demodata.tar.gz` with the following command. For detailed description of data files, see `README_dataset.txt` after decompression. 
```bash
tar -zxvf data.tar.gz 
```

Now you can run through the entire enhancer-gene association analysis using a single command:
```
 ./link_mc_rna.sh & ./link_atac_rna.sh
```
This will generate a result folder `demo_res` that includes integrated datasets, metacells, and correlations between enhancers and genes for mC-RNA and ATAC-RNA, respectively. For speed, this demo only randomly samples 10% cells from each dataset. 

To visualize the results, then run through the `visualize_links.ipynb` jupyter notebook. This will readily generate visualizations with single-commands.
![](./doc/result_dist.png)

The whole demo takes about 5 minutes to run through.

### Prepare your data ###

You need to prepare your data as in the demo in `.h5ad` ([AnnData](https://anndata.readthedocs.io/en/latest/)) format. Specifically, for each dataset you need:
- `counts_yourdataset.h5ad`: a count matrix (cell-by-gene for RNA; cell-by-enhancer for mC and ATAC) of the single-cell transcriptome/epigenome data.
- `gene_profiles_yourdataaset.h5ad`: a gene-level feature matrix (average gene-level DNA methylation or ATAC signals for epigenome data). This information is not directly used for enhancer-gene association, but only to integrate cells from different datasets to identify cross-dataset metacells.

In addition, you need annotation file tables (.tsv):
- a gene list
- an enhancer list
- a list of enhancer-gene pairs to examine (all pairs within ~1Mbp)

### API ###
