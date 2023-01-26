Link: 

This data package contains 3 datasets from the mouse primary visual cortex (Yao et al. 2021 Nature)
- rna: scRNA-seq 
- atac: snATAC-seq
- mc: snmC-seq

In total, there are 9 files in 3 groups. 

The enhancer- or gene-by-cell count matrices were located in the following .h5ad files (Scanpy/AnnData format): 
- counts_enh_atac.h5ad
- counts_enh_mc.h5ad
- counts_gene_rna.h5ad

The normalized gene-level profiles were located in the following .h5ad files: 
- profiles_hvgene_atac.h5ad
- profiles_hvgene_mc.h5ad
- profiles_hvgene_rna.h5ad

The annotation files, including the list of genes, enhancers, and their pairs within 1Mb, 
 were located in the following .tsv files:
- enhs_list.tsv
- genes_list.tsv
- enhancer_gene_pairs_1mbp.tsv