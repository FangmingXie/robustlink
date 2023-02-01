#!/bin/bash
# 1. run scFusion with different downsampling
# output: knn_xy matrices, and cells (i)

# 2. within-modal clustering ...
# for each cell set from (i) -> run Leiden clustering with different resolution (r)
# output: clustering_results.tsv.gz (i, r): (i) separate tables, (r) columns per table

# 3. correlation analysis ...
# for each combination of (i, r), run a correlation analysis
# output: corrs.pkl (i, r) 

## configurations
# input and output
out_dir="./demoresults"
data_dir="./demodata"
tolink="enhancer_gene_pairs_1mbp.tsv"
countdata_gene="counts_gene_rna.h5ad"
countdata_enh="counts_enh_mc.h5ad"
fusiondata_rna="profiles_hvgene_rna.h5ad"
fusiondata_mc="profiles_hvgene_mc.h5ad" 

# parameters
ka=30  # within modality k nearest neighbors
knn=30 # across modality k nearest neighbors
corr_type='spearmanr'
subsample_frac=0.1 # take 10% of all the data -- faster for demo
subsample_times=1
resolutions=(10) # Leiden clustering resolutions used to generate metacells -- just 1 for demo
num_metacell_limit=1001
## end of configuration

# prepare
study_tag="link_rna_mc_ka${ka}_knn${knn}"
echo $study_tag
scfusion_dir=${out_dir} # results of scFusion

# 1.
# # run SCF - repeatedly on subsets of cells (i) 
echo "STEP1..."
python -m robustlink scfusion \
	-i  ${data_dir} \
	-id ${fusiondata_mc} ${fusiondata_rna} \
	-fd ${fusiondata_mc} \
	-im "mc" "rna" \
	-o  ${out_dir} \
	-tag ${study_tag} \
	--ka_smooth $ka \
	--knn $knn \
	-s  ${subsample_frac} \
	-sn ${subsample_times}

# 2.
# run leiden clustering for each (i) 
echo "STEP2..."
python -m robustlink metacell \
	-i  "${data_dir}/${fusiondata_rna}" \
	-o  ${out_dir} \
	-tag ${study_tag} \
	-sn ${subsample_times} \
	-r  ${resolutions}

# # 3.
# # correlation analysis (i, r) - r is for resolution
echo "STEP3..."
for (( i=0; i<${subsample_times}; i++ )); do
	python -m robustlink corr_mc \
		--tolink         "${data_dir}/$tolink" \
		--countdata_gene "${data_dir}/${countdata_gene}" \
		--countdata_enh  "${data_dir}/${countdata_enh}" \
		--scfusion_dir   ${scfusion_dir} \
		--fusiondata_rna ${fusiondata_rna} \
		--fusiondata_mc  ${fusiondata_mc} \
		-tag ${study_tag} \
		-isub $i \
		--corr_type ${corr_type} \
		-o ${out_dir} \
		-n ${num_metacell_limit} \
		-f
done
