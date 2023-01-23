#!/bin/bash

# 1. run_scf_repeatedly_ with different knn and downsampling
# output: knn_xy matrices, and cells (i, knn)

# 2. within-modal clustering ...
# for each cell set from (i, knn) -> run Leiden clustering with different resolution r
# output: clustering result.tsv (i, knn, r): (i, knn) separate tables, r columns per table

# 3. correlation analysis ...
# for each combination of (i, knn, r), run a correlation analysis
# output: corrs.pkl (i, knn, r) 

date="211115" # still use the (finished) 1130 results
out_dir="./results"
data_dir="./newdata"

modx='rna'
mody='atac'
input_datasets="profiles_hvgene_atac.h5ad profiles_hvgene_rna.h5ad"
input_modalities="atac rna"
feature_datasets="profiles_hvgene_atac.h5ad"
input_dataset_rna="./newdata/profiles_hvgene_rna.h5ad"
modx_step3="profiles_hvgene_rna"
mody_step3="profiles_hvgene_atac"

ka=30
knn=30
corr_type='spearmanr'
subsample_frac=0.1 # take 10% of all the data -- faster for demo
subsample_times=1
resolutions=(10) # Leiden clustering resolutions used to generate metacells -- just 1 for demo
num_metacell_limit=1001

# generate (i, knn) knn_xy matrices
# modalities
# scf config template
nameTag="mop_${modx}_${mody}_ka${ka}_knn${knn}_${date}" # need to make sure they are consistent with the config_template
inputNameTag="mop_${modx}_${mody}_ka${ka}_knn{}_${date}"
echo $modx, $mody, $ka, $knn

# # 1.
# # # run SCF
# echo "STEP1..."
# python robustlink scf \
# 	-i ${data_dir} \
# 	-id ${input_datasets} \
# 	-im ${input_modalities} \
# 	-fd ${feature_datasets} \
# 	-o ${out_dir} \
# 	-on ${nameTag} \
# 	-s ${subsample_frac} \
# 	-sn ${subsample_times}

# # # 2.
# # run leiden clustering for each (i, knn) 
# echo "STEP2..."
# python robustlink metacell \
# 	-i  ${input_dataset_rna} \
# 	-o  ${out_dir} \
# 	--knns $knn \
# 	-sn ${subsample_times} \
# 	-tag ${inputNameTag} \
# 	-r ${resolutions}

# 3.
# # correlation analysis (i, knn, r)
echo "STEP3..."
for (( i=0; i<${subsample_times}; i++ )); do
	python robustlink corr_atac \
		-modx ${modx_step3} \
		-mody ${mody_step3} \
		-tag $nameTag \
		-isub $i \
		--corr_type ${corr_type} \
		-n ${num_metacell_limit} \
		-f
done