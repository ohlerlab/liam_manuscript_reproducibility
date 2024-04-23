# Software versions
# samtools-1.12 and bedtools2 (v2.29.1) need to be on PATH!
# wd: liam_manuscript_reproducibility/data/derived/Mimitou2021/ATAC/
conda activate lpr_scregseg

for i in "./../../data/original/Mimitou2021/ASAP_seq/LLL_CTRL/ATAC/GSM4732109_CD28_CD3_control_ASAP_fragments.tsv.gz ASAP_LLL_CTRL" "./../../data/original/Mimitou2021/ASAP_seq/LLL_STIM/ATAC/GSM4732111_CD28_CD3_stim_ASAP_fragments.tsv.gz ASAP_LLL_STIM"  "./../../data/original/Mimitou2021/DOGMA_seq/LLL_CTRL/ATAC/GSM5065524_LLL_CTRL_fragments.tsv.gz DOGMA_LLL_CTRL" "./../../data/original/Mimitou2021/DOGMA_seq/LLL_STIM/ATAC/GSM5065527_LLL_STIM_fragments.tsv.gz DOGMA_LLL_STIM" "./../../data/original/Mimitou2021/DOGMA_seq/DIG_CTRL/ATAC/GSM5065530_DIG_CTRL_fragments.tsv.gz DOGMA_DIG_CTRL" "./../../data/original/Mimitou2021/DOGMA_seq/DIG_STIM/ATAC/GSM5065533_DIG_STIM_fragments.tsv.gz DOGMA_DIG_STIM" ;
do
path_name_pair=( $i ); 
echo "${path_name_pair[1]}";
echo "${path_name_pair[0]}";

echo "scregseg make_tile \
          --regions tile1kb_${path_name_pair[1]}.bed \
          --binsize 1000 \
          --fragmentfile ${path_name_pair[0]} \
          --keep_chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"

scregseg make_tile \
          --regions tile1kb_${path_name_pair[1]}.bed \
          --binsize 1000 \
          --fragmentfile ${path_name_pair[0]} \
          --keep_chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22          

echo "scregseg fragments_to_counts \
          --fragmentfile ${path_name_pair[0]} \
          --regions tile1kb_${path_name_pair[1]}.bed \
          --with-fraglen \
          --counts countmatrix_${path_name_pair[1]}.h5ad"

scregseg fragments_to_counts \
          --fragmentfile ${path_name_pair[0]} \
          --regions tile1kb_${path_name_pair[1]}.bed \
          --with-fraglen \
          --counts countmatrix_${path_name_pair[1]}.h5ad
          
echo "scregseg subset \
          --incounts countmatrix_${path_name_pair[1]}.h5ad \
          --outcounts subsetted_countmatrix_${path_name_pair[1]}.h5ad \
          --subset qcontrolled_cells_${path_name_pair[1]}.csv"
          
scregseg subset \
          --incounts countmatrix_${path_name_pair[1]}.h5ad \
          --outcounts subsetted_countmatrix_${path_name_pair[1]}.h5ad \
          --subset qcontrolled_cells_${path_name_pair[1]}.csv

echo "scregseg filter \
          --incounts subsetted_countmatrix_${path_name_pair[1]}.h5ad \
          --outcounts filtered_countmatrix_${path_name_pair[1]}.h5ad \
          --minregioncount 1 \
          --trimcount 1"
          
scregseg filter \
          --incounts subsetted_countmatrix_${path_name_pair[1]}.h5ad \
          --outcounts filtered_countmatrix_${path_name_pair[1]}.h5ad \
          --minregioncount 1 \
          --trimcount 1
          
done