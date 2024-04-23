# Software versions
# samtools-1.12 and bedtools2 (v2.29.1) need to be on PATH!
# wd: liam_manuscript_reproducibility/data/derived/Mimitou2021/ATAC/
conda activate lpr_scregseg

scregseg seg_to_bed --storage scregseg_fi_all_samples_50_states --output all_datasets_informative_rarest_states_out_of_50.bed --method abundancethreshold --threshold 0.9 --no_bookended_merging --counts all_datasets_input_scregseg_second_pass.h5ad --max_state_abundance 0.015 


