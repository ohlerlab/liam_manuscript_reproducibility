# Software versions
# samtools-1.12 and bedtools2 (v2.29.1) need to be on PATH!
# wd: liam_manuscript_reproducibility/data/derived/Mimitou2021/ATAC/

conda activate lpr_scregseg

scregseg fit_segment --counts all_datasets_input_scregseg_second_pass.h5ad \
    --storage scregseg_fi_all_samples_50_states \
    --randomseed 32 7 943 301 72 580 714 --n_jobs 64 --nstates 50 --niter 300 

