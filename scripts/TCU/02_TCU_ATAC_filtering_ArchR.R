# Setup environment ----
# guix env: lpr_R (see session info at the bottom of this file)

# Set working directory
setwd("liam_manuscript_reproducibility/scripts/TCU")

# Filter fragment files from PBMC stimulation experiments (ASAP-seq and DOGMA-seq data) from Mimitou et al., 2021 ----
# Retain only barcodes with an TSS enrichment > 4 and > 1000 fragments

# Imports
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
# Needs to be added manually (only fixed in ArchR development version https://github.com/GreenleafLab/ArchR/issues/1197)
library(parallel)

# Paths relative to working directory
paths <- c("./../../data/derived/Mimitou2021/ASAP_seq/LLL_CTRL/GSM4732109_CD28_CD3_control_ASAP_fragments.tsv.gz",
           "./../../data/derived/Mimitou2021/ASAP_seq/LLL_STIM/GSM4732111_CD28_CD3_stim_ASAP_fragments.tsv.gz",
           "./../../data/derived/Mimitou2021/DOGMA_seq/LLL_CTRL/GSM5065524_LLL_CTRL_fragments.tsv.gz",
           "./../../data/derived/Mimitou2021/DOGMA_seq/LLL_STIM/GSM5065527_LLL_STIM_fragments.tsv.gz",           
           "./../../data/derived/Mimitou2021/DOGMA_seq/DIG_CTRL/GSM5065530_DIG_CTRL_fragments.tsv.gz",
           "./../../data/derived/Mimitou2021/DOGMA_seq/DIG_STIM/GSM5065533_DIG_STIM_fragments.tsv.gz"
           )
           

addArchRGenome("hg38")
addArchRThreads(threads = 30) 

sampleNames <- c("ASAP_LLL_CTRL", "ASAP_LLL_STIM", "DOGMA_LLL_CTRL", "DOGMA_LLL_STIM", "DOGMA_DIG_CTRL", "DOGMA_DIG_STIM")
ArrowFiles <- createArrowFiles(
  inputFiles = paths,
  sampleNames = sampleNames,
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = FALSE,
  addGeneScoreMat = FALSE
)

# Establish ArchR projects and save arrow files
path_arrow <- "./../../data/derived/Mimitou2021/ATAC/"

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = path_arrow,
  copyArrows = TRUE
)

# Remove files from working directory
file.remove("DOGMA_DIG_STIM.arrow")
file.remove("DOGMA_DIG_CTRL.arrow")
file.remove("DOGMA_LLL_STIM.arrow")
file.remove("DOGMA_LLL_CTRL.arrow")
file.remove("ASAP_LLL_CTRL.arrow")
file.remove("ASAP_LLL_STIM.arrow")

# Save cell ids/sample names of cells that passed filtering to csv for further preprocessing with ScregSeg 
# Indices of cells that remain after initial filtering based on TSS enrichment and fragment counts
all_cells <- proj$cellNames
  
write.csv(all_cells, paste0("./../../data/derived/Mimitou2021/ATAC/qcontrolled_cells_all_ATAC_samples", ".csv"), row.names = FALSE)

# Indices of cells of individual samples remaining after initial filtering based on TSS enrichment and fragment counts
samples <- sapply(strsplit(all_cells, "#"), `[`, 1)
barcodes <- sapply(strsplit(all_cells, "#"), `[`, 2)
all_cells_df <- data.frame(samples, barcodes)

for ( sample in sampleNames) { 
  tmp_df <- subset(all_cells_df, all_cells_df$samples == sample)
  write.csv(tmp_df$barcodes, paste0("./../../data/derived/Mimitou2021/ATAC/qcontrolled_cells_", sample, ".csv"), row.names = FALSE)
  nb_cells_remaining <- dim(tmp_df)[1]
  print(paste0("For sample: ", sample, " ", nb_cells_remaining, " cells remain after initial filtering with ArchR."))
}

# [1] "For sample: ASAP_LLL_CTRL 4477 cells remain after initial filtering with ArchR."
# [1] "For sample: ASAP_LLL_STIM 5836 cells remain after initial filtering with ArchR."
# [1] "For sample: DOGMA_LLL_CTRL 8031 cells remain after initial filtering with ArchR."
# [1] "For sample: DOGMA_LLL_STIM 6718 cells remain after initial filtering with ArchR."
# [1] "For sample: DOGMA_DIG_CTRL 11245 cells remain after initial filtering with ArchR."
# [1] "For sample: DOGMA_DIG_STIM 11054 cells remain after initial filtering with ArchR."


sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-unknown-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /gnu/store/bs9pl1f805ins80xaf4s3n35a0x2lyq3-openblas-0.3.9/lib/libopenblasp-r0.3.9.so
# 
# locale:
# [1] C
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] BSgenome.Hsapiens.UCSC.hg38_1.4.1 BSgenome_1.62.0                   rtracklayer_1.54.0                Biostrings_2.62.0                
# [5] XVector_0.34.0                    ArchR_1.0.0                       magrittr_2.0.1                    rhdf5_2.38.0                     
# [9] Matrix_1.3-4                      data.table_1.14.2                 SummarizedExperiment_1.24.0       Biobase_2.54.0                   
# [13] GenomicRanges_1.46.1              GenomeInfoDb_1.30.0               IRanges_2.28.0                    S4Vectors_0.32.3                 
# [17] BiocGenerics_0.40.0               MatrixGenerics_1.6.0              matrixStats_0.61.0                ggplot2_3.3.5                    
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.1.1         purrr_0.3.4              lattice_0.20-45          colorspace_2.0-2         vctrs_0.3.8              generics_0.1.1          
# [7] yaml_2.2.1               utf8_1.2.2               XML_3.99-0.8             rlang_0.4.12             pillar_1.6.4             glue_1.5.0              
# [13] withr_2.4.2              BiocParallel_1.28.1      GenomeInfoDbData_1.2.0   lifecycle_1.0.1          stringr_1.4.0            munsell_0.5.0           
# [19] gtable_0.3.0             restfulr_0.0.13          fansi_0.5.0              Rcpp_1.0.7               scales_1.1.1             DelayedArray_0.20.0     
# [25] Rsamtools_2.10.0         rjson_0.2.20             stringi_1.7.5            dplyr_1.0.7              BiocIO_1.4.0             grid_4.1.2              
# [31] tools_4.1.2              rhdf5filters_1.6.0       bitops_1.0-7             RCurl_1.95-0.1.2         tibble_3.1.6             crayon_1.4.2            
# [37] pkgconfig_2.0.3          ellipsis_0.3.2           Rhdf5lib_1.16.0          R6_2.5.1                 GenomicAlignments_1.30.0 compiler_4.1.2          