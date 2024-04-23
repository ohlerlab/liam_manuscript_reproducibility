# Download PBMC from a Healthy Donor - Granulocytes Removed Through Cell Sorting (10k) from 10x
Author: Pia Rautenstrauch

Date: 2023/25/10

## Data source
https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-1-0-0

## General notes
Fragment files are decompressed and recompressed, as ScregSeg (https://github.com/BIMSBbioinfo/scregseg/), which I use to generate ATAC features, sometimes fails on compressed CellRanger fragment files.

### PBMC from a Healthy Donor - No Cell Sorting (10k)
Workind directory: ```liam_manuscript_reproducibility/data/original/10x/10k_sorted```

#### ATAC
```
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
gunzip -f *fragments.tsv.gz
gzip -f *fragments.tsv
```

#### RNA
```
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
```