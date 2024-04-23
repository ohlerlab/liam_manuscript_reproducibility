# Download PBMC Multiome data whole cells from TEA-seq publication (Swanson et al., 2021)
Author: Pia Rautenstrauch

Date: 2023/17/10

## Reference
Swanson, E. et al. (2021) Simultaneous trimodal single-cell measurement of transcripts, epitopes, and chromatin accessibility using TEA-seq. eLife. DOI: 10.7554/eLife.63632 
    
## Data source
GEO- Series: GSE158013

Sample: GSM5123950 - Multiome RNA/ATAC of leukapheresis-purified, 0.01% digitonin perm, FACS neutrophil-depleted PBMCs, method comparison

## General notes
Fragment files are decompressed and recompressed, as ScregSeg (https://github.com/BIMSBbioinfo/scregseg/), which I use to generate ATAC features, sometimes fails on compressed CellRanger fragment files.

### Multiome perm cells
Working directory: ```liam_manuscript_reproducibility/data/original/Swanson2021/perm_cells```

#### ATAC
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5123nnn/GSM5123950/suppl/GSM5123950%5FX066%2DMP0C1W2%5Fleukopak%5Fperm%2Dcells%5Fmultiome%5F200M%5Fatac%5Ffiltered%5Ffragments%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5123nnn/GSM5123950/suppl/GSM5123950%5FX066%2DMP0C1W2%5Fleukopak%5Fperm%2Dcells%5Fmultiome%5F200M%5Fatac%5Ffiltered%5Fmetadata%2Ecsv%2Egz
gunzip -f *fragments.tsv.gz
gzip -f *fragments.tsv
```

#### RNA
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5123nnn/GSM5123950/suppl/GSM5123950%5FX066%2DMP0C1W2%5Fleukopak%5Fperm%2Dcells%5Fmultiome%5F200M%5Fcellranger%2Darc%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5123nnn/GSM5123950/suppl/GSM5123950%5FX066%2DMP0C1W2%5Fleukopak%5Fperm%2Dcells%5Fmultiome%5F200M%5Fcellranger%2Darc%5Fper%5Fbarcode%5Fmetrics%2Ecsv%2Egz
```