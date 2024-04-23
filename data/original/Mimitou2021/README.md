# Download DOGMA-seq and ASAP-seq data from T cell stimulation experiments (Mimitou et al., 2021)
Author: Pia Rautenstrauch

Date: 2021/11/30

## Reference
Mimitou, E.P. et al. (2021) Scalable, multimodal profiling of chromatin accessibility, gene expression and protein levels in single cells. Nat. Biotechnol. DOI: 10.1038/s41587-021-00927-2

## Data source 
GEO - Superseries: GSE156478

### DOGMA-seq T cell stimulation experiments
SubSeries: GSE166188: High throughput joint profiling of chromatin accessibility and protein levels in single cells [PBMC_Stim_Multiome]

### ASAP-seq T cell stimulation experiments
SubSeries: GSE156473: High throughput joint profiling of chromatin accessibility and protein levels in single cells [PBMC stimulation]


## General notes
Fragment files are decompressed and recompressed, as ScregSeg (https://github.com/BIMSBbioinfo/scregseg/), which I use to generate ATAC features, sometimes fails on compressed CellRanger fragment files.


## DOGMA-seq
### LLL_CTRL
#### ATAC
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/LLL_CTRL/ATAC```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065524/suppl/GSM5065524%5FLLL%5FCTRL%5Ffragments%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065524/suppl/GSM5065524%5FLLL%5FCTRL%5Fper%5Fbarcode%5Fmetrics%2Ecsv%2Egz
gunzip -f *_fragments.tsv.gz
gzip -f *_fragments.tsv
```

#### RNA
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/LLL_CTRL/RNA```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065525/suppl/GSM5065525%5FLLL%5FCTRL%5FGExp%5FATAC%5Ffiltered%5Fbarcodes%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065525/suppl/GSM5065525%5FLLL%5FCTRL%5FGExp%5FATAC%5Ffiltered%5Ffeatures%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065525/suppl/GSM5065525%5FLLL%5FCTRL%5FGExp%5FATAC%5Ffiltered%5Fmatrix%2Emtx%2Egz
````

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_GExp_ATAC_filtered_matrix.mtx.gz matrix.mtx.gz
ln -s *_GExp_ATAC_filtered_features.tsv.gz features.tsv.gz
ln -s *_GExp_ATAC_filtered_barcodes.tsv.gz barcodes.tsv.gz
```

#### ADT
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/LLL_CTRL/ADT```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065526/suppl/GSM5065526%5FLLL%5Fctrl%5FADT%5FallCounts%2Ebarcodes%2Etxt%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065526/suppl/GSM5065526%5FLLL%5Fctrl%5FADT%5FallCounts%2Emtx%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065526/suppl/GSM5065526%5FLLL%5Fctrl%5FADT%5FallCounts%2Eproteins%2Etxt%2Egz
```

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_ADT_allCounts.mtx.gz matrix.mtx.gz
ln -s *_ADT_allCounts.proteins.txt.gz features.tsv.gz
ln -s *_ADT_allCounts.barcodes.txt.gz barcodes.tsv.gz
```


### LLL_STIM
#### ATAC
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/LLL_STIM/ATAC```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065527/suppl/GSM5065527%5FLLL%5FSTIM%5Ffragments%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065527/suppl/GSM5065527%5FLLL%5FSTIM%5Fper%5Fbarcode%5Fmetrics%2Ecsv%2Egz
gunzip -f *_fragments.tsv.gz
gzip -f *_fragments.tsv
```

#### RNA
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/LLL_STIM/RNA```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065528/suppl/GSM5065528%5FLLL%5FSTIM%5FGExp%5FATAC%5Ffiltered%5Fbarcodes%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065528/suppl/GSM5065528%5FLLL%5FSTIM%5FGExp%5FATAC%5Ffiltered%5Ffeatures%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065528/suppl/GSM5065528%5FLLL%5FSTIM%5FGExp%5FATAC%5Ffiltered%5Fmatrix%2Emtx%2Egz
````

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_GExp_ATAC_filtered_matrix.mtx.gz matrix.mtx.gz
ln -s *_GExp_ATAC_filtered_features.tsv.gz features.tsv.gz
ln -s *_GExp_ATAC_filtered_barcodes.tsv.gz barcodes.tsv.gz
```

#### ADT
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/LLL_STIM/ADT```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065529/suppl/GSM5065529%5FLLL%5Fstim%5FADT%5FallCounts%2Ebarcodes%2Etxt%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065529/suppl/GSM5065529%5FLLL%5Fstim%5FADT%5FallCounts%2Emtx%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065529/suppl/GSM5065529%5FLLL%5Fstim%5FADT%5FallCounts%2Eproteins%2Etxt%2Egz
```

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_ADT_allCounts.mtx.gz matrix.mtx.gz
ln -s *_ADT_allCounts.proteins.txt.gz features.tsv.gz
ln -s *_ADT_allCounts.barcodes.txt.gz barcodes.tsv.gz
```


### DIG_CTRL
#### ATAC
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/DIG_CTRL/ATAC```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065530/suppl/GSM5065530%5FDIG%5FCTRL%5Ffragments%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065530/suppl/GSM5065530%5FDIG%5FCTRL%5Fper%5Fbarcode%5Fmetrics%2Ecsv%2Egz
gunzip -f *_fragments.tsv.gz
gzip -f *_fragments.tsv
```

#### RNA
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/DIG_CTRL/RNA```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065531/suppl/GSM5065531%5FDIG%5FCTRL%5FGExp%5FATAC%5Ffiltered%5Fbarcodes%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065531/suppl/GSM5065531%5FDIG%5FCTRL%5FGExp%5FATAC%5Ffiltered%5Ffeatures%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065531/suppl/GSM5065531%5FDIG%5FCTRL%5FGExp%5FATAC%5Ffiltered%5Fmatrix%2Emtx%2Egz
````

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_GExp_ATAC_filtered_matrix.mtx.gz matrix.mtx.gz
ln -s *_GExp_ATAC_filtered_features.tsv.gz features.tsv.gz
ln -s *_GExp_ATAC_filtered_barcodes.tsv.gz barcodes.tsv.gz
```

#### ADT
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/DIG_CTRL/ADT```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065532/suppl/GSM5065532%5FDIG%5Fctrl%5FADT%5FallCounts%2Ebarcodes%2Etxt%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065532/suppl/GSM5065532%5FDIG%5Fctrl%5FADT%5FallCounts%2Emtx%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065532/suppl/GSM5065532%5FDIG%5Fctrl%5FADT%5FallCounts%2Eproteins%2Etxt%2Egz
```

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_ADT_allCounts.mtx.gz matrix.mtx.gz
ln -s *_ADT_allCounts.proteins.txt.gz features.tsv.gz
ln -s *_ADT_allCounts.barcodes.txt.gz barcodes.tsv.gz
```


### DIG_STIM
#### ATAC
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/DIG_STIM/ATAC```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065533/suppl/GSM5065533%5FDIG%5FSTIM%5Ffragments%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065533/suppl/GSM5065533%5FDIG%5FSTIM%5Fper%5Fbarcode%5Fmetrics%2Ecsv%2Egz
gunzip -f *_fragments.tsv.gz
gzip -f *_fragments.tsv
```

#### RNA
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/DIG_STIM/RNA```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065534/suppl/GSM5065534%5FDIG%5FSTIM%5FGExp%5FATAC%5Ffiltered%5Fbarcodes%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065534/suppl/GSM5065534%5FDIG%5FSTIM%5FGExp%5FATAC%5Ffiltered%5Ffeatures%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065534/suppl/GSM5065534%5FDIG%5FSTIM%5FGExp%5FATAC%5Ffiltered%5Fmatrix%2Emtx%2Egz
````

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_GExp_ATAC_filtered_matrix.mtx.gz matrix.mtx.gz
ln -s *_GExp_ATAC_filtered_features.tsv.gz features.tsv.gz
ln -s *_GExp_ATAC_filtered_barcodes.tsv.gz barcodes.tsv.gz
```

#### ADT
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/DIG_STIM/ADT```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065535/suppl/GSM5065535%5FDIG%5Fstim%5FADT%5FallCounts%2Ebarcodes%2Etxt%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065535/suppl/GSM5065535%5FDIG%5Fstim%5FADT%5FallCounts%2Emtx%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/GSM5065535/suppl/GSM5065535%5FDIG%5Fstim%5FADT%5FallCounts%2Eproteins%2Etxt%2Egz
```

Create softlinks enabling 10X-like formatted input reading with scanpy
```
ln -s *_ADT_allCounts.mtx.gz matrix.mtx.gz
ln -s *_ADT_allCounts.proteins.txt.gz features.tsv.gz
ln -s *_ADT_allCounts.barcodes.txt.gz barcodes.tsv.gz
```

## ASAP-seq
### LLL_CTRL
#### ATAC
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/ASAP_seq/LLL_CTRL/ATAC```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4732nnn/GSM4732109/suppl/GSM4732109%5FCD28%5FCD3%5Fcontrol%5FASAP%5Ffragments%2Etsv%2Egz
gunzip -f *_fragments.tsv.gz
gzip -f *_fragments.tsv
```

#### ADT
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/ASAP_seq/LLL_CTRL/ADT```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4732nnn/GSM4732110/suppl/GSM4732110%5FCD28%5FCD3%5Fcontrol%5FASAP%5FADT%2Etsv%2Egz
```


### LLL_STIM
#### ATAC
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/ASAP_seq/LLL_STIM/ATAC```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4732nnn/GSM4732111/suppl/GSM4732111%5FCD28%5FCD3%5Fstim%5FASAP%5Ffragments%2Etsv%2Egz
gunzip -f *_fragments.tsv.gz
gzip -f *_fragments.tsv
```

#### ADT
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/ASAP_seq/LLL_STIM/ADT```
```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4732nnn/GSM4732112/suppl/GSM4732112%5FCD28%5FCD3%5Fstim%5FASAP%5FADT%2Etsv%2Egz
```

### Preprocessed R object to select same cells as in Mimitou et al., 2021 
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/ASAP_seq/from_asap_large_data_files```
Downloaded: ArchR_PBMCs_stim.rds Version ID: 015ade52 from https://osf.io/96va3/, navigate to: ```asap_large_data_files/pbmc_stim_data/output/archr_pbmc_stim```



### Preproccesed R object to get celltype annotations for DOGMA_seq LLL
Working directory: ```liam_manuscript_reproducibility/data/original/Mimitou2021/DOGMA_seq/multiome_pbmc_stim/output/
Downloaded: pbmc_LLL_processed.rds Version ID: 015b9cb8 from https://osf.io/6kr4v/download  