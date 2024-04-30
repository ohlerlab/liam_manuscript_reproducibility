Author: Pia Rautenstrauch

## Use case to directory mapping
The data folders follow a different naming convention than presented in the manuscript. For orientation, consider the following mapping:
- TCU: Mimitou2021
- ETCU: Mimitou2021 + 10x + Swanson2021 + MOFA+
- CU: neurips_competition
- CU_mosaic: neurips_competition

## Original data
References, links, etc., are provided in READMEs. 

## Derived data
You can derive all processed data from the raw data following the preprocessing scripts provided in ```./../scripts/<use_case>/```, respectively. 

Since it involves many steps to derive the preprocessed data for the TCU and ETCU use cases, we provide them here.

- ```liam_manuscript_reproducibility/data/derived/Mimitou2021/DOGMA_seq/extended_treatment_control_use_case_revisions_r1.h5ad.gz```
- ```liam_manuscript_reproducibility/data/derived/Mimitou2021/DOGMA_seq/preprocessed_DOGMA.h5ad.gz```

You need to ```gunzip``` them before use.

The data sets for the CU and CU_mosaic use cases can be easily recreated from the raw data following the scripts in ```./../scripts/CU/``` and ```./../scripts/CU_mosaic/```, respectively.
