# liam_manuscript_reproducibility
Author: Pia Rautenstrauch

This Zenodo record contains scripts, notebooks, data, and environment configuration files essential for reproducing the data preprocessing, analyses, and figures presented in the "Liam tackles complex multimodal single-cell data integration challenges" peer-reviewed publication. Data preprocessing and analyses for models presented in the preprint [[1]](#1) are equivalent. It also hosts intermediate data created that requires many preprocessing steps to recreate. If you merely search for an example, consider the tutorials at https://github.com/ohlerlab/liam, which hosts the software under active development.

## Notes on analysis scripts and notebooks
I reorganized the scripts for readability. As a result, some of the relative paths/absolute paths might be erroneous but should be recoverable from the directory structure and naming conventions. I omitted some cross-referenced files, notebooks, and analysis scripts that are not required to reproduce the paper's results.

### Directory structure
In every ```scripts/``` subfolder (TCU - treatment-control use case, ETCU - extended treatment-control use case, CU - competition use case, and CU_mosaic - mosaic use case) are:
- scripts/notebooks detailing the preprocessing steps to prepare the data
- scripts/notebooks showcasing model parameters and training of liam [1] and models we compare to
- evaluation scripts/notebooks 
- notebooks for generating the figures from the publication

For mapping the model names in the scripts to model names in the publication, consider the notebooks for generating publication figures, which detail the mapping. 

All scripts/notebooks mention which conda environment/environment they use in comments, with yml or config files detailed in ```environments/```.

## Notes on intermediate data
The directory ```data/original/``` describes how to obtain raw data/published data used in the manuscript.

We detail how we further processed this data in scripts and notebooks in the respective ```scripts/``` subdirectories.

Technically, you can recreate the input data for the TCU and ETCU use cases from the presented scripts and notebooks. However, this involves many steps, and we thus include the intermediate anndata objects.
- TCU: data/derived/Mimitou2021/DOGMA_seq/preprocessed_DOGMA.h5ad.gz
- ETCU: data/derived/Mimitou2021/DOGMA_seq/extended_treatment_control_use_case_revisions_r1.h5ad.gz

We provide intermediate outputs (embeddings) sufficient for generating publication figures for all use cases.

For CU and CU_mosaic: the embeddings are in respective folders ```scripts/<use_case>Predictions/```.

For TCU and ETCU: the embeddings can be loaded from ```models/TCU``` and ```models/ETCU``` in models.zip.

## Notes on different liam versions
For software from the publication under active development, see https://github.com/ohlerlab/liam.

For the legacy version of the software used in the presented scripts, see https://github.com/ohlerlab/liam_challenge_reproducibility. We performed the presented analyses with different versions of the legacy software. For details, see the GitHub repos release notes and the ```environments/``` README. All added functionality is theoretically downward compatible (except for a bug of the mosaic functionality, breaking single-modality models). Version v0.1.1 of [liam](https://github.com/ohlerlab/liam) under active development comprises all functionalities and a downward compatible bug fix of the mosaic functionality. Running liam v0.1.1 will thus result in equivalent results, as obtained with the distinct liam_challenge_reproducibility releases.

## References
<a id="1">[1]</a>
Rautenstrauch, P. & Ohler, U. (2022) [Liam tackles complex multimodal single-cell data integration challenges](https://www.biorxiv.org/content/10.1101/2022.12.21.521399v1). bioRxiv DOI: 10.1101/2022.12.21.521399
