# Directory content
Author: Pia Rautenstrauch

In contrast to the manuscript, the variable 'replicate' is called 'buffer' in these scripts and outputs. We initially used this terminology, as the original data authors generated the replicates using two distinct buffers (lysis conditions).

There might be discrepancies between the path in Jupyter notebook cells and their outputs as I reorganized the repository structure for readability. Please consider the paths from the notebook cells and the general organization of this repository.

Data preprocessing: 
- Overview: 
    - notebook 00

- ATAC preprocessing 
    - derivation of joint feature set for ATAC: notebooks/scripts 01 - 07

- GEX and ADT preprocessing and combination with ATAC into final model input data set:
    - notebook 08
    
Model training: notebook 09

Export of embeddings: notebook 10
    
Scripts/notebooks for evaluation (iLISI) and publication figures are in the subfolder ```publication_figures/```.

Evaluation scores: ```Evaluation/```.
