# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2022/18/05

### Train models on phase2 full data (Liam)
#### For five different random seeds
#### GEX only - sample 10 dims
#### ATAC only - sample 10 dims
#### combine into one concat model! (20 dims)


# Liam_challenge_reproducibility environment

# Imports
import liam_NeurIPS2021_challenge_reproducibility
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
import anndata as ad
import os
import scvi

print(torch.cuda.is_available())

# Setup model to param mapping
model_path_mapping = {}

seeds = [0, 994, 236, 71, 415]
for seed in seeds:
    model_path_mapping['BAVAE_sample_100_concat_seed_{}'.format(seed)] = ('BAVAE_sample_100_rna_only_seed_{}'.format(seed), 'BAVAE_sample_100_atac_only_seed_{}'.format(seed))

# Train models
for model in model_path_mapping.keys():
    print(model)
    scvi._settings.ScviConfig()
    
    # Load trained individual models 
    vae_rna = liam_NeurIPS2021_challenge_reproducibility.Liam.load("Models/{}".format(model_path_mapping[model][0]))
    vae_atac = liam_NeurIPS2021_challenge_reproducibility.Liam.load("Models/{}".format(model_path_mapping[model][1]))

    # Setup input object
    input = vae_atac.adata.copy()
    input.obsm["embedding_atac"] = vae_atac.get_latent_representation()
    input.obsm["embedding_gex"] = vae_rna.get_latent_representation()
    
    # Concatenate embeddings
    input.obsm["embedding"] = np.concatenate((input.obsm["embedding_atac"],
                                              input.obsm["embedding_gex"]), axis=1)

    adata = ad.AnnData(
        X=input.X,
        obs=input.obs,
        uns={
            'dataset_id': input.uns['dataset_id'],
            'method_id': model,
        },
        obsm={'embedding' : input.obsm['embedding']}
    )

    evaluation = ad.AnnData(
        X=input.obsm['embedding'],
        obs=input.obs,
        uns={
            'dataset_id': input.uns['dataset_id'],
            'method_id': model,
        }
    )

    adata.write_h5ad("Embeddings/adata_{}.h5ad".format(model), compression="gzip")
    evaluation.write_h5ad("Predictions/{}.prediction.h5ad".format(model), compression="gzip")

    del vae_atac
    del vae_rna
    del input
    
# For reproducibility across all trained models
print('! nvidia-smi')
print(os.system('nvidia-smi')) 
print('torch.cuda.get_device_name()')
print(torch.cuda.get_device_name())
print('torch.version.cuda')
print(torch.version.cuda)
print('torch.cuda.get_device_capability()')
print(torch.cuda.get_device_capability())
print('torch.cuda.get_device_properties(torch.device)')
print(torch.cuda.get_device_properties(torch.device))
print('torch.cuda.get_arch_list()')
print(torch.cuda.get_arch_list())