# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2022/22/07

### Train MultiVI with default params on phase2 data
#### MultiVI with ATAC 10% subsample (5 different subsamples)
#### MultiVI with ATAC 25% subsample (5 different subsamples)
# Liam_challenge_reproducibility environment
## Following the tutorial, using extra categorical covariates for batch despite having only one modality
### https://docs.scvi-tools.org/en/stable/tutorials/notebooks/MultiVI_tutorial.html (as of 2022-15-04)

# Imports
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
import anndata as ad
import os
import scvi
import scipy

print("CUDA available?", torch.cuda.is_available())

print("scvi version:", scvi.__version__)

# Setup model to param mapping

model_param_mapping = {}

seeds = [8831, 234, 11, 9631, 94]

for seed in seeds:
    # Note that this seed is not the random seed for the model training, but the random seed for the ATAC data subsampling!
    model_param_mapping['MultiVI_10_seed_{}'.format(seed)] = {}
    model_param_mapping['MultiVI_10_seed_{}'.format(seed)]['path_to_data'] = './../../data/derived/neurips_competition/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2_10_subsample_seed_{}.h5ad'.format(seed)
    model_param_mapping['MultiVI_25_seed_{}'.format(seed)] = {}
    model_param_mapping['MultiVI_25_seed_{}'.format(seed)]['path_to_data'] = './../../data/derived/neurips_competition/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2_25_subsample_seed_{}.h5ad'.format(seed)

# Train models
for model in model_param_mapping.keys():
    print(model)    
    # All model have the same seed (0)
    scvi._settings.ScviConfig()

    # Load data
    ad_mod1 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad")
    ad_mod2 = ad.read_h5ad(model_param_mapping[model]['path_to_data'])
    solution = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_solution.h5ad")

    # Setup input object
    input = ad.AnnData(scipy.sparse.hstack([ad_mod1.layers["counts"], ad_mod2.X]).tocsr(), ad_mod1.obs, pd.concat((ad_mod1.var, ad_mod2.var), axis=0))
    input.obs["sample"] = input.obs["batch"]
    input.obs["donor"] = input.obs["batch"].apply(lambda x: x.split("d")[1])
    input.obs["site"] = input.obs["batch"].apply(lambda x: x.split("d")[0])
    input.obs["cell type"] = solution.obs["cell_type"][input.obs.index]
    del ad_mod1
    del ad_mod2
    del solution

    input = scvi.data.organize_multiome_anndatas(input)

    print("Filtering to retain only features present in more than 1% of the cells.")
    print(input.shape)
    sc.pp.filter_genes(input, min_cells=int(input.shape[0] * 0.01))
    print(input.shape)

    # Setting up data.
    scvi.model.MULTIVI.setup_anndata(input, batch_key='modality', categorical_covariate_keys=['sample'])
    # Note difference here: batch_key='modality', categorical_covariate_keys=['sample']
    # important to specify them as list! 
    vae = scvi.model.MULTIVI(input,
        n_genes=(input.var['feature_types']=='GEX').sum(),
        n_regions=(input.var['feature_types']=='ATAC').sum(),
                            )

    # Training models
    vae.train()


    input.obsm["embedding"] = vae.get_latent_representation()

    # Save trained model
    vae.save("Models/{}".format(model), save_anndata=True)

    input.uns['dataset_id'] = model

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
    
    # For reproducibility
    print('Model: {}.'.format(model))
    print("Model's state_dict:")
    for param_tensor in vae.module.state_dict():
        print(param_tensor, '\t', vae.module.state_dict()[param_tensor].size())
    
    del vae
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