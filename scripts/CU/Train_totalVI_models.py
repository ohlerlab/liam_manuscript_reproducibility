# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2022/17/05

### Train totalVI model  on phase2 full data 
# Liam_challenge_reproducibility environment
## Following the totalVI tutorial
### https://docs.scvi-tools.org/en/stable/tutorials/notebooks/totalVI.html (as of 2022-15-04)

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

print("CUDA available?", torch.cuda.is_available())

print("scvi version:", scvi.__version__)


seeds = [0, 994, 236, 71, 415]

for seed in seeds:
    # Setup model to param mapping
    model = 'totalVI_seed_{}'.format(seed)
    print(model)
    scvi._settings.ScviConfig(seed=seed)

    # Load data
    ad_mod1 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_cite_phase2/openproblems_bmmc_cite_phase2.censor_dataset.output_mod1.h5ad")
    ad_mod2 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_cite_phase2/openproblems_bmmc_cite_phase2.censor_dataset.output_mod2.h5ad")
    solution = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_cite_phase2/openproblems_bmmc_cite_phase2.censor_dataset.output_solution.h5ad")
    
    # Setup input object
    input = ad_mod1.copy()
    input.layers["counts"] = input.layers["counts"].copy().tocsr()
    input.X = input.layers["counts"].copy().tocsr()
    input.obsm["ADT"] = ad_mod2.X.todense()
    input.obs["sample"] = input.obs["batch"]
    input.obs["donor"] = input.obs["batch"].apply(lambda x: x.split("d")[1])
    input.obs["site"] = input.obs["batch"].apply(lambda x: x.split("d")[0])
    input.obs["cell type"] = solution.obs["cell_type"][input.obs.index]
    del ad_mod1
    del ad_mod2  
    del solution

    sc.pp.normalize_total(input, target_sum=1e4)
    sc.pp.log1p(input)
    input.raw = input

    sc.pp.highly_variable_genes(
        input,
        n_top_genes=4000,
        flavor="seurat_v3",
        batch_key="sample",
        subset=True,
        layer="counts"
    )

    # Setting up data and model
    scvi.model.TOTALVI.setup_anndata(
        input,
        protein_expression_obsm_key="ADT",
        layer="counts",
        batch_key="sample"
    )

    vae = scvi.model.TOTALVI(input, latent_distribution="normal")

    # Training totalVI model
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
    del adata
    del evaluation

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