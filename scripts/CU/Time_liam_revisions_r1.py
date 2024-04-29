# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2023/26/10

### Train models on phase2 full data (Liam)
#### BAVAE - sample

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

# Import for timing
from timeit import default_timer as timer
from datetime import timedelta

print(torch.cuda.is_available())

# Setup model to param mapping
model_param_mapping = {}

model_param_mapping['BAVAE_sample_100_seed_0_timing_t4'] = {}
model_param_mapping['BAVAE_sample_100_seed_0_timing_t4']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'sample'}
model_param_mapping['BAVAE_sample_100_seed_0_timing_t4']['Liam_params'] = {'adversarial_training': True, 'n_latent': 20}

# Train models
for model in model_param_mapping.keys():
    print(model)
    scvi._settings.ScviConfig()
    
    # Load data
    ad_mod1 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad")
    ad_mod2 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2.h5ad")
    solution = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_solution.h5ad")


    # Setup input object
    input = ad_mod1.copy()
    input.X = input.layers["counts"].tocsr()
    input.obsm["ATAC"] = ad_mod2.X.tocsr()
    input.obs["sample"] = input.obs["batch"]
    input.obs["donor"] = input.obs["batch"].apply(lambda x: x.split("d")[1])
    input.obs["site"] = input.obs["batch"].apply(lambda x: x.split("d")[0])
    input.obs["cell type"] = solution.obs["cell_type"][input.obs.index]
    del ad_mod1
    del ad_mod2
    del solution


    # Setting up data.
    liam_NeurIPS2021_challenge_reproducibility.Liam.setup_anndata(
            input,
        **model_param_mapping[model]['setup_anndata_params'],
        )

    vae = liam_NeurIPS2021_challenge_reproducibility.Liam(input, **model_param_mapping[model]['Liam_params'])


    # Training Liam models
    start = timer()
    
    vae.train(train_size=0.95, validation_size=0.05,
                  batch_size=128, early_stopping=True, save_best=True, early_stopping_patience=10)
    
    end = timer()
    print('Wallclock training run time.')
    print(timedelta(seconds=end-start))
    print()
    
    print('{}: {}'.format(model, timedelta(seconds=end-start)), file=open('Time_Wallclock_training_run_time_revisions_r1.txt', 'a'))

    input.obsm["embedding"] = vae.get_latent_representation()

    # Save trained model
    vae.save("Models/{}".format(model), save_anndata=True)



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