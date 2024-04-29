# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2022/16/06

### Train models on phase2 full data (Liam)
#### BAVAE ATAC 10% subsample - ATAC only  - sample
#### BAVAE ATAC 10% subsample - sample 
#### BAVAE concat ATAC 10% subsample - sample 

#### BAVAE ATAC 25% subsample - ATAC only  - sample
#### BAVAE ATAC 25% subsample - sample 
#### BAVAE concat ATAC 25% subsample - sample 


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

# Train 'de novo' models
# Setup model to param mapping
model_param_mapping = {}

seeds = [8831, 234, 11, 9631, 94]

for seed in seeds:
    # Note that this seed is not the random seed for the model training, but the random seed for the ATAC data subsampling!
    model_param_mapping['BAVAE_sample_10_atac_only_seed_{}'.format(seed)] = {}
    model_param_mapping['BAVAE_sample_10_atac_only_seed_{}'.format(seed)]['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'sample'}
    model_param_mapping['BAVAE_sample_10_atac_only_seed_{}'.format(seed)]['Liam_params'] = {'adversarial_training': True, 'n_latent': 10, 'atac_only': True}
    model_param_mapping['BAVAE_sample_10_atac_only_seed_{}'.format(seed)]['path_to_data'] = './../../data/derived/neurips_competition/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2_10_subsample_seed_{}.h5ad'.format(seed)

    model_param_mapping['BAVAE_sample_10_seed_{}'.format(seed)] = {}
    model_param_mapping['BAVAE_sample_10_seed_{}'.format(seed)]['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'sample'}
    model_param_mapping['BAVAE_sample_10_seed_{}'.format(seed)]['Liam_params'] = {'adversarial_training': True, 'n_latent': 20}
    model_param_mapping['BAVAE_sample_10_seed_{}'.format(seed)]['path_to_data'] = './../../data/derived/neurips_competition/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2_10_subsample_seed_{}.h5ad'.format(seed)
    
    model_param_mapping['BAVAE_sample_25_atac_only_seed_{}'.format(seed)] = {}
    model_param_mapping['BAVAE_sample_25_atac_only_seed_{}'.format(seed)]['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'sample'}
    model_param_mapping['BAVAE_sample_25_atac_only_seed_{}'.format(seed)]['Liam_params'] = {'adversarial_training': True, 'n_latent': 10, 'atac_only': True}
    model_param_mapping['BAVAE_sample_25_atac_only_seed_{}'.format(seed)]['path_to_data'] = './../../data/derived/neurips_competition/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2_25_subsample_seed_{}.h5ad'.format(seed)

    model_param_mapping['BAVAE_sample_25_seed_{}'.format(seed)] = {}
    model_param_mapping['BAVAE_sample_25_seed_{}'.format(seed)]['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'sample'}
    model_param_mapping['BAVAE_sample_25_seed_{}'.format(seed)]['Liam_params'] = {'adversarial_training': True, 'n_latent': 20}
    model_param_mapping['BAVAE_sample_25_seed_{}'.format(seed)]['path_to_data'] = './../../data/derived/neurips_competition/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2_25_subsample_seed_{}.h5ad'.format(seed)


# Train models
for model in model_param_mapping.keys():
    print(model)
    scvi._settings.ScviConfig()
    
    # Load data
    ad_mod1 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad")
    ad_mod2 = ad.read_h5ad(model_param_mapping[model]['path_to_data'])
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
    vae.train(train_size=0.95, validation_size=0.05,
                  batch_size=128, early_stopping=True, save_best=True, early_stopping_patience=10)

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

    
# Concatenate models
# Setup model to param mapping
model_param_mapping = {}

seeds = [8831, 234, 11, 9631, 94]

for seed in seeds:
    model_param_mapping['BAVAE_sample_10_concat_seed_{}'.format(seed)] = ('BAVAE_sample_100_rna_only_seed_0', 'BAVAE_sample_10_atac_only_seed_{}'.format(seed))
    
# Train models
for model in model_param_mapping.keys():
    print(model)
    scvi._settings.ScviConfig()
    
    vae_rna = liam_NeurIPS2021_challenge_reproducibility.Liam.load("Models/{}".format(model_param_mapping[model][0]))
    vae_atac = liam_NeurIPS2021_challenge_reproducibility.Liam.load("Models/{}".format(model_param_mapping[model][1]))  
    
    input = vae_atac.adata
    input.obsm["embedding_atac"] = vae_atac.get_latent_representation()
    input.obsm["embedding_gex"] = vae_rna.get_latent_representation()
    
    # Concatenate embeddings
    input.obsm["embedding"] = np.concatenate((input.obsm["embedding_atac"], input.obsm["embedding_gex"]), axis=1)
      
    # Save adata for concat model
    os.mkdir("Models/{}".format(model))
    input.write("Models/{}/adata.h5ad".format(model))
    
  
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
    
    del vae_rna
    del vae_atac
    del input
    del adata

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