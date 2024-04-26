# Extended treatment-control use case 
### Save both full models and predictions in format as required for NeurIPS data
### Environment: dynamic_LIAM_challenge_reproducibility
### Author: Pia Rautenstrauch
### Date: 24th of January 2024

# Imports
import liam_NeurIPS2021_challenge_reproducibility
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
import anndata as ad
import scvi
from sklearn.metrics import silhouette_score, silhouette_samples

torch.cuda.is_available()

def assign_rep(x):
    if "DIG" in x:
        return "1"
    elif "LLL" in x:
        return "2"
    elif "10k" in x:
        return "3"
    else:
        return "4"

model_param_mapping = {}

model_param_mapping['BAVAE_sample_100_extended_use_case_x10'] = {}
model_param_mapping['BAVAE_sample_100_extended_use_case_x10']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'batch'}
model_param_mapping['BAVAE_sample_100_extended_use_case_x10']['Liam_params'] = {'adversarial_training': True, 'n_latent': 20, 'factor_adversarial_loss': 10.0}

model_param_mapping['BAVAE_sample_100_extended_use_case_x5'] = {}
model_param_mapping['BAVAE_sample_100_extended_use_case_x5']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'batch'}
model_param_mapping['BAVAE_sample_100_extended_use_case_x5']['Liam_params'] = {'adversarial_training': True, 'n_latent': 20, 'factor_adversarial_loss': 5.0}

model_param_mapping['BAVAE_sample_100_extended_use_case'] = {}
model_param_mapping['BAVAE_sample_100_extended_use_case']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'batch'}
model_param_mapping['BAVAE_sample_100_extended_use_case']['Liam_params'] = {'adversarial_training': True, 'n_latent': 20}

# also model ATAC lib size (but batch independent)
model_param_mapping['VAE_100_extended_use_case'] = {}
model_param_mapping['VAE_100_extended_use_case']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': None}
model_param_mapping['VAE_100_extended_use_case']['Liam_params'] = {'n_latent': 20, 'dispersion_gex': 'gene', 'dispersion_atac': 'constant'}

model_param_mapping['BAVAE_sample_100_extended_use_case_x50'] = {}
model_param_mapping['BAVAE_sample_100_extended_use_case_x50']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'batch'}
model_param_mapping['BAVAE_sample_100_extended_use_case_x50']['Liam_params'] = {'adversarial_training': True, 'n_latent': 20, 'factor_adversarial_loss': 50.0}

model_param_mapping['BAVAE_sample_100_extended_use_case_x25'] = {}
model_param_mapping['BAVAE_sample_100_extended_use_case_x25']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'batch'}
model_param_mapping['BAVAE_sample_100_extended_use_case_x25']['Liam_params'] = {'adversarial_training': True, 'n_latent': 20, 'factor_adversarial_loss': 25.0}

model_param_mapping['BAVAE_sample_100_extended_use_case_x100'] = {}
model_param_mapping['BAVAE_sample_100_extended_use_case_x100']['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'batch'}
model_param_mapping['BAVAE_sample_100_extended_use_case_x100']['Liam_params'] = {'adversarial_training': True, 'n_latent': 20, 'factor_adversarial_loss': 100.0}

for model in model_param_mapping.keys():
    print(model)
    scvi._settings.ScviConfig()
    input = ad.read_h5ad("./../../data/derived/Mimitou2021/DOGMA_seq/extended_treatment_control_use_case_revisions_r1.h5ad")

    input.obs["sample"] = input.obs["batch"]
    input.obs['replicate'] = input.obs['batch'].apply(lambda x: assign_rep(x))
    
    liam_NeurIPS2021_challenge_reproducibility.Liam.setup_anndata(
            input,
            **model_param_mapping[model]['setup_anndata_params']
        )

    vae = liam_NeurIPS2021_challenge_reproducibility.Liam(input, **model_param_mapping[model]['Liam_params'])

    vae.train(train_size=0.95, validation_size=0.05,
                  batch_size=128, early_stopping=True, save_best=True, early_stopping_patience=10)

    input.obsm["embedding"] = vae.get_latent_representation()

    sc.pp.neighbors(input, use_rep="embedding")
    sc.tl.umap(input)
    sc.tl.leiden(input, key_added="leiden_embedding")    

    vae.save("./../../models/ETCU/{}".format(model), save_anndata=True)
    
    del vae
    del input
    