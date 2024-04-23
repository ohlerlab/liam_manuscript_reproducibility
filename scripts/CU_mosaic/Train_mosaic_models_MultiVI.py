# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2023/02/21

### Train MultiVI with default params on mosaic use cases based on phase2 full data
# dynamic_Liam_challenge_reproducibility environment


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

# Seeds (either for random data partition (training with random seed 0, or for training))
seeds = [8831, 234, 11, 9631, 94]

for use_case in ['mosaic_full', 'mosaic_a', 'mosaic_b']: 
    for seed in seeds:
        model = 'MultiVI_{}_seed_{}'.format(use_case, seed)
        print(model)
        
        # Load data
        if use_case == 'mosaic_full':
            scvi._settings.ScviConfig(seed=seed)
            # note that in the not subsampled control use case a and b are identical!
            ad_mod1 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod1.h5ad")
            ad_mod2 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod2.h5ad")
        
        elif use_case in ['mosaic_a']:
            scvi._settings.ScviConfig()
            ad_mod1 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod1.h5ad")
            ad_mod2 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod2_10_subsample_seed_{}.h5ad".format(seed))
        
        elif use_case in ['mosaic_b']:
            scvi._settings.ScviConfig()
            ad_mod1 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_b.output_mod1.h5ad")
            ad_mod2 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_b.output_mod2_10_subsample_seed_{}.h5ad".format(seed))
         
        # Setup input object
        input = ad.AnnData(scipy.sparse.hstack([ad_mod1.layers["counts"], ad_mod2.X]).tocsr(), ad_mod1.obs, pd.concat((ad_mod1.var, ad_mod2.var), axis=0))
        del ad_mod1
        del ad_mod2

        input_paired = input[input.obs['mod'] == 'paired'].copy()
        input_rna = input[input.obs['mod'] == 'gex'].copy()
        input_atac = input[input.obs['mod'] == 'atac'].copy()
        
        del input

        input = scvi.data.organize_multiome_anndatas(input_paired, input_rna, input_atac)

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
        vae.save("./../../models/CU_mosaic/{}".format(model), save_anndata=True)

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
