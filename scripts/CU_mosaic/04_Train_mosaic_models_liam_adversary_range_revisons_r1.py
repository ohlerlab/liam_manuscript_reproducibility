# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2023/17/10

### Train liam on mosaic use cases based on phase2 full data
# dynamic_Liam_challenge_reproducibility environment


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

# simple greedy grid search 
adversarial_loss_factors = [1.0, 10.0, 25.0]

for factor in adversarial_loss_factors:
    # Setup param mapping
    model_param_mapping = {}
    model_param_mapping['setup_anndata_params'] = {'chrom_acc_obsm_key': 'ATAC', 'batch_key': 'sample'}
    model_param_mapping['Liam_params'] = {'adversarial_training': True, 'n_latent': 20, 'factor_adversarial_loss': factor}

    # Seeds (either for random data partition (training with random seed 0, or for training))
    seeds = [8831, 234, 11, 9631, 94]

    for use_case in ['mosaic_full', 'mosaic_a', 'mosaic_b']:
        for seed in seeds:
            model = 'liam_{}_x{}_seed_{}'.format(use_case, factor, seed)
            print(model)

            # Load data
            if use_case == 'mosaic_full':
                scvi._settings.ScviConfig(seed=seed)
                # note that in the not subsampled control use case a and b are identical!
                ad_mod1 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod1.h5ad")
                ad_mod2 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod2.h5ad")

            elif use_case in ['mosaic_a', 'mosaic_atac_a', 'mosaic_rest_a']:
                scvi._settings.ScviConfig()
                ad_mod1 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod1.h5ad")
                ad_mod2 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_a.output_mod2_10_subsample_seed_{}.h5ad".format(seed))

            elif use_case in ['mosaic_b', 'mosaic_atac_b', 'mosaic_rest_b']:
                scvi._settings.ScviConfig()
                ad_mod1 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_b.output_mod1.h5ad")
                ad_mod2 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_b.output_mod2_10_subsample_seed_{}.h5ad".format(seed))

            elif use_case == 'paired_full':
                scvi._settings.ScviConfig(seed=seed)
                ad_mod1 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_full.output_mod1.h5ad")
                ad_mod2 = ad.read_h5ad("./../../data/derived/neurips_competition/mosaic_use_case_full.output_mod2.h5ad")

            # Setup input object
            input = ad_mod1.copy()
            input.X = input.layers["counts"].tocsr()
            input.obsm["ATAC"] = ad_mod2.X.tocsr()
            del ad_mod1
            del ad_mod2

            if 'mosaic_atac' in use_case:
                input = input[(input.obs['site'] == 's1') | (input.obs['site'] == 's3') ].copy()

            if 'mosaic_rest' in use_case:
                input = input[(input.obs['site'] == 's2') | (input.obs['site'] == 's4') ].copy()
            # Setting up data.
            liam_NeurIPS2021_challenge_reproducibility.Liam.setup_anndata(
                input,
            **model_param_mapping['setup_anndata_params'],
            )

            vae = liam_NeurIPS2021_challenge_reproducibility.Liam(input, **model_param_mapping['Liam_params'])


            # Training Liam models
            vae.train(train_size=0.95, validation_size=0.05,
                          batch_size=128, early_stopping=True, save_best=True, early_stopping_patience=10)

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
