# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2022/17/05

### Train Winner CITE online model on phase2 full data
# Liam_challenge_reproducibility environment

# Source: https://github.com/openproblems-bio/neurips2021_multimodal_topmethods/blob/main/src/joint_embedding/methods/Guanlab-dengkw/run/script.py

# Added random_state to TruncatedSVD calls

# Imports
import anndata as ad
import numpy as np
import scanpy as sc


from sklearn.decomposition import TruncatedSVD

seeds = [0, 994, 236, 71, 415]
for seed in seeds:

    model = 'Winner_CITE_online_seed_{}'.format(seed)
    print(model)
    
    # Load data
    ad_mod1 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_cite_phase2/openproblems_bmmc_cite_phase2.censor_dataset.output_mod1.h5ad")
    ad_mod2 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_cite_phase2/openproblems_bmmc_cite_phase2.censor_dataset.output_mod2.h5ad")
    solution = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_cite_phase2/openproblems_bmmc_cite_phase2.censor_dataset.output_solution.h5ad")
    
    ad_mod1.obs["sample"] = ad_mod1.obs["batch"]
    ad_mod1.obs["donor"] = ad_mod1.obs["batch"].apply(lambda x: x.split("d")[1])
    ad_mod1.obs["site"] = ad_mod1.obs["batch"].apply(lambda x: x.split("d")[0])
    ad_mod1.obs["cell type"] = solution.obs["cell_type"][ad_mod1.obs.index]

    
    def normalize(arr):
        arr_sd = np.std(arr, axis=1).reshape(-1, 1)
        arr_mean = np.mean(arr, axis=1).reshape(-1, 1)
        return (arr - arr_mean) / arr_sd
    
    # if mod1_type == "GEX" and mod2_type == "ADT":
    n_mod1 = 73
    n_mod2 = 27

    # logging.info('Performing dimensionality reduction on modality 1 values...')
    embedder_mod1 = TruncatedSVD(n_components=n_mod1, random_state=seed)
    mod1_pca = embedder_mod1.fit_transform(ad_mod1.X)
    mod1_obs = ad_mod1.obs
    mod1_uns = ad_mod1.uns

    # logging.info('Performing dimensionality reduction on modality 2 values...')
    embedder_mod1 = TruncatedSVD(n_components=n_mod2, random_state=seed)
    mod2_pca = embedder_mod1.fit_transform(ad_mod2.X)

    # logging.info('Concatenating datasets')
    pca_combined = np.concatenate([normalize(mod1_pca), normalize(mod2_pca)], axis=1)


    # Save embedding
    adata = ad.AnnData(
        X=ad_mod1.X,
        obs=mod1_obs,
        uns={
            'dataset_id': mod1_uns['dataset_id'],
            'method_id': model,
        },
        obsm={'embedding': pca_combined}
    )

    evaluation = ad.AnnData(
        X=pca_combined,
        obs=mod1_obs,
        uns={
            'dataset_id': mod1_uns['dataset_id'],
            'method_id': model,
        }
    )


    adata.write_h5ad("Embeddings/adata_{}.h5ad".format(model), compression="gzip")
    evaluation.write_h5ad("Predictions/{}.prediction.h5ad".format(model), compression="gzip")
    
    # This is probably not necessary, but to be on the safe side delete all variables at the end of each loop to avoid carrying over of information from one run to another
    del ad_mod1
    del ad_mod2
    del solution
    