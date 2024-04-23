# Evaluate all models for Multiome data mosaic use case
# Author: Pia Rautenstrauch
# Created: 2023/02/26

# Environment: scib_v1.0.1_min
# Code adapted from https://github.com/openproblems-bio/neurips2021_multimodal_viash/tree/main/src/joint_embedding/metrics 
# folders with names of respective metrics, scripts called script.py

# Imports
import os
import anndata as ad
import numpy as np
import pandas as pd
import scib
import scanpy as sc
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

scib.__version__

# List of identifiers of trained models
models = []

seeds = [8831, 234, 11, 9631, 94]

for use_case in [#'mosaic_full', 'mosaic_a', 'mosaic_b', 
    'mosaic_atac_a', 'mosaic_rest_a', 'mosaic_atac_b', 'mosaic_rest_b', 'paired_full']:
    for seed in seeds:
        model = 'liam_{}_seed_{}'.format(use_case, seed)
        models += [model]
        
for use_case in ['mosaic_full', 'mosaic_a', 'mosaic_b']: 
    for seed in seeds:
        model = 'MultiVI_{}_seed_{}'.format(use_case, seed)
        models += [model]


# Collect computed scores, nested dict is simple to convert to pd.DataFrame
score_dict = {}
for model in models:
    np.random.seed(61)

    # Initialize nested dict
    score_dict[model] = {}
    
    # Read embeddings
    embedding = ad.read_h5ad("./Predictions/{}.prediction.h5ad".format(model))

    # Remove _paired suffix (from MultiVI models)
    if "MultiVI" in model:
        embedding.obs.index = embedding.obs.index.map(lambda x: x.split('_')[0])

    # How many dimensions does embedding have
    score_dict[model]['dims'] = embedding.X.shape[1]
    
    # Compute neighbors 
    embedding.obsm['X_emb'] = embedding.X.copy()

    sc.pp.neighbors(embedding, use_rep='X_emb')
    sc.tl.umap(embedding)
    
    # Load metadata (cell type information, etc.)
    solution = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_solution.h5ad")
    
    # Make sure order is the same as for the embedding
    solution = solution[embedding.obs.index,:].copy()
    
    # Transfer batch information and cell type labels to embedding as done for neurips evaluation
    embedding.obs['sample'] = solution.obs['batch']
    embedding.obs['site'] = solution.obs['batch'].apply(lambda x: x.split("d")[0])
    embedding.obs['donor'] = solution.obs['batch'].apply(lambda x: x.split("d")[1])
    solution.obs['sample'] = solution.obs['batch']
    embedding.obs['cell type'] = solution.obs['cell_type']
    
    ### Rename cell_type for solution
    solution.obs.rename({'cell_type': 'cell type'}, axis=1, inplace=True)
    
    # make sure these variables are categorical
    embedding.obs['sample'] = embedding.obs['sample'].astype('category')
    embedding.obs['site'] = embedding.obs['site'].astype('category')
    embedding.obs['donor'] = embedding.obs['donor'].astype('category')
    solution.obs['sample'] = embedding.obs['sample'].astype('category')
    embedding.obs['cell type'] = embedding.obs['cell type'].astype('category')
    solution.obs['cell type'] = solution.obs['cell type'].astype('category')

    
    # Compute scores
    ## Level of evaluation: batch/sample
    ### graph iLISI and cLISI on variable batch
    score_dict[model]['iLISI_batch'], score_dict[model]['cLISI_full'] =  scib.me.lisi.lisi_graph(embedding, batch_key='sample', label_key='cell type', multiprocessing=True)
    
    ### asw_batch
    score = scib.me.silhouette_batch(
        embedding,
        batch_key='sample',
        group_key='cell type',
        embed='X_emb',
        verbose=False
    )

    score_dict[model]['asw_batch'] = score
    
    ### asw_label
    score = scib.me.silhouette(
        embedding,
        group_key='cell type',
        embed='X_emb'
    )
    
    score_dict[model]['asw_label'] = score

    ### cc_cons
    organism = solution.uns['organism']
    
    score = scib.me.cell_cycle(
        adata_pre=solution,
        adata_post=embedding,
        batch_key='sample',
        embed='X_emb',
        recompute_cc=True,
        organism=organism
    )
    
    score_dict[model]['cc_cons'] = score

    ### graph_conn
    score = scib.me.graph_connectivity(
        embedding,
        label_key='cell type'
    )

    score_dict[model]['graph_conn'] = score
      
    ### nmi
    scib.cl.opt_louvain(
        embedding,
        label_key='cell type',
        cluster_key='cluster',
        plot=False,
        inplace=True,
        force=True
    )

    score = scib.me.nmi(
        embedding,
        group1='cluster',
        group2='cell type'
    )
    
    score_dict[model]['nmi'] = score

    ### ti_cons_batch_mean
    adt_atac_trajectory = 'pseudotime_order_ATAC' if 'pseudotime_order_ATAC' in solution.obs else 'pseudotime_order_ADT'
    obs_keys = solution.obs_keys()

    if 'pseudotime_order_GEX' in obs_keys:
        score_rna = scib.me.trajectory_conservation(
            adata_pre=solution,
            adata_post=embedding,
            label_key='cell type',
            batch_key='sample',
            pseudotime_key='pseudotime_order_GEX'
        )
    else:
        score_rna = np.nan

    if adt_atac_trajectory in obs_keys:
        score_adt_atac = scib.me.trajectory_conservation(
            adata_pre=solution,
            adata_post=embedding,
            label_key='cell type',
            batch_key='sample',
            pseudotime_key=adt_atac_trajectory
        )
    else:
        score_adt_atac = np.nan

    score_mean = (score_rna + score_adt_atac) / 2
    score_dict[model]['ti_cons_batch_gex'] = score_rna
    score_dict[model]['ti_cons_batch_adt_atac'] = score_adt_atac
    score_dict[model]['ti_cons_batch_mean'] = score_mean


    ## Level of evaluation: site
    ### graph iLISI and cLISI on variable site
    score_dict[model]['iLISI_site'] =  scib.me.lisi.ilisi_graph(embedding, batch_key='site', multiprocessing=True)

    ### asw_site
    score = scib.me.silhouette_batch(
        embedding,
        batch_key='site',
        group_key='cell type',
        embed='X_emb',
        verbose=False
    )

    score_dict[model]['asw_site'] = score
    
    
    ## Level of evaluation: modality (mod)
    ### graph iLISI on variable modality (mod)
    if len(embedding.obs['mod'].unique()) > 1:
        score_dict[model]['iLISI_modality'] =  scib.me.lisi.ilisi_graph(embedding, batch_key='mod', multiprocessing=True)

        ### asw_batch_modality
        score = scib.me.silhouette_batch(
            embedding,
            batch_key='mod',
            group_key='cell type',
            embed='X_emb',
            verbose=False
        )

        score_dict[model]['asw_batch_modality'] = score

        ### Save UMAP visualization as png 
        sc.pl.umap(
            embedding,
            color=['cell type', 'sample', 'site', 'mod'],
            size=3,
            frameon=False,
            ncols=4,
            legend_loc=None,
            wspace=0.0,
            save="_{}_multiome.png".format(model)
        )
    
    else:
        score_dict[model]['iLISI_modality'] = np.nan
        score_dict[model]['asw_batch_modality'] = np.nan
        
        ### Save UMAP visualization as png 
        sc.pl.umap(
            embedding,
            color=['cell type', 'sample', 'site', 'mod'],
            size=3,
            frameon=False,
            ncols=4,
            legend_loc=None,
            wspace=0.0,
            save="_{}_multiome.png".format(model)
        )
    

    ## Level of evaluation: donor 1 (d1) 
    ### Subset data to donor 1
    subset = embedding[embedding.obs['donor'] == "1", :].copy()
    if len(subset.obs['sample'].unique()) > 1:
        ### graph iLISI and cLISI on batch of donor 1 only 
        #### compute new neighborhood graph on only subset of data
        sc.pp.neighbors(subset, use_rep="X_emb")
        score_dict[model]['iLISI_d1'], score_dict[model]['cLISI_d1'] =  scib.me.lisi.lisi_graph(subset, batch_key='sample', label_key='cell type',multiprocessing=True)

        ### asw_d1
        score = scib.me.silhouette_batch(
            subset,
            batch_key='sample',
            group_key='cell type',
            embed='X_emb',
            verbose=False
        )

        score_dict[model]['asw_batch_d1'] = score

        ### Save UMAP visualization
        sc.pl.umap(
            subset,
            color=['cell type', 'sample', 'mod'],
            size=3,
            frameon=False,
            ncols=3,
            legend_loc=None,
            wspace=0.0,
            save="_{}_d1_subset_multiome.png".format(model)
        )
    else:
        score_dict[model]['iLISI_d1'], score_dict[model]['cLISI_d1'] = np.nan, np.nan
        score_dict[model]['asw_batch_d1'] = np.nan
        
        ### Save UMAP visualization
        sc.pl.umap(
            subset,
            color=['cell type', 'sample', 'mod'],
            size=3,
            frameon=False,
            ncols=3,
            legend_loc=None,
            wspace=0.0,
            save="_{}_d1_subset_multiome.png".format(model)
        )
   

    pd.DataFrame(score_dict).to_csv("Evaluation/scores/Multiome/batch_removal_scores_{}.csv".format(model), index=True)

# Store legend once for figures   
sc.pl.umap(
    embedding,
    color=['cell type', 'sample', 'site', 'mod'],
    size=3,
    frameon=False,
    ncols=4,
    wspace=1,
    save="_{}_legend_multiome.png".format(model)
)

sc.pl.umap(
        subset,
        color=['cell type', 'sample', 'mod'],
        size=3,
        frameon=False,
        ncols=3,
        wspace=1,
        save="_{}_d1_subset_legend_multiome.png".format(model)
)
    

pd.DataFrame(score_dict).to_csv("Evaluation/scores/Multiome/batch_removal_scores_all_models.csv", index=True)
