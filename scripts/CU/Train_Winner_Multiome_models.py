# NeurIPS 2021 Competition- Multimodal Single-Cell Data Integration
# Author: Pia Rautenstrauch
# Created: 2022/17/05

### Train Winner Multiome online model on phase2 full data
# submission_170825 environment

# Source: https://github.com/openproblems-bio/neurips2021_multimodal_topmethods/blob/main/src/joint_embedding/methods/submission_170825/run/script.py

# Imports
import anndata as ad
import pandas as pd
from tensorflow.keras.layers import Input, Dense, Dropout
from tensorflow.keras.layers import concatenate
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow import keras
import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
#from keras import backend as K
from tensorflow.keras.constraints import Constraint
import tensorflow.keras.backend as K

from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model
import warnings
from numpy.random import seed
import tensorflow as tf

print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))


seeds = [0, 994, 236, 71, 415]
for single_seed in seeds:
    seed(single_seed)
    tf.compat.v1.random.set_random_seed(single_seed)

    model = 'Winner_Multiome_online_seed_{}'.format(single_seed)
    print(model)

    ad_mod1 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad")
    ad_mod2 = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2.h5ad")
    solution = ad.read_h5ad("./../../data/original/neurips_competition/openproblems_bmmc_multiome_phase2/openproblems_bmmc_multiome_phase2.censor_dataset.output_solution.h5ad")

    
    ad_mod1.obs["sample"] = ad_mod1.obs["batch"]
    ad_mod1.obs["donor"] = ad_mod1.obs["batch"].apply(lambda x: x.split("d")[1])
    ad_mod1.obs["site"] = ad_mod1.obs["batch"].apply(lambda x: x.split("d")[0])
    ad_mod1.obs["cell type"] = solution.obs["cell_type"][ad_mod1.obs.index]
    
    # high variable gene calculation
    min_cells = int(ad_mod2.shape[0] * 0.03)
    sc.pp.highly_variable_genes(ad_mod1, batch_key ='batch', subset = True)
    sc.pp.filter_genes(ad_mod2, min_cells=min_cells)

    ad_mod_1 = ad_mod1[:, ad_mod1.var.highly_variable]

    ## Convert to  csv for AE training
    scRNAseq1 = ad_mod_1.X.toarray()
    scRNAseq2 = ad_mod2.X.toarray()
    
    

    class WeightsOrthogonalityConstraint(Constraint):
        def __init__(self, encoding_dim, weightage = 1.0, axis = 0):
            self.encoding_dim = encoding_dim
            self.weightage = weightage
            self.axis = axis

        def weights_orthogonality(self, w):
            if(self.axis==1):
                w = K.transpose(w)
            if(self.encoding_dim > 1):
                m = K.dot(K.transpose(w), w) - K.eye(self.encoding_dim)
                return self.weightage * K.sqrt(K.sum(K.square(m)))
            else:
                m = K.sum(w ** 2) - 1.
                return m

        def __call__(self, w):
            return self.weights_orthogonality(w)


    # Input Layer
    ncol_scRNAseq1 = scRNAseq1.shape[1]
    input_dim_scRNAseq1 = Input(shape = (ncol_scRNAseq1, ), name = "scRNAseq1")
    ncol_scRNAseq2 = scRNAseq2.shape[1]
    input_dim_scRNAseq2 = Input(shape = (ncol_scRNAseq2, ), name = "scRNAseq2")

    encoding_dim_scRNAseq1 = 64
    encoding_dim_scRNAseq2 = 64

    dropout_scRNAseq1 = Dropout(0.1, name = "Dropout_scRNAseq1")(input_dim_scRNAseq1)
    dropout_scRNAseq2 = Dropout(0.1, name = "Dropout_scRNAseq2")(input_dim_scRNAseq2)

    encoded_scRNAseq1 = Dense(encoding_dim_scRNAseq1, activation = 'relu', name = "Encoder_scRNAseq1", use_bias=True, kernel_regularizer=WeightsOrthogonalityConstraint(64, weightage=1., axis=0))(dropout_scRNAseq1) #300 #prv 256 
    encoded_scRNAseq2 = Dense(encoding_dim_scRNAseq2, activation = 'relu', name = "Encoder_scRNAseq2", use_bias=True, kernel_regularizer=WeightsOrthogonalityConstraint(64, weightage=1., axis=0))(dropout_scRNAseq2)

    merge = concatenate([encoded_scRNAseq1,  encoded_scRNAseq2])

    bottleneck = Dense(64, kernel_initializer = 'uniform', activation = 'linear', name = "Bottleneck")(merge) #50

    merge_inverse = Dense(encoding_dim_scRNAseq1 + encoding_dim_scRNAseq2, activation = 'relu', name = "Concatenate_Inverse")(bottleneck)

    decoded_scRNAseq1 = Dense(ncol_scRNAseq1, activation = 'relu', name = "Decoder_scRNAseq1")(merge_inverse) #sigmoid

    decoded_scRNAseq2 = Dense(ncol_scRNAseq2, activation = 'relu', name = "Decoder_scRNAseq2")(merge_inverse)

    autoencoder = Model([input_dim_scRNAseq1, input_dim_scRNAseq2],  [decoded_scRNAseq1, decoded_scRNAseq2])

    opt = Adam(lr=0.0001)
    autoencoder.compile(optimizer = opt, loss={'Decoder_scRNAseq1': 'mean_squared_error', 'Decoder_scRNAseq2': 'mean_squared_error'}) #loss_weights = [1., 1.]
    autoencoder.summary()

    es = EarlyStopping(monitor='val_loss', mode='min', verbose=1,patience=20)
    # Autoencoder training
    estimator = autoencoder.fit([scRNAseq1, scRNAseq2], [scRNAseq1, scRNAseq2], epochs = 600, batch_size = 32, validation_split = 0.2, shuffle = True, verbose = 1, callbacks=[es]) #prev 64 BS prev 32


    encoder = Model([input_dim_scRNAseq1, input_dim_scRNAseq2], bottleneck)
    bottleneck_representation = encoder.predict([scRNAseq1, scRNAseq2])

    embd = pd.DataFrame(bottleneck_representation)
    #embd  = scipy.sparse.csr_matrix(RNA_ATAC_Latent.values)

    mod1_obs = ad_mod1.obs
    mod1_uns = ad_mod1.uns

    # Save embedding
    adata = ad.AnnData(
        X=ad_mod1.X,
        obs=mod1_obs,
        uns={
            'dataset_id': ad_mod1.uns['dataset_id'],
            'method_id': model,
        },
        obsm={'embedding' : embd.values}
    )

    evaluation = ad.AnnData(
        X=embd.values,
        obs=mod1_obs,
        uns={
            'dataset_id': ad_mod1.uns['dataset_id'],
            'method_id': model,
        }
    )


    adata.write_h5ad("Embeddings/adata_{}.h5ad".format(model), compression="gzip")
    evaluation.write_h5ad("Predictions/{}.prediction.h5ad".format(model), compression="gzip")
    
    # This is probably not necessary, but to be on the safe side delete all variables at the end of each loop to avoid carrying over of information from one run to another
    del ad_mod1
    del ad_mod2
    del solution
    del model
    del min_cells
    del ad_mod_1
    del scRNAseq1
    del scRNAseq2
    del ncol_scRNAseq1
    del input_dim_scRNAseq1
    del ncol_scRNAseq2
    del input_dim_scRNAseq2
    del dropout_scRNAseq1
    del dropout_scRNAseq2 
    del encoded_scRNAseq1 
    del encoded_scRNAseq2 
    del merge 
    del bottleneck 
    del merge_inverse 
    del decoded_scRNAseq1
    del decoded_scRNAseq2 
    del autoencoder 
    del opt 
    del estimator 
    del encoder 
    del bottleneck_representation 
    del embd
    del mod1_obs 
    del mod1_uns
    del adata 
    del evaluation