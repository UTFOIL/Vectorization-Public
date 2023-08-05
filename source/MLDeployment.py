# MLDeployment
# for deploying the the Machine Curator on a vector set SAM 10/5/22

import MLLibrary as mll
import sys
# import keras

# feature__file = r"E:\2P imaging\2021_Chronic_Imaging\Doug week 5\batch_220518-185057\curations\vertices_220518-185057_features_Doug_Fused_Raw_w5.mat"
feature__file = sys.argv[1] # ??? the r"..." formating no longer needed when string is passed from outside ???

# model_file = "laplacian-log_radius.h5"
model_file    = sys.argv[2]

# mymodel = keras.models.load_model( model_file )
VectorDataset2 = mll.DeploymentData( model_file, feature__file )