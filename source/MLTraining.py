import MLLibrary as mll
import sys

# SAM 8/19/2022

# hardcoded base features # ???? make this an input as well ??????
base_feature_names=["laplacian","log_radius"] # base features # ~!!!!! any other choice of features triggers a restriction of the training data to just the vertices labeled as foreground by the base network

# # # feature_names=["laplacian","log_radius","gradient_1","gradient_2"] # base + gradients1-2
# # foreground_feature_names = base_feature_names + ["gradient_1","vesselness"] # !!!!!! make this input
# # foreground_feature_names = ["log_radius","gradient_1"] AUC: 0.50 ?? rerun ??
# # foreground_feature_names = ["gradient_1"] AUC: 0.52
# # foreground_feature_names = ["gradient_1","vesselness"] AUC: 0.68
# # foreground_feature_names = ["vesselness"] # AUC: 0.68
# # foreground_feature_names = ["vesselness","log_radius"] # AUC: 0.69
# foreground_feature_names = ["gradient__"] # AUC = 0.55
# # foreground_feature_names = ["gradient_N"] # AUC = 0.59
# # foreground_feature_names = ["gradient_N","vesselness"] # 0.69
# # foreground_feature_names = ["laplacian","log_radius","gradient_N","vesselness"]  # 0.75
# foreground_feature_names = ["gradient__","log_radius"] # T1: AUC = 0.66 T2: 0.77
# foreground_feature_names = ["gradient__","log_radius","vesselness"] # T2: AUC = 0.77
# foreground_feature_names = ["gradient__","log_radius","blobness"] # T2: AUC = 0.81
# foreground_feature_names = ["gradient__","log_radius","blobness","vesselness"] # T2: AUC = 0.80
# foreground_feature_names = ["gradient_N","log_radius","blobness"] # T2: AUC = 0.82 ************ 3 feature winner
# foreground_feature_names = ["gradient_N","log_radius","blobness","laplacian"] # T2: AUC = 0.82
# foreground_feature_names = ["gradient__","log_radius","blobness","laplacian"] # T2: AUC = 0.83
foreground_feature_names = ["gradient_N","log_radius","blobness_1","blobness_2"] # T2: AUC =

# foreground_feature_names = base_feature_names + sys.argv[1]

# # Training Set # 2 
# curation_file_1 = r"E:\2P imaging\2021_Chronic_Imaging\Brett week 2\New folder\Truncated\batch_230109-123958\curations\vertices_230109-123958_Brett_Fused_Raw_w2"
# curation_file_2 = r"E:\2P imaging\2021_Chronic_Imaging\Brett week 3\truncated\batch_221106-224630\curations\vertices_221106-224630_Brett_Fused_Raw_w3"
# curation_file_3 = r"E:\2P imaging\2021_Chronic_Imaging\Brett week 4\batch_221018-164034\curations\vertices_221018-164034_Brett_Fused_Raw_w4"

# Training Set # 3
curation_file_1 = r"E:\Anna\batch_230707-232540\curations\vertices_230707-232540_block_0"
curation_file_2 = r"E:\Anna\batch_230708-142756\curations\vertices_230708-142756_block_0"

# # # !!!!!! make this input
# curation_file_1 = r"E:\2P imaging\2021_Chronic_Imaging\Doug week 5\batch_220518-185057\curations\vertices_220518-185057_Doug_Fused_Raw_w5.mat"
# curation_file_2 = r"E:\2P imaging\2021_Chronic_Imaging\Doug week 3\batch_220511-223915\curations\vertices_220511-223915_Doug_Fused_Raw_w3.mat"
# # # curation_file_3 = r"E:\2P imaging\2021_Chronic_Imaging\Doug week 1\batch_220126-151318\curations\vertices_220126-151318_Doug_Fused_Raw_w1.mat"

# curation_file = r"E:\2P imaging\2021_Chronic_Imaging\Doug week 5\batch_220518-185057\curations\vertices_220518-185057_Doug_Fused_Raw_w5.mat"
# # # feature__file = r"E:\2P imaging\2021_Chronic_Imaging\Doug week 5\batch_220518-185057\curations\vertices_220518-185057_features_Doug_Fused_Raw_w5.mat"
# # curation_file = sys.argv[2]

 #VectorDataset1=TrainingData("/content/vertices_220511-223915_Doug_Fused_Raw_w3.mat","/content/vertices_220511-223915_features_Doug_Fused_Raw_w3.mat")

# training_set_name = 'TrainingSet1'
# training_set_name = 'TrainingSet2'
training_set_name = 'TrainingSet3'

# mll.GenerateMachineCurator(      base_feature_names,training_set_name,curation_file_1,curation_file_2,curation_file_3) # etc. up to curation_file_N
# mll.GenerateMachineCurator(foreground_feature_names,training_set_name,curation_file_1,curation_file_2,curation_file_3) 

# mll.GenerateMachineCurator(      base_feature_names,training_set_name,curation_file_1,curation_file_2) 
mll.GenerateMachineCurator(foreground_feature_names,training_set_name,curation_file_1,curation_file_2) 

# mll.GenerateMachineCurator(      base_feature_names,training_set_name,curation_file_1) 
# mll.GenerateMachineCurator(foreground_feature_names,training_set_name,curation_file_1) 
