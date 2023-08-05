# MLLibrary
# Machine curation library SAM 10/6/22


# hardcode logradius and laplacian features and have model predict on these features.


#when predictions are written, have some record of which features were used to train, and which dataset were used to train. perhaps title.

#try deploying model on features from a different dataset.

#try the first three derivative features, you always want laplacian, try perhaps shape factor feature (vesselness).

# !!!!! these are inputs now for the training file !!!!!
# !wget https://www.dropbox.com/s/6o92zh6tni8fk7q/vertices_220511-223915_features_Doug_Fused_Raw_w3.mat
# !wget https://www.dropbox.com/s/kwg5l0brcc290sy/vertices_220511-223915_Doug_Fused_Raw_w3.mat

# !wget https://www.dropbox.com/s/saop8zw3zc35a0o/vertices_220126-151318_Doug_Fused_Raw_w1.mat
# !wget https://www.dropbox.com/s/g1nze20vcjmv83u/vertices_220126-151318_features_Doug_Fused_Raw_w1.mat

# #USING CORRECT .mat file
# #using extra features
# !wget https://www.dropbox.com/s/uix7f9imyiusi3c/vertices_220518-185057_Doug_Fused_Raw_w5.mat

# !wget https://www.dropbox.com/s/v10t1cf5dj3483b/vertices_220518-185057_features_Doug_Fused_Raw_w5.mat

import numpy as np
import sys
sys.path.append("/home/me/mypy")
import os
import scipy
import math 
import statistics
import keras
import pandas as pd
from keras import layers 


from keras.models import Sequential
from keras.layers import Dense, Activation,Dropout
import tensorflow as tf

import matplotlib.cm
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn import metrics

from scipy import stats
import scipy.io

#import matlab.engine
#eng = matlab.engine.start_matlab()
#tf = eng.isprime(37)
#print(tf)


model_directory = "MLModels"

class DeploymentData:

  def __init__(self,model_file,feature_file):
  # def predictonfeatures(self,model_file,feature_file):

    # self.model_file    = model_file
    # self.feature_file  = feature_file
    # self.feature_names = feature_names

    mymodel = keras.models.load_model(model_file)

    feature_names = mymodel.name.split('-') # feature_names[0] = trainnig_set_name
    feature_names.pop(0) 

    #print(self._feature_dict.keys())
    feature_file_name =             feature_file.split(os.sep)[  -1].split(".")[0] # ??? will this line work for other OS's >>> replace "/" with a call to the system "file seperator" parameter
    # feature_file_path=feature_file.split(os.sep)[0]
    # for directory in feature_file.split(os.sep)[1:-1]:
    #     feature_file_path += os.sep + directory
    feature_file_path = os.sep.join(feature_file.split(os.sep)[0:-1]          )

    feature_dict = scipy.io.loadmat(feature_file)

    prediction_directory = feature_file_path + os.sep + feature_file_name + "_" + model_directory

    if not os.path.isdir(prediction_directory):
      os.mkdir(prediction_directory)
    # prediction_file_path = feature_file_path + os.sep + feature_file_name + "_" + mymodel.name + "_predictions.mat"
    prediction_file_path = prediction_directory + os.sep + mymodel.name + "_predictions.mat"

    feature_array = np.zeros((len(feature_dict[feature_names[0]]),
                              len(             feature_names    ) ))

    feature_idx = -1

    for feature_name in feature_names:
      feature_idx += 1
      feature_array[:,feature_idx]=np.asarray(feature_dict[feature_name]).squeeze()

    # feature_array = np.asarray(feature_dict[feature_names[0]]) # ?????? how to improve the readability of this FOR loop ????
    # for feature_name in feature_names[1:]:
    #   feature_array = np.concatenate((feature_array,np.asarray(feature_dict[feature_name])), axis=1 )
    # # feature_array = np.concatenate((feature_array,np.asarray(feature_dict[feature_names])), axis=1)

    # laplacian=np.asarray(feature_dict["laplacian"])
    # logradius=np.asarray(feature_dict["log_radius"])
    # gradient_1=np.asarray(feature_dict["gradient_1"])
    # gradient_2=np.asarray(feature_dict["gradient_2"])

    # feature_array=np.concatenate((laplacian,logradius,gradient_1, gradient_2),axis=1)
    predictions = mymodel.predict(feature_array)

    # self._feature_dict["predictions"]=predictions.squeeze()
    scipy.io.savemat( prediction_file_path, {'predictions':predictions})

class TrainingData:

  def __init__(self,curation_file):
    # args is one or more curation file paths
    
    # self.curation_dict = { }
    # self._feature_dict = { }
    # self.foregrnd_dict = { }
    # for curation_file in curation_files:

    file_path_directories = curation_file.split(os.sep)
    
    file_path_stem = os.sep.join(file_path_directories[0:-1])
    file_name      =             file_path_directories[  -1].split(".")[0]

    # # 'vertices_' = file_name[0:9]
    # timeStamp = file_name[9 :22]#9+13
    name_prefix = file_name[ 0:22]
    name_suffix = file_name[23:  ]
    
    file_path_stem += os.sep + name_prefix
    
    # curation_file = file_path_stem +     "_"    + name_suffix
    _feature_file = file_path_stem +  "_features_"  + name_suffix
    foregrnd_file = file_path_stem + "_foreground_" + name_suffix

    self.curation_dict = scipy.io.loadmat(curation_file)
    self._feature_dict = scipy.io.loadmat(_feature_file)
    self.foregrnd_dict = scipy.io.loadmat(foregrnd_file)

    # temp_curation_dict = scipy.io.loadmat(curation_file)
    # temp__feature_dict = scipy.io.loadmat(_feature_file)
    # temp_foregrnd_dict = scipy.io.loadmat(foregrnd_file)

    # for                  key  in temp_curation_dict.keys():
    #   self.curation_dict[key] += temp_curation_dict[key]

    # for                  key  in temp__feature_dict.keys():
    #   self._feature_dict[key] += temp__feature_dict[key]
      
    # for                  key  in temp_foregrnd_dict.keys():
    #   self.foregrnd_dict[key] += temp_foregrnd_dict[key]
    
    # self.train_name = file_name
    # # self.train_name   =                  curation_file.split(os.sep)[-1].split(".")[0] # train name should be a nickname for many datasets, input the nickname, keep track of all the files contributing to the training in a list of strings
    # # # self.train_name   =                  curation_file.split("/")[-1].split(".")[0] # train name should be a nickname for many datasets, input the nickname, keep track of all the files contributing to the training in a list of strings
          
  def get_features(self):
      # print( "Curator Data Names:",list(self.mat.keys()))
      print( "Feature Data Column Headers:",list(self._feature_dict.keys()))

  def ChooseFeatures(self, feature_names ):
      curator_label_names = ["displayed_vertices","true_vertices"]

      curator_data_names  = ["vertex_scale_subscripts",
                             "vertex_energies",
                             "vertex_space_subscripts"]+curator_label_names

      self.is_training_base_model = '-'.join(feature_names) == 'laplacian-log_radius' # base model features

      max_num_of_added_vertices=self.curation_dict["max_number_of_added_vertices"][0][0]
      
      # if self.is_training_base_model:
      #   chosen_vertex_indices  =self._feature_dict["chosen_vertex_indices"].squeeze() # use all the vertices (error type dominated by background/low-contrast vertices)
      #                                                                                 # remove indices that did not make it into the curator because of some initial volume exclusion algorithm
      # else:
      #   # chosen_vertex_indices        =np.subtract(chosen_vertex_indices        ,np.ones(len(chosen_vertex_indices        ))).astype(int)
      #   # 
      #   # chosen_vertex_indices_curator = self.foregrnd_dict["chosen_vertex_indices"].squeeze()[chosen_vertex_indices] #  restricted to just the vertices that made it into the curator
      #   chosen_vertex_indices  =self.foregrnd_dict["chosen_vertex_indices"].squeeze()                        # pull chosen vertices from just those labeled as foreground by the base network
        
      chosen_vertex_indices  =self._feature_dict["chosen_vertex_indices"].squeeze() 

        # chosen_vertex_indices_curator=np.subtract(chosen_vertex_indices_curator,np.ones(len(chosen_vertex_indices_curator))).astype(int) # subtract one for MATLAB to Python indexing conversion
      chosen_vertex_indices          =np.subtract(chosen_vertex_indices        ,np.ones(len(chosen_vertex_indices        ))).astype(int) # ""                                                 ""

      self.feature_names      =      feature_names # variable: set by network architect
      self.curator_data_names = curator_data_names # constant: set by output of curator      
      
      self.DataDict={}
      for i in self.curator_data_names: 
          self.DataDict[i]=self.curation_dict[i][max_num_of_added_vertices+1:] # remove manually added vertices + background code (first 1001? places)
      for j in self.feature_names:
          self.DataDict[j]=self._feature_dict[j][chosen_vertex_indices] # remove indices that did not make it into the curator because of some initial volume exclusion algorithm
      if not self.is_training_base_model:
          # self.DataDict[curator_label_names[0]].pop
          # self.DataDict[curator_label_names[0]]=self.foregrnd_dict[curator_label_names[0]] 
          self.DataDict.update({curator_label_names[0]:np.multiply(self.foregrnd_dict['foreground_vertices'],
                                                                   self.curation_dict[ 'displayed_vertices'][max_num_of_added_vertices+1:] )}) # overwrite the displayed vertices to exclude those labeled as background by the base network
          # my_dict.update({'key1': 'value1', 'key2': 'value2'})
          # my_dict.update(key1='value1', key2='value2')
          # my_dict = { **my_dict, 'key1': 'value1', 'key2': 'value2'}
      return self.DataDict

#want to pass x_train list, but dont know how to implement my idea. 

#deciding which features to use
  def format_trainingdata(self,testsplit=0.8 ):

        #xtrainlist=vertex_energies,log_radius

        #deciding which features to use. maybe try for loop.

        # featurematrix=np.asarray(list(zip(self.DataDict["laplacian"], self.DataDict["log_radius"], self.DataDict["gradient_1"], self.DataDict["gradient_2"], self.DataDict["displayed_vertices"],self.DataDict["true_vertices"]))).squeeze()
        # featurematrix  = 
        # training_data_names = self.feature_names + ["displayed_vertices","true_vertices"]
        training_data_names = self.feature_names + self.curator_data_names[-2:]

        featurematrix = np.zeros((len(self.DataDict[training_data_names[0]]),
                                  len(              training_data_names    ) ))

        feature_idx = -1

        for data_name in training_data_names:
          feature_idx += 1
          featurematrix[:,feature_idx]=np.asarray(self.DataDict[data_name]).squeeze()

        # featurematrix = np.array(np.asarray(self.DataDict[training_data_names[0]])) # create empty array to begin concatenating onto (if this doesn't work, just initialize the array to the first feature) # SAM 10/1/22
        # for data_name in training_data_names[1:]: # concatenate feature columns into an array, soft-coding what used to be hard-coded # SAM 10/1/22
        #     featurematrix=np.concatenate((featurematrix,np.asarray(self.DataDict[data_name])),axis=1)

        # print("Length before deleting",len(featurematrix))

        #filter based on 3rd parameter or coloumn, so displayed vertices. 
        # ????????? should be using the second from last feature (-2 not 2 index) to call out the displayed vertices ??????????????? SAM 10/1/22
        featurematrix = np.delete(featurematrix, np.where((featurematrix[:, -2] == 0))[0], axis=0)
        # featurematrix = np.delete(featurematrix, np.where((featurematrix[:, 2] == 0))[0], axis=0)
        # print("Length after deleting",len(featurematrix))

        # ??????? these lines appeared to cancel each other out - lines commented out - oh I see now, they are deleting the displayed vertices column: replaced by one line ?????????? SAM 10/1/22
        featurematrix = np.delete(featurematrix, -2, axis=1)
        # cleaned_laplacians=np.expand_dims(featurematrix[:,0], 0).T
        # cleaned_log_radius=np.expand_dims(featurematrix[:,1], 0).T
        # cleaned_gradient_1=np.expand_dims(featurematrix[:,2], 0).T
        # cleaned_gradient_2=np.expand_dims(featurematrix[:,3], 0).T
        # cleaned_true_vertices=np.expand_dims(featurematrix[:,5], 0).T


        # featurematrix=np.concatenate((cleaned_laplacians,cleaned_log_radius, cleaned_gradient_1, cleaned_gradient_2, cleaned_true_vertices), axis=1)
        np.random.shuffle(featurematrix) 
        # from: https://numpy.org/doc/stable/reference/random/generated/numpy.random.shuffle.html
        # 
        # Multi-dimensional arrays are only shuffled along the first axis:
        # 
        # arr = np.arange(9).reshape((3, 3))
        # np.random.shuffle(arr)
        # arr
        # array([[3, 4, 5], # random
        #       [6, 7, 8],
        #       [0, 1, 2]])

# possibly change 2 to -1. 

        x=featurematrix[:,:-1]
        y=featurematrix[:, -1]

        #now to split data into training and test sets.

        thesplitvalue=int(testsplit*len(x))

        xtrain=x[:thesplitvalue]
        ytrain=y[:thesplitvalue]

        xtest=x[thesplitvalue:]
        ytest=y[thesplitvalue:]

        return xtrain, ytrain, xtest, ytest

def GenerateMachineCurator(feature_names,training_set_name,*args):
  
  xtrain = np.zeros((0,len(feature_names)))
  x_test = np.zeros((0,len(feature_names)))
  ytrain = np.zeros((0 ))
  y_test = np.zeros((0 ))
  
  for curation_file in args:
    VectorDataset1=TrainingData(curation_file)
    # VectorDataset1=TrainingData(curation_file)
    # # VectorDataset1=mll.TrainingData(curation_file, feature__file)
    # # # VectorDataset1=TrainingData(curation_file, feature_file, foreground_file)
    # # #VectorDataset1=TrainingData("/content/vertices_220518-185057_Doug_Fused_Raw_w5.mat","/content/vertices_220518-185057_features_Doug_Fused_Raw_w5.mat")
    # # #VectorDataset1=TrainingData("/content/vertices_220126-151318_Doug_Fused_Raw_w1.mat","/content/vertices_220126-151318_features_Doug_Fused_Raw_w1.mat")
    # # #VectorDataset1=TrainingData("/content/vertices_220511-223915_Doug_Fused_Raw_w3.mat","/content/vertices_220511-223915_features_Doug_Fused_Raw_w3.mat")

    # VectorDataset1.get_features()
    # thedata=VectorDataset1.ChooseFeatures(feature_names)
    VectorDataset1.ChooseFeatures(feature_names)
    xtrain_temp,ytrain_temp,x_test_temp,y_test_temp=VectorDataset1.format_trainingdata(0.7)

    xtrain = np.concatenate(( xtrain, xtrain_temp ), axis = 0 )
    ytrain = np.concatenate(( ytrain, ytrain_temp ), axis = 0 )
    x_test = np.concatenate(( x_test, x_test_temp ), axis = 0 )
    y_test = np.concatenate(( y_test, y_test_temp ), axis = 0 )


  mymodel=WorkingModel(training_set_name,feature_names,0.01)

  # mymodel.name = feature_names[0]

  # mymodel.feature_names = feature_names

  thehistory=train_model(mymodel,xtrain,ytrain,5,10)

  model_directory = "MLModels"

  if not os.path.isdir(model_directory):
    os.mkdir(model_directory)
    
  generate_plots(thehistory,model_directory,mymodel,x_test,y_test)

  #instead of the self.name
  #trained

  model_file = model_directory + os.sep + mymodel.name+".h5"

  mymodel.save( model_file )

  # VectorDataset2 = mll.DeploymentData( mymodel, mymodel.feature_names )
  # # VectorDataset1.predictonfeatures(mymodel,feature__file) # !!!!! ****** mymodel.Predict(model_file,feature_file)
  # # # VectorDataset1.predictonfeatures(mymodel,"/content/vertices_220126-151318_features_Doug_Fused_Raw_w1.mat")
  # # # # VectorDataset1.predictonfeatures(mymodel,"/content/vertices_220511-223915_features_Doug_Fused_Raw_w3.mat")
  # # # # VectorDataset1.predictonfeatures(mymodel,"/content/vertices_220518-185057_features_Doug_Fused_Raw_w5.mat")

  # mymodel.save(mymodel.name+".h5",feature__file)
  # mymodel.save("mymodel.h5","/content/vertices_220126-151318_features_Doug_Fused_Raw_w1.mat")

  #self,model,PredictFeatureFileAddress

def SingleNeuronModel(input_feature_length,learning_rate,loss_function="binary_crossentropy"):
    model = keras.Sequential([
    layers.Dense(units=1, input_shape=[input_feature_length],activation="sigmoid")])
    model.compile(loss=loss_function, optimizer=tf.keras.optimizers.Adam(learning_rate), metrics=['accuracy'])
    model.summary() 
    return model

# def WorkingModel(feature_names,learning_rate,loss_function="binary_crossentropy"):
#       input_feature_length = len(feature_names)
      
#       # model_name = feature_names[0]
#       # for name in feature_names[1:]:
#       #   model_name += '-' + name
#       model_name = '-'.join(feature_names)

#       # Multilayer Perceptron
#       from keras.models import Model
#       from keras.layers import Input
#       from keras.layers import Dense
#       visible = Input(shape=(input_feature_length,))
#       hidden1 = Dense(10, activation='relu'   )(visible)
#       hidden2 = Dense(20, activation='relu'   )(hidden1)
#       hidden3 = Dense(10, activation='relu'   )(hidden2)
#       output  = Dense(1,  activation='sigmoid')(hidden3)
#       model   = Model(inputs=visible, outputs=output, name=model_name)
#       model.compile(loss=loss_function, optimizer=tf.keras.optimizers.Adam(learning_rate), metrics=['accuracy'])
#       model.summary()

#       return model

def WorkingModel(training_set_name,feature_names,learning_rate,loss_function="binary_crossentropy"):
      input_feature_length = len(feature_names)
      
      # # model_name = feature_names[0]
      # # for name in feature_names[1:]:
      # #   model_name += '-' + name
      # model_name = 'MLModel_'+training_set_name+'-'+'-'.join(feature_names)
      model_name = training_set_name+'-'+'-'.join(feature_names)

      # Multilayer Perceptron
      from keras.models import Model
      from keras.layers import Input
      from keras.layers import Dense
      visible = Input(shape=(input_feature_length,))
      # hidden1 = Dense(5, activation='relu'   )(visible)
      # output  = Dense(1,  activation='sigmoid')(hidden1)
      output  = Dense(1,  activation='sigmoid')(visible)
      model   = Model(inputs=visible, outputs=output, name=model_name)
      model.compile(loss=loss_function, optimizer=tf.keras.optimizers.Adam(learning_rate), metrics=['accuracy'])
      model.summary()

      return model

def train_model(model,xtrain,ytrain,epochs,batch_size):
    from sklearn.utils.class_weight import compute_class_weight
    class_weights = compute_class_weight(class_weight = "balanced", classes= np.unique(ytrain), y= ytrain)
    weight = {i : class_weights[i] for i in range(2)}
    history=model.fit(xtrain,ytrain,epochs=epochs, batch_size=batch_size, class_weight=weight)

    return history

def cutoff_youdens_j(fpr,tpr,thresholds):
    j_scores = (tpr+ 1 -fpr)/2
    j_ordered = sorted(zip(j_scores,thresholds))
    return j_ordered[-1][:]

def generate_plots(history,model_directory,model,xtest,ytest):
  inputtitle= model_directory + os.sep + model.name+"_"
  loss=history.history["loss"]
  accuracy=history.history["accuracy"]
  epochs=range(1,len(loss)+1)
  #plt.plot(epochs,loss,"y",label="loss values")
  #plt.savefig("loss")

  #plt.plot(epochs,accuracy,"r",label="accuracy values")
  #plt.savefig("accuracy")

  fig, ax = plt.subplots()
  ax.plot(epochs, loss)
  ax.set_title('loss')
  plt.show()
  plt.savefig(inputtitle+"loss_plot")

  fig, ax = plt.subplots()
  ax.plot(epochs, accuracy)
  ax.set_title('accuracy')
  plt.show()
  plt.savefig(inputtitle+"accuracy_plot")

  labels=ytest
  predictions=model.predict(xtest)

  from sklearn.metrics import roc_curve
  from sklearn.metrics import auc

  fpr_rf, tpr_rf, thresholds_rf = roc_curve(list(labels), list(predictions))
  auc_rf = auc(fpr_rf, tpr_rf)

  plt.figure(20)
  plt.plot([0, 1], [0, 1], 'k--')
  plt.plot(fpr_rf, tpr_rf, label='Model (area = {:.3f})'.format(auc_rf))
  plt.xlabel('False positive rate')
  plt.ylabel('True positive rate')
  plt.title('ROC curve')
  plt.legend(loc='best')
  plt.show()
  plt.savefig(inputtitle+"ROC_curve")
  # Zoom in view of the upper left corner.
  plt.figure(2)
  plt.xlim(0, 0.2)
  plt.ylim(0.8, 1)
  plt.plot([0, 1], [0, 1], 'k--')
  plt.plot(fpr_rf, tpr_rf, label='Model (area = {:.3f})'.format(auc_rf))
  plt.xlabel('False positive rate')
  plt.ylabel('True positive rate')
  plt.title('ROC curve (zoomed in at top left)')
  plt.legend(loc='best')
  plt.show()
  plt.savefig(inputtitle+"top_lefthand_ROC_curve")

  print("Threshold at max. balanced accuracy: ",cutoff_youdens_j(fpr_rf, tpr_rf, thresholds_rf)[1])
  print(             "Max. Balanced Accuracy: ",cutoff_youdens_j(fpr_rf, tpr_rf, thresholds_rf)[0])