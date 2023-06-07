# -*- coding: utf-8 -*-
"""
Created on Wed May 24 19:41:10 2023

@author: Omisore
"""
import numpy as np
from scipy import io
import pandas as pd
from sklearn.metrics import classification_report, confusion_matrix 
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from tensorflow.keras.utils import to_categorical
from sklearn.model_selection import cross_validate
from sklearn.model_selection import KFold
from losses import figPlot

log_dir = "logs/"

#Load the preprocessed data from the matlab file
C = io.loadmat('warped_signals-6_Operators-17_15_trials.mat')
MM_Data = C['SpectSegMusc']

# extract dataset
Warped_Instance_Data = MM_Data[:,1]
Warped_Labels = MM_Data[:,2]
Warped_Labels = np.stack(MM_Data[:,2])

#Set for network
Warped_Instance = np.stack(Warped_Instance_Data, axis=0)
Warped_Instance = Warped_Instance.reshape(Warped_Instance.shape[0], Warped_Instance.shape[2])
#Warped_Labels = to_categorical(Warped_Labels-1)
Warped_Labels = Warped_Labels.reshape(Warped_Labels.shape[0], Warped_Labels.shape[2])



ModelName= 'BioCAS-MLP.v.1'

n_folds = 5
epochs = 20
batch_size = 128
model_history = []

train_x, test_x, train_y, test_y = train_test_split(Warped_Instance, Warped_Labels, 
                                                    random_state = 101, test_size=0.25)

kf = KFold(n_splits=10)
model = MLPClassifier(solver='adam', hidden_layer_sizes=(20, 5, 2),
                      batch_size=128, activation='relu', alpha=1e-5,
                      learning_rate='adaptive', learning_rate_init=1e-3,
                      max_iter=400, shuffle=True, random_state=1,
                      tol=1e-4, early_stopping=True,)

for train_indices, test_indices in kf.split(train_x):
    model.fit(train_x[train_indices], train_y[train_indices].ravel())
    print(model.score(train_x[test_indices], train_y[test_indices].ravel()))

model_prediction = model.predict(test_x)

figPlot(test_y.ravel(), model_prediction, ax_fz=30, fz=30, saveAs = 'Output-MLP.png', title='MLP')

#OR-2
#model = MLPClassifier()
#cv_results = cross_validate(model, Warped_Instance, Warped_Labels, cv=10,return_train_score=False, scoring=model.score) 
#print("Fit scores: {}".format(cv_results['test_score']))




"""
plt.title('MLP: Accuracies vs Epochs')
for i in range(n_folds):
    plt.plot(model_history[i].history['acc'], label='Training Fold '+str(i))
    plt.plot(model_history[i].history['val_acc'], label='Training Fold '+str(i))
plt.legend()
plt.show()
"""