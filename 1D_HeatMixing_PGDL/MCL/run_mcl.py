#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:38:51 2022

@author: robert
"""
from sklearn import datasets
from sklearn import metrics
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pandas as pd
import os

import keras.backend as K
from keras.layers import Input, Dense
from keras.models import Model
from keras.losses import mse
import numpy as np

# Some random training data
features = np.random.rand(100,20)
labels_1 = np.random.rand(100,4)
labels_2 = np.random.rand(100,1)

# Input layer, one hidden layer
input_layer = Input((20,))
dense_1 = Dense(128)(input_layer)

# Two outputs
output_1 = Dense(4)(dense_1)
output_2 = Dense(1)(dense_1)

# Two additional 'inputs' for the labels
label_layer_1 = Input((4,))
label_layer_2 = Input((1,))

# Instantiate model, pass label layers as inputs
model = Model(inputs=[input_layer, label_layer_1, label_layer_2], outputs=[output_1, output_2])

# Construct your custom loss as a tensor
loss = K.mean(mse(output_1, label_layer_1) * mse(output_2, label_layer_2))

# Add loss to model
model.add_loss(loss)

# Compile without specifying a loss
model.compile(optimizer='sgd')

dummy = np.zeros((100,))
model.fit([features, labels_1, labels_2], dummy, epochs=2)


os.chdir('/home/robert/Projects/LakeModeling')


inputdata=pd.read_csv("1D_HeatMixing_PGDL/output/meteorology_input.csv")
print(inputdata)

targetdata = pd.read_csv("1D_HeatMixing_PGDL/output/temp_diff01.csv")
print(targetdata)

inputdata_df = inputdata.iloc[1: , :]
inputdata_df.append(targetdata[:-1], ignore_index = True)
inputdata_df.head

targetdata_df = targetdata.iloc[1: , :]

#X = dataset.data; y = dataset.target
#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30)