# https://www.allaboutcircuits.com/technical-articles/how-to-create-a-multilayer-perceptron-neural-network-in-python/
import pandas
import numpy as np
import pandas as pd
import os
import math

def logistic(x):
    return 1.0/(1 + np.exp(-x))

def logistic_deriv(x):
    return logistic(x) * (1 - logistic(x))

# learning rate
LR = 1   

# dimensionality of input layer
I_dim = 9
# dimensionality of hidden layer
H_dim = 4

# epoch count
epoch_count = 1

#np.random.seed(1)
weights_ItoH = np.random.uniform(-1, 1, (I_dim, H_dim))
weights_HtoO = np.random.uniform(-1, 1, H_dim)

preActivation_H = np.zeros(H_dim)
postActivation_H = np.zeros(H_dim)

#training_data = pandas.read_excel('MLP_Tdata.xlsx')
#target_output = training_data.output
#training_data = training_data.drop(['output'], axis=1)
#training_data = np.asarray(training_data)
#training_count = len(training_data[:,0])

#validation_data = pandas.read_excel('MLP_Vdata.xlsx')
#validation_output = validation_data.output
#validation_data = validation_data.drop(['output'], axis=1)
#validation_data = np.asarray(validation_data)
#validation_count = len(validation_data[:,0])


os.chdir('/home/robert/Projects/LakeModeling')

inputdata=pd.read_csv("1D_HeatMixing_PGDL/output/meteorology_input.csv")
print(inputdata)

targetdata = pd.read_csv("1D_HeatMixing_PGDL/output/temp_diff01.csv")
print(targetdata)

inputdata = inputdata.drop(['time'], axis=1)
targetdata = targetdata.drop(['time'], axis=1)


inputdata_df = inputdata.iloc[1: , :]

inputdata_df.append(targetdata[:-1], ignore_index = True)
inputdata_df.head

targetdata_df = targetdata.iloc[1: , :]

# 70:30 train to test
train_len = math.floor(len(inputdata_df) * 0.7)

training_data = np.asarray(inputdata_df.iloc[0:train_len,:])
target_output = np.asarray(targetdata_df.iloc[0:train_len,:])
training_count = len(training_data[:,0])

validation_data = np.asarray(inputdata_df.iloc[(1+train_len):,:])
validation_output = np.asarray(targetdata_df.iloc[(1+train_len):,:])
validation_count = len(validation_data[:,0])


#####################
#training
#####################
for epoch in range(epoch_count):
    for sample in range(training_count):
        for node in range(H_dim):
            preActivation_H[node] = np.dot(training_data[sample,:], weights_ItoH[:, node])
            postActivation_H[node] = logistic(preActivation_H[node])
            
        preActivation_O = np.dot(postActivation_H, weights_HtoO)
        postActivation_O = logistic(preActivation_O)
        
        FE = postActivation_O - target_output[sample]
        
        for H_node in range(H_dim):
            S_error = FE * logistic_deriv(preActivation_O)
            gradient_HtoO = S_error * postActivation_H[H_node]
                       
            for I_node in range(I_dim):
                input_value = training_data[sample, I_node]
                gradient_ItoH = S_error * weights_HtoO[H_node] * logistic_deriv(preActivation_H[H_node]) * input_value
                
                weights_ItoH[I_node, H_node] -= LR * gradient_ItoH
                
            weights_HtoO[H_node] -= LR * gradient_HtoO

#####################
#validation
#####################            
correct_classification_count = 0
for sample in range(validation_count):
    for node in range(H_dim):
        preActivation_H[node] = np.dot(validation_data[sample,:], weights_ItoH[:, node])
        postActivation_H[node] = logistic(preActivation_H[node])
            
    preActivation_O = np.dot(postActivation_H, weights_HtoO)
    postActivation_O = logistic(preActivation_O)
        
    if postActivation_O > 0.5:
        output = 1
    else:
        output = 0     
        
    if output == validation_output[sample]:
        correct_classification_count += 1

print('Percentage of correct classifications:')
print(correct_classification_count*100/validation_count)
