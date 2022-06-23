#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:27:05 2022

@author: cal
"""

import os
import pandas as pd
import tensorflow as tf
from keras.layers import Activation
from keras.utils.generic_utils import get_custom_objects

print(tf.__version__)


# ==========================
# Read in data
# ==========================

# potentially helpful link for transfer learning:
    # https://www.tensorflow.org/guide/keras/transfer_learning#do_a_round_of_fine-tuning_of_the_entire_model

os.chdir('C:/Users/ladwi/Documents/Projects/R/LakeModeling')

inputdata=pd.read_csv("1D_HeatMixing_PGDL/output/meteorology_input.csv", parse_dates=['time']) # missing buoyancy?
inputdata['day_of_year'] = inputdata['time'].dt.dayofyear
inputdata['time_of_day'] = inputdata['time'].dt.hour / 24
buoy = pd.read_csv("1D_HeatMixing_PGDL/output/buoyancy.csv", parse_dates=['time'])
inputdata = inputdata.join(buoy.iloc[:, 1:]).drop(["Area_m2", "n2S-2_25"], axis=1) # note: if you don't drop area, things blow up

t_diffvals = pd.read_csv("1D_HeatMixing_PGDL/output/temp_diff01.csv", parse_dates=['time'])
d_diffvals = pd.read_csv("1D_HeatMixing_PGDL/output/diff.csv", parse_dates=['time'])
t_mixvals = pd.read_csv("1D_HeatMixing_PGDL/output/temp_mix02.csv", parse_dates=['time'])
t_convals = pd.read_csv("1D_HeatMixing_PGDL/output/temp_conv03.csv", parse_dates=['time'])
t_icevals = pd.read_csv("1D_HeatMixing_PGDL/output/temp_total04.csv", parse_dates=['time'])

## QUESTIONS
# Normalizing inputs prior to running model?
# How to set change in learning rate in Adam optimizer?
# Which metric is "Relative L2 Error" in tensorflow?
# How is depth incorporated as an input? What are shapes of inputs and outputs?
# How are individual MLPs combined into one for fine tuning?
# How is fine tuning loss calculated?


# inputdata_df = inputdata.iloc[1: , :]
# inputdata_df.append(targetdata[:-1], ignore_index = True)
# inputdata_df.head

# targetdata_df = targetdata.iloc[1: , :]

# ==========================
# Setup input/outputs
# ==========================
m1_inputs = inputdata.iloc[1:, 1:].join(t_icevals.iloc[:-1, 1:])
m1_outputs = t_diffvals.iloc[1:, 1:]
m2_inputs = t_diffvals.iloc[1:, 1:].join(inputdata.iloc[1:,:][["ShearVelocity_mS-1", "ShearStress_Nm-2", "day_of_year", "time_of_day"]])
m2_outputs = t_mixvals.iloc[1:,1:]
m3_inputs = t_mixvals.iloc[1:,1:].join(inputdata.iloc[1:,:][["day_of_year", "time_of_day"]])
m3_outputs = t_convals.iloc[1:,1:]
m4_inputs = t_convals.iloc[1:,1:].join(inputdata.iloc[1:,:][["day_of_year", "time_of_day"]])
m4_outputs = t_icevals.iloc[1:,1:]

# ==========================
# Individual models for each step
# ==========================
train_end = 1276
# m1
m1 = tf.keras.Sequential([
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(25)
])
m1.compile(optimizer="adam",
           loss="mse",
           metrics=["mae", "mse"])
m1.fit(m1_inputs.iloc[:train_end], m1_outputs.iloc[:train_end], epochs=1000)
m1.evaluate(m1_inputs.iloc[:train_end], m1_outputs.iloc[:train_end])
m1.evaluate(m1_inputs.iloc[train_end:-1], m1_outputs.iloc[train_end:-1]) # about 2-3x higher RMSE than Arka

# m2
m2 = tf.keras.Sequential([
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(25)
])
m2.compile(optimizer="adam",
           loss="mse",
           metrics=["mae", "mse"])
m2.fit(m2_inputs.iloc[:train_end], m2_outputs.iloc[:train_end], epochs=1000)
m2.evaluate(m2_inputs.iloc[:train_end], m2_outputs.iloc[:train_end])
m2.evaluate(m2_inputs.iloc[train_end:-1], m2_outputs.iloc[train_end:-1]) # very similar RMSE to Arka

# m3
m3 = tf.keras.Sequential([
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(25)
])
m3.compile(optimizer="adam",
           loss="mse",
           metrics=["mae", "mse"])
m3.fit(m3_inputs.iloc[:train_end], m3_outputs.iloc[:train_end], epochs=1000)
m3.evaluate(m3_inputs.iloc[:train_end], m3_outputs.iloc[:train_end])
m3.evaluate(m3_inputs.iloc[train_end:-1], m3_outputs.iloc[train_end:-1]) # train RMSE ~10x higher, test RMSE ~4x higher than Arka

# m4
m4 = tf.keras.Sequential([
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(32, activation="gelu"),
    tf.keras.layers.Dense(25)
])
m4.compile(optimizer="adam",
           loss="mse",
           metrics=["mae", "mse"])
m4.fit(m4_inputs.iloc[:train_end], m4_outputs.iloc[:train_end], epochs=1000)
m4.evaluate(m4_inputs.iloc[:train_end], m4_outputs.iloc[:train_end])
m4.evaluate(m4_inputs.iloc[train_end:-1], m4_outputs.iloc[train_end:-1]) # train RMSE ~ 50x higher, test RMSE ~3x higher

# ==========================
# combined models sequentially and get outputs
# ==========================

# ==========================
# "stack" models and fine tune together
# ==========================


# ==========================
# fine tune above on observations

# ==========================