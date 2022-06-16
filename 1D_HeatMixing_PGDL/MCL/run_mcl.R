cat("\f")
rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(keras)
library(tidyverse)
library(lubridate)

set.seed(123)

input_df <- read_csv('../output/meteorology_input.csv')
output_df <- read_csv('../output/temp_total04.csv')

input_df <- input_df %>%
  select(-time)
output_df <- output_df %>%
  select(-time)

input <- as.matrix(cbind(input_df[2:nrow(input_df),], output_df[1:(nrow(output_df)-1),]))
target <- as.matrix(output_df[2:nrow(output_df),])

# network <- keras_model_sequential() %>%
#   layer_dense(units = 15, activation = "sigmoid", input_shape = c(34)) %>%
#   layer_dense(units = 25)


model = keras_model_sequential() %>% 
  layer_dense(units=200, activation="relu", input_shape=c(34)) %>% 
  layer_dropout(rate = 0.2) %>% 
  layer_dense(units=32, activation = "relu") %>%
  layer_dropout(rate = 0.1) %>% 
  layer_dense(units=1, activation="linear")

summary(network)

network %>% compile(optimizer = 'adam', 
                    loss = "mse", metrics = c("mean_absolute_error"))

idx <- floor(nrow(input) * 0.6)
train_features <- input[1:idx,]
train_target <- target[1:idx,]

test_features <- input[(idx+1):nrow(input),]
test_target <- target[(idx+1):nrow(input),]

network %>% fit(train_features, train_target, epochs = 100, validation_split = 0.4)

history <- network %>% fit(train_features, train_target, epochs = 100)
plot(history)

scores = network %>% evaluate(train_features, train_target, verbose = 0)
print(scores)

pred_train <- network %>% predict(train_features, batch_size = 128)

pred_test <- network %>% predict(test_features, batch_size = 128)

# Confusion matrix
table(test_target, classes)
