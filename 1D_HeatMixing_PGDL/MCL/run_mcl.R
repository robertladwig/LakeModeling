cat("\f")
rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(keras)
library(tidyverse)
library(lubridate)

set.seed(123)

input_df <- read_csv('../output/meteorology_input.csv')
output_df <- read_csv('../output/temp_total04.csv')
buoy <- read_csv('../output/buoyancy.csv')
observed <- read_csv('../output/observed_temp.csv')

buoy <- buoy %>%
  select(-time)
input_df <- cbind(input_df,buoy) %>%
  select(-c(time,Area_m2,'n2S-2_25'))
output_df <- output_df %>%
  select(-time)
observed <- observed %>%
  select(-time)

input <- as.matrix(cbind(input_df[2:nrow(input_df),], output_df[1:(nrow(output_df)-1),]))
target <- as.matrix(output_df[2:nrow(output_df),])
target <- as.matrix(observed[2:nrow(observed),])

# network <- keras_model_sequential() %>%
#   layer_dense(units = 15, activation = "sigmoid", input_shape = c(34)) %>%
#   layer_dense(units = 25)

idx <- floor(nrow(input) * 0.6)
train_features <- input[1:idx,]
train_target <- target[1:idx,]

test_features <- input[(idx+1):nrow(input),]
test_target <- target[(idx+1):nrow(input),]

normalize <- function(x){
  num <- x - min(x)
  denom <- max(x) - min(x)
  return(num/denom)
}

# normalise stuff?
# train_features <- as.matrix(apply(train_features, 2, normalize))
# train_target <- as.matrix(apply(train_target, 2, normalize))
# test_features <- as.matrix(apply(test_features, 2, normalize))
# test_target <- as.matrix(apply(test_target, 2, normalize))

network = keras_model_sequential() %>% 
  layer_dense(units=32, activation="gelu", input_shape=c(57)) %>% 
  # layer_dropout(rate = 0.2) %>% 
  layer_dense(units=32, activation = "gelu") %>%
  # layer_dropout(rate = 0.1) %>% 
  layer_dense(units=25, activation="linear")

summary(network)

network %>% compile(optimizer = 'adam', 
                    loss = "mse", metrics = c("mean_absolute_error"))

network %>% fit(train_features, train_target, epochs = 1000, validation_split = 0.4)

history <- network %>% fit(train_features, train_target, epochs = 1000)
plot(history)

scores = network %>% evaluate(train_features, train_target, verbose = 0)
print(scores)

pred_train <- network %>% predict(train_features, batch_size = 128)

pred_test <- network %>% predict(test_features, batch_size = 128)

val.scores = network %>% evaluate(test_features, test_target, verbose = 0)
print(val.scores)
