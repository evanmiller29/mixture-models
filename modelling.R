library(MASS)
library(ggplot2)
library(dplyr)

#### Reading in the data

op <- options(digits = 7)

rm(list = ls())

basepath <- 'C:/Users/evanm_000/Documents/GitHub/mixture-models'
setwd(basepath)

data <- read.table("EvanExercise.dat", header=FALSE)
colnames(data) <- 'population'

#### Getting an overview of what the data looks like

ggplot(data, aes(x = population)) + geom_density() + xlim(0, 1.1) +theme_minimal() + ggtitle('Overview of population')
colnames(data) <- NULL

set.seed(737)
pop = as.matrix(data)

#### Starting points are as below. Will use K means to generate more reasonable estimates in future
#### Reducing the size of the starting points to make them more feasible

#### Using k means to select a starting point for the algo:

#kmeans <- kmeans(data, 2, iter.max = 1000, nstart = 1)

alpha1 <- 0.5
beta1 <- 0.25
alpha2 <- 0.25
beta2 <- 0.5
pi <- 0.5

init <- c(alpha1, beta1, alpha2, beta2, pi)

fit_beta <- function(subpop, shp1, shp2, scale){
  
  ### subpop: the number under invesitgation
  ### shp1: the initial alpha estimate for the underlying beta dist
  ### shp2: the initial beta estimate for the underlying beta dist
  
  beta_dist <- fitdistr(subpop, dbeta, list(shape1 = shp1 / scale, shape2 = shp2 / scale), lower=0.01)
  dist <- c(beta_dist$estimate['shape1'], beta_dist$estimate['shape2'])
  
  return(dist)
  
}

comp_params <- function(x, init_params, n_iter = 100, scale){
  
  ### x: the observations for which responsibilties will be calculated
  ### init_params: the initial values for the beta distributions
  ###             1 = alpha parameter for first dist
  ###             2 = beta parameter for first dist
  ###             3 = alpha parameter for second dist
  ###             4 = beta parameter for second dist
  ###             5 = the class split
  ### n_iter: the number of times the estimation process will be run to generate results. Default value = 100
  ### scale: used to reduce the initial estimates of alpha and beta to be reasonable starting points
  
  alpha1 <- init_params[1]
  beta1 <- init_params[2]
  alpha2 <- init_params[3]
  beta2 <- init_params[4]
  pi <- init_params[5]
  
  ### Generating probabilities that an observation would come from a given class
  
  initial_pop1 <- dbeta(x, alpha1, beta1)
  initial_pop2 <- dbeta(x, alpha2, beta2)
  
  resp <- matrix(ncol = length(x))
  
  ### Normalising class responsibilities to ensure that pi + (1 - pi) = 1
  
  for (i in 1:length(x)){
   
   resp[i] = (pi * initial_pop2[i]) / ((1 - pi) * initial_pop1[i] + pi * initial_pop2[i])
   
  }
  
  ### Using randomly allocating the number of 

  results <- matrix(nrow = n_iter, ncol = 5)
  class_matrix <- matrix(nrow = n_iter, ncol = length(pop))
  df <- NULL
  
  for (i in 1:n_iter){
      
        class_matrix[i, ] <- sapply(1:length(resp), function(y) rbinom(n = 1, size = 1, prob = resp[y]))
        # Add in drawing curves etc
        class_alloc <- class_matrix[i, ]
        class_alloc[class_alloc == 0.5] <- rbinom(n = 1, size = 1, prob = 0.5)
        df <- data.frame(pop = x, class_type = class_alloc)
         
        pop_0 <- as.numeric(df$pop[df$class_type == 0])
        pop_1 <- as.numeric(df$pop[df$class_type == 1])
         
        dist_0 <- fit_beta(pop_0, alpha1, beta1, scale)
        dist_1 <- fit_beta(pop_1, alpha2, beta2, scale)
         
         results[i, 1] <- dist_0['shape1'] ### alpha parameter estimate for latent class 0
         results[i, 2] <- dist_0['shape2'] ### beta parameter estimate for latent class 0
         results[i, 3] <- dist_1['shape1'] ### alpha parameter estimate for latent class 1
         results[i, 4] <- dist_1['shape2'] ### beta parameter estimate for latent class 1
         results[i, 5] <- sum(resp) / length(x) ### Class membership split
  }
  
#   class_alloc <- apply(class_matrix, 2, median)
#  class_alloc[class_alloc == 0.5] <- rbinom(n = 1, size = 1, prob = 0.5)
   
  ### Dealing with classes that have medians of 0.5 (balanced counts of 0's and 1's). Should be very unlikey outside testing
  
#   df <- data.frame(pop = x, class = class_alloc)
#   
#   pop_0 <- as.numeric(df$pop[df$class == 0])
#   pop_1 <- as.numeric(df$pop[df$class == 1])
#   
#   dist_0 <- fit_beta(pop_0, alpha1, beta1, scale)
#   dist_1 <- fit_beta(pop_1, alpha2, beta2, scale)
#       
#   shape1_0 <- dist_0['shape1'] ### alpha parameter estimate for latent class 0
#   shape2_0 <- dist_0['shape2'] ### beta parameter estimate for latent class 0
#   shape1_1 <- dist_1['shape1'] ### alpha parameter estimate for latent class 1
#   shape2_1 <- dist_1['shape2'] ### beta parameter estimate for latent class 1
#   pi_updated <- sum(resp) / length(x) ### Class membership split
#     
#   output <- c(shape1_0, shape2_0, shape1_1, shape2_1, pi_updated)
   
  print('Exiting program')

  output <- apply(results, 2, median)
  return(output)

}

library(distr)
 
myMix <- UnivarMixingDistribution(Beta(shape1=20, shape2=0.5), 
                                   Beta(shape1=1.5, shape2=20),
                                   mixCoeff=c(0.6, 0.4))

rmyMix <- r(myMix)
x_test <- rmyMix(1000)

alpha1 <- 0.5
beta1 <- 0.25
alpha2 <- 0.25
beta2 <- 0.5
pi <- 0.5

init <- c(alpha1, beta1, alpha2, beta2, pi)

est1 <- comp_params(x_test, init, 150, 1)
est2 <- comp_params(x_test, est1, 1000, 1)
est3 <- comp_params(x_test, est2, 1000, 1)
est4 <- comp_params(x_test, est3, 1000, 1)

alpha1 <- 0.5
beta1 <- 0.25
alpha2 <- 0.25
beta2 <- 0.5
pi <- 0.5

init <- c(alpha1, beta1, alpha2, beta2, pi)

myMix <- UnivarMixingDistribution(Beta(shape1=10, shape2=4), 
                                  Beta(shape1=7, shape2=10),
                                  mixCoeff=c(0.8, 0.2))

rmyMix <- r(myMix)
x_test <- rmyMix(1000)

est1 <- comp_params(x_test, init, 500, 10)
est2 <- comp_params(x_test, est1, 5000, 10)
est3 <- comp_params(x_test, est2, 8000, 10)
est4 <- comp_params(x_test, est3, 8000, 10)
est5 <- comp_params(x_test, est4, 8000, 10)

