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

alpha1 <- 0.5
beta1 <- 0.125
alpha2 <- 0.5
beta2 <- 0.25
pi <- 0.5

init <- c(alpha1, beta1, alpha2, beta2, pi)

initial <- pop[1]

fit_beta <- function(df, shp1, shp2, class){
  
  ### df: the dataframe of results containing:
  ###     obs: the original observations
  ###     resp: the class responsiveness
  ### params: the parameters required for fitting the beta distribution
  
  pop_class <- as.matrix(df$obs[df$class == class])
  beta_dist <- fitdistr(pop_class, dbeta, list(shape1 = shp1, shape2 = shp2))
  
  dist <- c(beta_dist$estimate['shape1'], beta_dist$estimate['shape2'])
  
  return(dist)
  
}

comp_resp <- function(x, init_params, n_iter){
  
  ### x: the observations for which responsibilties will be calculated
  ### init_params: the initial values for the beta distributions
  ###             1 = alpha parameter for first dist
  ###             2 = beta parameter for first dist
  ###             3 = alpha parameter for second dist
  ###             4 = beta parameter for second dist
  ###             5 = the class split
  ### n_iter: the number of times the estimation process will be run to generate results
  
  alpha1 <- init_params[1]
  beta1 <- init_params[2]
  alpha2 <- init_params[3]
  beta2 <- init_params[4]
  pi <- init_params[5]

  resp <- NULL
  initial_pop1 <- NULL
  initial_pop2 <- NULL
  
  initial_pop1 <- sapply(1:length(x), function(y) pbeta(x[y], alpha1, beta1, ncp = 0, lower.tail = TRUE, log.p = FALSE))
  initial_pop2 <- sapply(1:length(x), function(y) pbeta(x[y], alpha2, beta2, ncp = 0, lower.tail = TRUE, log.p = FALSE))
  
  for (i in 1:length(x)){
    
    #initial_pop1[i] <- pbeta(x[i], alpha1, beta1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    #initial_pop2[i] <- pbeta(x[i], alpha2, beta2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      
    resp[i] = (pi * initial_pop2[i]) / ((1 - pi) * initial_pop1[i] + pi * initial_pop2[i])
    
    #### Iterate over these draws 8000-10000 times to get median parameters for next iteration
    #### Keep going until we get convergence
  }
  
  nrow <- NULL
  class_alloc <- NULL
  
  shape1_0 <- NULL
  shape2_0 <- NULL
  
  shape1_1 <- NULL
  shape2_1 <- NULL
  
  for (i in 1:n_iter){
    
    #### Now want to add bernoulli random draws to split into two classes. Split conditional on resp
#    for (j in 1:length(x)){
#      class_alloc[j] <- rbinom(n = 1, size = 1, prob = resp[j])
    
#    }
    
    class_alloc <- sapply(1:length(resp), function(x) rbinom(n = 1, size = 1, prob = resp[x]))
    #### Do I fit the beta distributions based on the observations or the responsiveness measure for the classes?
        
    results <- data.frame(obs = x, resp = resp, class = class_alloc)

    dist_0 <- fit_beta(results, alpha2, beta2, 0)
    dist_1 <- fit_beta(results, alpha1, beta1, 1)
        
    nrow[i] <- i
    
    shape1_0[i] <- dist_0[[1]]
    shape2_0[i] <- dist_0[[2]]
    
    shape1_1[i] <- dist_1[[1]]
    shape2_1[i] <- dist_1[[2]]
      
    
  }
  
  shape_total_draws <- data.frame(row = nrow, class0_shape1 = shape1_0, class0_shape2 = shape2_0,
                                  class1_shape1 = shape1_1, class1_shape2 = shape2_1
                                  )
  
  return(shape_total_draws)
  
}

init_dists <- comp_resp(pop, init, 100)
init_dists[[1]]

summary(as.factor(init_dists$class))

summary(df)