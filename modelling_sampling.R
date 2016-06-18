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

set.seed(737)
pop = as.matrix(data)

comp_params <- function(x, init_params, n_iter = 100, run_all){
  
  ### output: 
  ### x: the observations for which responsibilties will be calculated
  ### init_params: the initial values for the beta distributions
  ###             1 = alpha parameter for first dist
  ###             2 = beta parameter for first dist
  ###             3 = alpha parameter for second dist
  ###             4 = beta parameter for second dist
  ###             5 = the class split
  ### n_iter: the number of times the estimation process will be run to generate results. Default value = 100
  ### run_all: 
  ###         - TRUE: Runs the much larger sampling method for estimating the beta parameters
  ###         - FALSE: Fits a simple beta distribution to the median class memberships from the bernoulli draws. Good for getting a baseline
  
  alpha1 <- init_params[1]
  beta1 <- init_params[2]
  alpha2 <- init_params[3]
  beta2 <- init_params[4]
  pi <- init_params[5]
  
  ### Generating probabilities that an observation would come from a given class
  
  initial_pop1 <- dbeta(x, alpha1, beta1)
  initial_pop2 <- dbeta(x, alpha2, beta2)
  
  resp <- matrix(ncol = length(x))
  
  ### Normalising class responsibilities to ensure that resp + (1 - resp) = 1
  
  for (i in 1:length(x)){
   
   resp[i] = (pi * initial_pop2[i]) / ((1 - pi) * initial_pop1[i] + pi * initial_pop2[i])
   
  }

  results <- matrix(nrow = n_iter, ncol = 5)
  class_matrix <- matrix(nrow = n_iter, ncol = length(pop))
  
  for (i in 1:n_iter){
      
        class_matrix[i, ] <- sapply(1:length(resp), function(y) rbinom(n = 1, size = 1, prob = resp[y]))
        
        if (run_all == TRUE){
        
          class_alloc <- class_matrix[i, ]
          class_alloc[class_alloc == 0.5] <- rbinom(n = 1, size = 1, prob = 0.5)
          df <- data.frame(pop = x, class_type = class_alloc)
           
          pop_0 <- as.numeric(df$pop[df$class_type == 0])
          pop_1 <- as.numeric(df$pop[df$class_type == 1])
           
          dist_0 <- fit_beta(pop_0, alpha1, beta1)
          dist_1 <- fit_beta(pop_1, alpha2, beta2)
           
          results[i, 1] <- dist_0['shape1'] ### alpha parameter estimate for latent class 0
          results[i, 2] <- dist_0['shape2'] ### beta parameter estimate for latent class 0
          results[i, 3] <- dist_1['shape1'] ### alpha parameter estimate for latent class 1
          results[i, 4] <- dist_1['shape2'] ### beta parameter estimate for latent class 1
        } 
  }
  
  if (run_all == FALSE){
  
    class_alloc <- apply(class_matrix, 2, median)
    class_alloc[class_alloc == 0.5] <- rbinom(n = 1, size = 1, prob = 0.5)
     
    ## Dealing with classes that have medians of 0.5 (balanced counts of 0's and 1's). Should be very unlikey outside testing
    
    df <- data.frame(pop = x, class = class_alloc)
    
    pop_0 <- as.numeric(df$pop[df$class == 0])
    pop_1 <- as.numeric(df$pop[df$class == 1])
    
    dist_0 <- fit_beta(pop_0, alpha1, beta1)
    dist_1 <- fit_beta(pop_1, alpha2, beta2)
        
    shape1_0 <- dist_0['shape1'] ### alpha parameter estimate for latent class 0
    shape2_0 <- dist_0['shape2'] ### beta parameter estimate for latent class 0
    shape1_1 <- dist_1['shape1'] ### alpha parameter estimate for latent class 1
    shape2_1 <- dist_1['shape2'] ### beta parameter estimate for latent class 1
    pi_updated <- sum(resp) / length(x) ### Class membership split
      
    output <- c(shape1_0, shape2_0, shape1_1, shape2_1, pi_updated)
  
  }
  
  if (run_all == TRUE){
  
    results[ ,5] <- sum(resp) / length(x) ### Class membership split
    output <- apply(results, 2, median)
  
  }
  print('Exiting program')
  return(output)

}

fit_beta <- function(subpop, shp1, shp2){
  
  ### output: returns estimated parameter estimates for the supplied sub-population
  ### subpop: the number under invesitgation
  ### shp1: the initial alpha estimate for the underlying beta dist
  ### shp2: the initial beta estimate for the underlying beta dist
  ### note: BFGS method was defined as the out of the box estimator kept giving errors. See here for more information:
  ### http://www.inside-r.org/r-doc/stats/optim
  
  
  beta_dist <- fitdistr(subpop, dbeta, list(shape1 = shp1, shape2 = shp2), method = "BFGS")
  dist <- c(beta_dist$estimate['shape1'], beta_dist$estimate['shape2'])
  
  return(dist)
  
}

graph_results <- function(population, estimates){
  
  ### ouptput: a graph comparing the population distribution to the estimated distribution. 
  ### Population: a matrix representing the data that is to be analysed
  ### estimates: results from the mixture model fitted on the population data
  ### Note that as this code does random draws multiple runs will give different estimated curved
  
  library(distr)
  library(reshape2)
  
  myMix <- UnivarMixingDistribution(Beta(shape1=estimates[1], shape2=estimates[2]), 
                                    Beta(shape1=estimates[3], shape2=estimates[4]),
                                    mixCoeff=c(1 - estimates[5], estimates[5]))
  
  rmyMix <- r(myMix)
  x_test <- rmyMix(8000)
  
  
  #### Getting an overview of what the data looks like
  df <- data.frame(pop = population, estimate = x_test)
  
  colnames(df)[1] <- 'pop'
  
  data_adj <- melt(df)
  ggplot(data_adj ,aes(x=value, fill=variable)) + geom_density(alpha=0.25) + theme_minimal() + ggtitle("Plotting actual vs estimated")
}

alpha1 <- 0.5
beta1 <- 0.25
alpha2 <- 0.25
beta2 <- 0.5
pi <- 0.5

initial_estimate <- comp_params(pop, init, 5000, FALSE)
graph_results(pop, initial_estimate)

### First iteration - small number of random draws

refined_est <- comp_params(pop, initial_estimate, 500, TRUE)
graph_results(pop, refined_est, 'est_1.png')

refined_est_1 <- comp_params(pop, refined_est, 500, TRUE)
graph_results(pop, refined_est_1)

refined_est_2 <- comp_params(pop, refined_est_1, 500, TRUE)
graph_results(pop, refined_est_2)

refined_est_3 <- comp_params(pop, refined_est_2, 500, TRUE)
graph_results(pop, refined_est_3)

refined_est_4 <- comp_params(pop, refined_est_3, 500, TRUE)
graph_results(pop, refined_est_4)

refined_est_5 <- comp_params(pop, refined_est_4, 500, TRUE)
graph_results(pop, refined_est_5)

### Trying out a deeper search of the parameter space

initial_estimate <- comp_params(pop, init, 5000, FALSE)
graph_results(pop, initial_estimate)

refined_est_1_1000 <- comp_params(pop, initial_estimate, 1000, TRUE)
graph_results(pop, refined_est_1_1000)

refined_est_2_1000 <- comp_params(pop, refined_est_1_1000, 1000, TRUE)
graph_results(pop, refined_est_2_1000)

refined_est_3_1000 <- comp_params(pop, refined_est_2_1000, 1000, TRUE)
graph_results(pop, refined_est_3_1000)

refined_est_4_1000 <- comp_params(pop, refined_est_3_1000, 1000, TRUE)
graph_results(pop, refined_est_4_1000)

refined_est_5_1000 <- comp_params(pop, refined_est_4_1000, 1000, TRUE)
graph_results(pop, refined_est_5_1000)

refined_est_6_1000 <- comp_params(pop, refined_est_5_1000, 1000, TRUE)
graph_results(pop, refined_est_5_1000)

refined_est_7_1000 <- comp_params(pop, refined_est_6_1000, 1000, TRUE)
graph_results(pop, refined_est_7_1000)

refined_est_8_1000 <- comp_params(pop, refined_est_7_1000, 1000, TRUE)
graph_results(pop, refined_est_8_1000)